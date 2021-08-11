#include "equations.hpp"
#include "math.hpp"
#include "micro/demag_field.hpp"
#include "micro/exchange_field.hpp"
#include "micro/external_field.hpp"
#include "util/arg_parser.hpp"
#include <boost/numeric/odeint.hpp>
#include <iostream>
#include <utility>

// Necessary template specialization for af::array
namespace boost::numeric::odeint {
template <> struct vector_space_norm_inf<af::array> {
    typedef double result_type; // typedef definition is explicitly needed here for odeint internals.
    result_type operator()(const af::array& x) const {
        return af::max(af::max(af::max(af::max(af::abs(x), 0), 1), 2), 3).as(f64).scalar<double>();
    }
};
} // namespace boost::numeric::odeint

using namespace magnumaf;

// // C++20 concept
// template <typename T> concept HasHfield = requires(T field_term, const State& state) { field_term.H_in_Apm(state); };

// Type Erasure:
class HfieldType {
  public:
    // ctor: copy-then-move idiom
    template <typename T> HfieldType(T term) : pimpl_(std::make_unique<Model<T>>(std::move(term))) {}

    ~HfieldType() = default;
    HfieldType(HfieldType const& s) : pimpl_(s.pimpl_->clone()) {}
    HfieldType(HfieldType&& s) = default;
    HfieldType& operator=(const HfieldType& s) {
        HfieldType tmp(s);
        std::swap(pimpl_, tmp.pimpl_);
        return *this;
    }
    HfieldType& operator=(HfieldType&& s) = default;

  private:
    friend af::array Expose_H_in_Apm(HfieldType const& term, State const& state) {
        return term.pimpl_->wrap_H_in_Apm(state);
    }
    friend double Expose_Energy_in_J(HfieldType const& term, State const& state) {
        return term.pimpl_->wrap_Energy_in_J(state);
    }

    struct Concept {
        virtual ~Concept() = default;
        [[nodiscard]] virtual af::array wrap_H_in_Apm(State const&) const = 0;
        [[nodiscard]] virtual double wrap_Energy_in_J(State const&) const = 0;
        [[nodiscard]] virtual std::unique_ptr<Concept> clone() const = 0; // Prototype design pattern
    };

    // ctor: perfect-forwarding
    // enable with c++20 // template <HasHfield T>
    template <typename T> struct Model final : public Concept {
        explicit Model(T&& term) : term_(std::forward<T>(term)) {}
        [[nodiscard]] af::array wrap_H_in_Apm(State const& state) const override { return term_.H_in_Apm(state); }
        [[nodiscard]] double wrap_Energy_in_J(State const& state) const override { return term_.Energy_in_J(state); }
        [[nodiscard]] std::unique_ptr<Concept> clone() const override { return std::make_unique<Model>(*this); }
        T term_;
    };

    std::unique_ptr<Concept> pimpl_;
};

using HfieldTypes = std::vector<HfieldType>;

/// Calculate effective field by accumulating all h(state) terms in container.
/// Expects non-empty container with at least one element
/// Value semantics:
template <typename T> af::array value_Heff_in_Apm(const T& fieldterms, const State& state) {
    const auto H_init = Expose_H_in_Apm(fieldterms[0], state);
    const auto accum = [&state](const auto& sum, const auto& elem) { return sum + Expose_H_in_Apm(elem, state); };
    return std::accumulate(std::begin(fieldterms) + 1, std::end(fieldterms), H_init, accum);
}

int main(int argc, char** argv) {
    std::cout << "Start" << std::endl;

    const auto [outdir_tmp, posargs] = ArgParser(argc, argv).outdir_posargs;
    const auto outdir = outdir_tmp;

    constexpr double x_magaf = 5.e-7, y = 1.25e-7, z = 3.e-9;
    constexpr int nx = 100, ny = 25, nz = 1;
    constexpr Mesh mesh{nx, ny, nz, x_magaf / nx, y / ny, z / nz};

    double alpha = 1.;
    constexpr double Ms = 8e5;
    constexpr double A = 1.3e-11;

    auto fieldterms = HfieldTypes{DemagField(mesh, true, true, 0), ExchangeField(A)};

    constexpr auto type = af::dtype::f64;
    af::array m = af::constant(0, mesh::dims_v(mesh), type);
    m(0, af::span, af::span, 1) = 1;
    m(af::seq(1, af::end - 1), af::span, af::span, 0) = 1;
    m(-1, af::span, af::span, 1) = 1;

    State state(mesh, Ms, m);
    state.write_vti(outdir / "minit");

    for (const auto& i : fieldterms) {
        std::cout << "E=" << Expose_Energy_in_J(i, state) << '\n';
    }

    std::ofstream stream(outdir / "m.dat"); // clears old file
    stream.precision(12);
    std::vector<std::size_t> steps;
    namespace ode = boost::numeric::odeint;

    using stepper_type = ode::runge_kutta_dopri5<af::array, double, af::array, double, ode::vector_space_algebra>;
    constexpr double eps_abs = 1e-6;
    constexpr double eps_rel = 1e-6;
    const auto stepper = make_controlled(eps_abs, eps_rel, stepper_type());

    constexpr double start_time = 0.0;
    constexpr double middle_time = 1e-9;
    constexpr double end_time = 2e-9;
    constexpr double dt = 1e-11; // Initial step for ADS control

    // not-m-normalizing LLG, has to be done in observer!:
    auto llg_regular = [&alpha, &state, &fieldterms](const af::array& m_in, af::array& dxdt, const double t) {
        state.t = t;
        state.m = m_in;
        const auto H_eff_in_Apm = value_Heff_in_Apm(fieldterms, state);
        dxdt = equations::LLG(alpha, state.m, H_eff_in_Apm);
    };

    struct observe_m {
        std::filesystem::path outdir_;
        explicit observe_m(std::filesystem::path outdir) : outdir_(std::move(outdir)) {}

        // this renorms m, thus taking by ref
        // NOT handling zero vals in m
        void operator()(af::array& m, double t) const {
            // Possible renorm location, when observer is called after every step, i.e. adaptive mode:
            // ... to renorm after every step, only works with integrate_adaptive()!
            m = util::normalize(m); // NOTE: This is a hack!

            const auto [mx, my, mz] = math::mean_3d<double>(m);

            // V2
            // // TODO this is slow!:
            // af::eval(m); // TODO This prevents slowing down whenusing .as(f64)
            // const auto [mx, my, mz] = math::mean_3d<double>(m.as(f64)); // This is slooow

            // V3
            // const auto mean = af::mean(af::mean(af::mean(m, 0), 1), 2).as(f64);
            // const auto mx = mean(0, 0, 0, 0).scalar<double>();
            // const auto my = mean(0, 0, 0, 1).scalar<double>();
            // const auto mz = mean(0, 0, 0, 2).scalar<double>();

            std::cout << t << "\t" << mx << "\t" << my << "\t" << mz << std::endl;
            std::ofstream stream(outdir_ / "m.dat", std::ios::app);
            stream.precision(12);
            stream << t << "\t" << mx << "\t" << my << "\t" << mz << std::endl;
        }
    };

    // Integrate
    steps.push_back(integrate_adaptive(stepper, llg_regular, m, start_time, middle_time, dt, observe_m{outdir}));
    std::cout << "steps=" << steps << std::endl;

    // Setting external field
    af::array external = af::constant(0.0, nx, ny, nz, 3, type);
    external(af::span, af::span, af::span, 0) = -24.6e-3 / constants::mu0;
    external(af::span, af::span, af::span, 1) = +4.3e-3 / constants::mu0;
    fieldterms.emplace_back(ExternalField(external));
    alpha = 0.02;

    steps.push_back(integrate_adaptive(stepper, llg_regular, m, middle_time, end_time, dt, observe_m{outdir}));

    std::cout << "steps=" << steps << std::endl;
    std::cout << "accum_steps=" << std::accumulate(steps.begin(), steps.end(), 0) << std::endl;
}
