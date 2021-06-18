#include "magnum_af.hpp"
#include <algorithm>
#include <boost/numeric/odeint.hpp>
#include <iostream>
#include <utility>

#include <variant>

// Necessary template specialization for af::array
// Adapted form boost/numeric/odeint/algebra/vector_space_algebra.hpp
// and boost/numeric/odeint/external/vexcl/vexcl_norm_inf.hpp
namespace boost::numeric::odeint {
template <> struct vector_space_norm_inf<af::array> {
    typedef double result_type; // typedef definition is explicitly needed here for odeint internals.
    result_type operator()(const af::array& x) const {
        return af::max(af::max(af::max(af::max(af::abs(x), 0), 1), 2), 3).as(f64).scalar<double>();
    }
};
} // namespace boost::numeric::odeint

namespace magnumafcpp {

// Oneliner:
template <typename T> af::array Heff_in_Apm_onliner(const T& fieldterms, const State& state) {
    return std::accumulate(std::begin(fieldterms) + 1, std::end(fieldterms),
                           std::visit([&state](auto&& vari) { return vari.H_in_Apm(state); }, fieldterms[0]),
                           [&state](const auto& sum, const auto& elem) {
                               return sum + std::visit([&state](auto&& vari) { return vari.H_in_Apm(state); }, elem);
                           });
}

template <typename T> af::array Heff_in_Apm(const T& fieldterms, const State& state) {
    const auto initial_H = std::visit([&state](auto&& vari) { return vari.H_in_Apm(state); }, fieldterms[0]);
    const auto accumul_H = [&state](const auto& sum, const auto& elem) {
        return sum + std::visit([&state](auto&& vari) { return vari.H_in_Apm(state); }, elem);
    };
    return std::accumulate(std::begin(fieldterms) + 1, std::end(fieldterms), initial_H, accumul_H);
}
} // namespace magnumafcpp

// // Possible lamda implementation:
// const auto Heff_in_Apm_lamda = [](const FieldtermsVari& fieldterms, const State& state) {
//     const auto initial_H = std::visit([&state](auto&& vari) { return vari.H_in_Apm(state); }, fieldterms[0]);
//     const auto accumul_H = [&state](const auto& sum, const auto& elem) {
//         return sum + std::visit([&state](auto&& vari) { return vari.H_in_Apm(state); }, elem);
//     };
//     return std::accumulate(std::begin(fieldterms) + 1, std::end(fieldterms), initial_H, accumul_H);
// };

int main(int argc, char** argv) {
    std::cout << "Start" << std::endl;

    using namespace magnumafcpp;
    const auto [outdir_tmp, posargs] = ArgParser(argc, argv).outdir_posargs;
    const auto outdir = outdir_tmp;

    constexpr double x_magaf = 5.e-7, y = 1.25e-7, z = 3.e-9;
    constexpr int nx = 100, ny = 25, nz = 1;
    constexpr Mesh mesh{nx, ny, nz, x_magaf / nx, y / ny, z / nz};

    double alpha = 1.;
    constexpr double Ms = 8e5;
    constexpr double A = 1.3e-11;

    using FieldtermVari = std::variant<DemagField, ExchangeField, ExternalField>;
    using FieldtermsVari = std::vector<FieldtermVari>;

    FieldtermsVari fieldterms_vari{DemagField(mesh, true, true, 0), ExchangeField(A)};
    // FieldtermsVari fieldterms_vari{DemagField(mesh, true, true, 0), ExchangeField(A)};

    // FieldtermsVari fieldterms_vari;
    // fieldterms_vari.push_back(DemagField(mesh, true, true, 0));
    // fieldterms_vari.push_back(ExchangeField(A));

    // FieldtermsVari fieldterms_vari;
    // auto dmag = DemagField(mesh, true, true, 0);
    // auto exch = ExchangeField(A);
    // fieldterms_vari.push_back(dmag);
    // fieldterms_vari.push_back(exch);

    constexpr auto type = af::dtype::f64;
    af::array m = af::constant(0, mesh::dims_v(mesh), type);
    m(0, af::span, af::span, 1) = 1;
    m(af::seq(1, af::end - 1), af::span, af::span, 0) = 1;
    m(-1, af::span, af::span, 1) = 1;

    State state(mesh, Ms, m);
    state.write_vti(outdir / "minit");

    // not-m-normalizing LLG, has to be done in observer!:
    auto llg_regular = [&alpha, &state, &fieldterms_vari](const af::array& m_in, af::array& dxdt, const double t) {
        state.t = t;
        state.m = m_in;
        const auto H_eff_in_Apm = Heff_in_Apm(fieldterms_vari, state);
        dxdt = equations::LLG(alpha, state.m, H_eff_in_Apm);
    };

    struct observe_m {
        std::filesystem::path outdir_;
        observe_m(std::filesystem::path outdir) : outdir_(std::move(outdir)) {}

        // this renorms m, thus taking by ref
        // NOT handling zero vals in m
        void operator()(af::array& m, double t) {
            // Possible renorm location, when observer is called after every step, i.e. adaptive mode:
            // ... to renorm after every step, only works with integrate_adaptive()!
            m = util::normalize(m); // NOTE: This is a hack!
            const auto mean = af::mean(af::mean(af::mean(m, 0), 1), 2).as(f64);
            const auto mx = mean(0, 0, 0, 0).scalar<double>();
            const auto my = mean(0, 0, 0, 1).scalar<double>();
            const auto mz = mean(0, 0, 0, 2).scalar<double>();
            std::cout << t << "\t" << mx << "\t" << my << "\t" << mz << std::endl;
            std::ofstream stream(outdir_ / "m.dat", std::ios::app);
            stream.precision(12);
            stream << t << "\t" << mx << "\t" << my << "\t" << mz << std::endl;
        }
    };

    std::ofstream stream(outdir / "m.dat"); // clears old file
    stream.precision(12);
    std::vector<std::size_t> steps;
    namespace ode = boost::numeric::odeint;

    typedef ode::runge_kutta_dopri5<af::array, double, af::array, double, ode::vector_space_algebra> stepper_type;
    constexpr double eps_abs = 1e-6;
    constexpr double eps_rel = 1e-6;
    const auto stepper = make_controlled(eps_abs, eps_rel, stepper_type());

    constexpr double start_time = 0.0;
    constexpr double middle_time = 1e-9;
    constexpr double end_time = 2e-9;
    constexpr double dt = 1e-11; // Initial step for ADS control

    steps.push_back(integrate_adaptive(stepper, llg_regular, m, start_time, middle_time, dt, observe_m{outdir}));
    std::cout << "steps=" << steps << std::endl;

    // Setting external field
    af::array external = af::constant(0.0, nx, ny, nz, 3, type);
    external(af::span, af::span, af::span, 0) = -24.6e-3 / constants::mu0;
    external(af::span, af::span, af::span, 1) = +4.3e-3 / constants::mu0;
    fieldterms_vari.push_back(ExternalField(external));
    alpha = 0.02;

    steps.push_back(integrate_adaptive(stepper, llg_regular, m, middle_time, end_time, dt, observe_m{outdir}));

    std::cout << "steps=" << steps << std::endl;
    std::cout << "accum_steps=" << std::accumulate(steps.begin(), steps.end(), 0) << std::endl;
}
