#include "state.hpp"
#include "math.hpp"
#include "util/af_overloads.hpp" // for 'os << af::array'
#include "util/color_string.hpp"
#include "util/util.hpp"
#include "vtk_io.hpp"
#include <iomanip>
#include <stdexcept>
#include <utility>

namespace magnumaf {

/// Overloaded '+' operator adds an af::array to af::array this->m
State State::operator+(const af::array& a) const {
    // TODO: This line calls assignment operator of State Members which causes: warning: implicitly-declared ‘af::dim4&
    // af::dim4::operator=(const af::dim4&)’ is deprecated [-Wdeprecated-copy], as arrayfire does not explicitly define
    // default operators.
    State result = *this;
    result.m += a;
    return result;
    // This would confict with const qualifier:
    // this->m += a;
    // return *this;
}

std::ostream& operator<<(std::ostream& os, const State& state) {
    const auto mean_m = state.mean_m_as_afarray();
    os << state.t << '\t' << mean_m;
    return os;
}

af::array State::mean_m_as_afarray() const {
    if (Ms_field.isempty() or n_cells_ == 0) {
        return af::mean(af::mean(af::mean(m, 0), 1), 2);
    } else {
        return af::sum(af::sum(af::sum(m, 0), 1), 2) / n_cells_;
    }
}

std::array<double, 3> State::mean_m() const { return math::vec_components<double>(mean_m_as_afarray().as(f64)); }

auto State::mean_M_as_afarray() const -> af::array {
    af::array nzero = !af::iszero(util::vecnorm(m));
    const auto no_nonzero_cells = util::afvalue_u32(af::sum(af::sum(af::sum(nzero, 0), 1), 2));
    if (Ms_field.isempty()) {
        return af::sum(af::sum(af::sum(Ms * m, 0), 1), 2) / no_nonzero_cells;
    } else {
        return af::sum(af::sum(af::sum(get_Ms_field_in_vec_dims() * m, 0), 1), 2) / no_nonzero_cells;
    }
}

auto State::mean_M() const -> std::array<double, 3> {
    return math::vec_components<double>(mean_M_as_afarray().as(f64));
}

af::array State::get_Ms_as_field() const {
    if (Ms_field.isempty()) {
        return af::constant(Ms, m.dims(0), m.dims(1), m.dims(2), 1, m.type());
    } else {
        return Ms_field;
    }
}
af::array State::get_Ms_as_field_in_vector_dims() const {
    return af::tile(get_Ms_as_field(), 1, 1, 1, 3);
} //!< return Ms tiled to dims [nx, ny, nz, 3].

void State::set_Ms_field_if_m_minvalnorm_is_zero(const af::array& m, af::array& Ms_field) {
    // Initializes Ms_field if any entry of initial m has zero norm
    if (util::min_4d(util::vecnorm(m)) == 0) {
        if (verbose) {
            printf("%s in state.cpp: initial m has values with zero norm, "
                   "building Ms_field array\n",
                   color_string::info());
        }
        af::array nzero = !af::iszero(util::vecnorm(m));
        n_cells_ = util::afvalue_u32(af::sum(af::sum(af::sum(nzero, 0), 1), 2));
        if (Ms == 0) {
            printf("Wraning: State::set_Ms_field: State.Ms is used but set to "
                   "zero. It appears that you are using a legacy constuctor. "
                   "Please pass Ms in constructor!\n");
        }
        Ms_field = af::constant(this->Ms, nzero.dims(),
                                m.type()); // TODO this yields probem as Ms is not set in constuctor!
        Ms_field *= nzero;
    }
}

State::State(Mesh mesh, double Ms, af::array m, bool verbose, bool mute_warning)
    : mesh(mesh), m(std::move(m)), Ms(Ms), verbose(verbose), mute_warning(mute_warning) {
    util::normalize_inplace(this->m);
    set_Ms_field_if_m_minvalnorm_is_zero(this->m, this->Ms_field);
}

af::array check_Ms_field_dims(af::array Ms_field) {
    if (Ms_field.dims(3) == 1) {
        return Ms_field;
    } else if (Ms_field.dims(3) == 3) {
        std::cout << "Note: Legacy State.Ms_field dims with [nx, ny, nz , 3] used. Dim 3 will be cropped to 1. Now "
                     "please use [nx, ny, nz , 1]."
                  << std::endl;
        return Ms_field(af::span, af::span, af::span, 0);
    } else {
        throw std::invalid_argument("Ms_field dimension expected to be [nx, ny, nz , 1]");
    }
}

State::State(Mesh mesh, af::array Ms_field_in, af::array m_in, bool verbose, bool mute_warning)
    : mesh(mesh), m(std::move(m_in)), Ms_field(check_Ms_field_dims(std::move(Ms_field_in))), verbose(verbose),
      mute_warning(mute_warning) {
    util::normalize_inplace(this->m);
}

// No Mesh:
State::State(af::array m, double Ms, bool verbose, bool mute_warning)
    : State(Mesh{0, 0, 0, 0, 0, 0}, Ms, std::move(m), verbose, mute_warning) {}

State::State(af::array m, af::array Ms_field, bool verbose, bool mute_warning)
    : State(Mesh{0, 0, 0, 0, 0, 0}, std::move(Ms_field), std::move(m), verbose, mute_warning) {}

// Wrapping:
State::State(Mesh mesh, double Ms, long int m, bool verbose, bool mute_warning)
    : State(mesh, Ms, util::pywrap::make_copy_form_py(m), verbose, mute_warning) {}

// Wrapping only, memory management to be done by python:
State::State(Mesh mesh, long int Ms_field_ptr, long int m, bool verbose, bool mute_warning)
    : State(mesh, util::pywrap::make_copy_form_py(Ms_field_ptr), util::pywrap::make_copy_form_py(m), verbose,
            mute_warning) {}

void State::Normalize() { this->m = util::normalize(this->m); }

void State::set_m(long int aptr) {
    m = util::pywrap::make_copy_form_py(aptr);
    util::normalize_inplace(this->m);
}

long int State::get_m_addr() const { return util::pywrap::send_copy_to_py(m); }

void State::set_Ms_field(long int aptr) { Ms_field = util::pywrap::make_copy_form_py(aptr); }

long int State::wrapping_get_Ms_field() const { return util::pywrap::send_copy_to_py(Ms_field); }

void State::write_vti(const std::string& outputname) const { vti_writer_micro(m.as(f64), mesh, outputname); }
void State::_vti_writer_atom(std::string outputname) const { vti_writer_atom(m.as(f64), mesh, std::move(outputname)); }
void State::_vti_reader(const std::string& inputname) { vti_reader(m, mesh, inputname); }

double State::meani(const int i) const {
    auto m = mean_m();
    return m[i];
}

} // namespace magnumaf
