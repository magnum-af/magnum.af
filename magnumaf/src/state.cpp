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

unsigned get_number_of_nonzero_cells_form_vectorfield(const af::array& a) {
    af::array one_if_val_is_zero_else_zero = !af::iszero(util::vecnorm(a));
    const auto number_of_nonzero_cells =
        af::sum(af::sum(af::sum(one_if_val_is_zero_else_zero, 0), 1), 2).scalar<unsigned>();
    return number_of_nonzero_cells;
}

af::array State::mean_m_as_afarray() const {
    auto number_of_nonzero_cells = get_number_of_nonzero_cells_form_vectorfield(m);
    auto mean_m = af::sum(af::sum(af::sum(m, 0), 1), 2) / number_of_nonzero_cells;
    return mean_m;
}

std::array<double, 3> State::mean_m() const { return math::vec_components<double>(mean_m_as_afarray().as(f64)); }

auto State::mean_M_as_afarray() const -> af::array {
    auto number_of_nonzero_cells = get_number_of_nonzero_cells_form_vectorfield(m);
    if (Ms_field.isempty()) {
        return af::sum(af::sum(af::sum(Ms * m, 0), 1), 2) / number_of_nonzero_cells;
    } else {
        return af::sum(af::sum(af::sum(get_Ms_field_in_vec_dims() * m, 0), 1), 2) / number_of_nonzero_cells;
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
    : mesh(mesh), m(util::normalize(m)), Ms(Ms), verbose(verbose), mute_warning(mute_warning) {
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
    : mesh(mesh), m(util::normalize(m_in)), Ms_field(check_Ms_field_dims(std::move(Ms_field_in))), verbose(verbose),
      mute_warning(mute_warning) {}

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

void State::set_m(long int aptr) { m = util::normalize(util::pywrap::make_copy_form_py(aptr)); }

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
