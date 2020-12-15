#include "state.hpp"
#include "func.hpp"
#include "misc.hpp"
#include "vtk_IO.hpp"
#include <iomanip>

namespace magnumafcpp {

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
    auto [mx, my, mz] = state.mean_m();
    os << state.t << "\t" << mx << "\t" << my << "\t" << mz;
    return os;
}

af::array State::mean_m_as_afarray() const {
    if (Ms_field.isempty() or n_cells_ == 0) {
        return af::mean(af::mean(af::mean(m, 0), 1), 2);
    } else {
        return af::sum(af::sum(af::sum(m, 0), 1), 2) / n_cells_;
    }
}

std::array<double, 3> State::mean_m() const {
    af::array mean_m_af = mean_m_as_afarray().as(f64);
    const double mx = mean_m_af(0, 0, 0, 0).scalar<double>();
    const double my = mean_m_af(0, 0, 0, 1).scalar<double>();
    const double mz = mean_m_af(0, 0, 0, 2).scalar<double>();
    return {mx, my, mz};
}

af::array State::get_Ms_field() const {
    if (Ms_field.isempty()) {
        return af::constant(Ms, m.dims(0), m.dims(1), m.dims(2), 1, m.type());
    } else {
        // TODO enable after changing Ms_field.dims(4) to 1:
        // TODO return Ms_field;
        if (Ms_field.dims(3) == 1) {
            return Ms_field;
        }
        // Current workaround:
        else {
            return Ms_field(af::span, af::span, af::span, 0);
        }
    }
}
af::array State::get_Ms_field_in_vector_dims() const {
    return af::tile(get_Ms_field(), 1, 1, 1, 3);
} //!< return Ms tiled to dims [nx, ny, nz, 3].

void State::set_Ms_field_if_m_minvalnorm_is_zero(const af::array& m, af::array& Ms_field) {
    // Initializes Ms_field if any entry of initial m has zero norm
    if (min_4d(vecnorm(m)) == 0) {
        if (verbose) {
            printf("%s in state.cpp: initial m has values with zero norm, "
                   "building Ms_field array\n",
                   Info());
        }
        af::array nzero = !af::iszero(vecnorm(m));
        n_cells_ = afvalue_u32(af::sum(af::sum(af::sum(nzero, 0), 1), 2));
        if (Ms == 0)
            printf("Wraning: State::set_Ms_field: State.Ms is used but set to "
                   "zero. It appears that you are using a legacy constuctor. "
                   "Please pass Ms in constructor!\n");
        Ms_field = af::constant(this->Ms, nzero.dims(),
                                m.type()); // TODO this yields probem as Ms is not set in constuctor!
        Ms_field *= nzero;
        Ms_field = af::tile(Ms_field, 1, 1, 1, 3);
    }
}

// void State::check_nonequispaced_discretization(){
//    if ( this->material.A != 0 && this->material.Ku1 != 0) { // TODO implement
//    better way of checking
//        const double max_allowed_cellsize =
//        sqrt(this->material.A/this->material.Ku1); const double max_dz =
//        *std::max_element(nonequimesh.z_spacing.begin(),
//        nonequimesh.z_spacing.end()); if (verbose && (this->mesh.dx >
//        max_allowed_cellsize || this->mesh.dy > max_allowed_cellsize || max_dz
//        > max_allowed_cellsize )){
//            if( ! mute_warning) printf("%s State::check_discretization: cell
//            size is too large (greater than sqrt(A/Ku1)\n", Warning());
//        }
//    }
//}

// void State::check_m_norm(double tol) { // allowed norm is 1 or 0 (for no Ms_field)
//    af::array one_when_value_is_zero = af::iszero(vecnorm(m));
//    double meannorm = af::mean(af::mean(af::mean(af::mean(vecnorm(m) + 1. * one_when_value_is_zero, 0), 1), 2), 3)
//                          .as(f64)
//                          .scalar<double>();
//    if ((std::fabs(meannorm - 1.) > tol) && (this->mute_warning == false)) {
//        printf("%s State::check_m_norm: non-zero parts of the magnetization are "
//               "not normalized to 1! Results won't be physically meaningfull.\n",
//               Warning());
//    }
//}

// long int State::get_m_addr(){
//    u_out = this->m.copy();
//    return (long int) m_out.get();
//}
//

// Micromagnetic:
State::State(Mesh mesh, double Ms, af::array m, bool verbose, bool mute_warning)
    : mesh(mesh), m(m), Ms(Ms), verbose(verbose), mute_warning(mute_warning) {
    normalize_inplace(this->m);
    set_Ms_field_if_m_minvalnorm_is_zero(this->m, this->Ms_field);
    // check_discretization();
}

State::State(Mesh mesh, af::array Ms_field, af::array m, bool verbose, bool mute_warning)
    : mesh(mesh), m(m), Ms_field(Ms_field.dims(3) == 1 ? af::tile(Ms_field, 1, 1, 1, 3) : Ms_field), verbose(verbose),
      mute_warning(mute_warning) {
    if (Ms_field.dims(3) == 3) {
        printf("%s State: You are using legacy dimension [nx, ny, nz, 3] for "
               "Ms, please now use scalar field dimensions [nx, ny, nz, 1].\n",
               Warning());
    }

    normalize_inplace(this->m);
}

// No Mesh:
State::State(af::array m, double Ms, bool verbose, bool mute_warning)
    : m(m), Ms(Ms), verbose(verbose), mute_warning(mute_warning) {
    normalize_inplace(this->m);
    set_Ms_field_if_m_minvalnorm_is_zero(this->m, this->Ms_field);
}

State::State(af::array m, af::array Ms_field, bool verbose, bool mute_warning)
    : m(m), Ms_field(Ms_field.dims(3) == 1 ? af::tile(Ms_field, 1, 1, 1, 3) : Ms_field), verbose(verbose),
      mute_warning(mute_warning) {
    if (Ms_field.dims(3) == 3) {
        printf("%s State: You are using legacy dimension [nx, ny, nz, 3] for "
               "Ms, please now use scalar field dimensions [nx, ny, nz, 1].\n",
               Warning());
    }

    normalize_inplace(this->m);
}

// Wrapping:
State::State(Mesh mesh, double Ms, long int m, bool verbose, bool mute_warning)
    : mesh(mesh), m(*(new af::array(*((void**)m)))), Ms(Ms), verbose(verbose), mute_warning(mute_warning) {
    normalize_inplace(this->m);
    set_Ms_field_if_m_minvalnorm_is_zero(this->m, this->Ms_field);
    // check_discretization();
}

// Wrapping only, memory management to be done by python:
State::State(Mesh mesh, long int Ms_field_ptr, long int m, bool verbose, bool mute_warning)
    : mesh(mesh), m(*(new af::array(*((void**)m)))),
      Ms_field((*(new af::array(*((void**)Ms_field_ptr)))).dims(3) == 1
                   ? af::tile(*(new af::array(*((void**)Ms_field_ptr))), 1, 1, 1, 3)
                   : *(new af::array(*((void**)Ms_field_ptr)))),
      verbose(verbose),

      mute_warning(mute_warning) {
    if ((*(new af::array(*((void**)Ms_field_ptr)))).dims(3) == 3) {
        printf("%s State: You are using legacy dimension [nx, ny, nz, 3] for "
               "Ms, please now use scalar field dimensions [nx, ny, nz, 1].\n",
               Warning());
    }
    normalize_inplace(this->m);
}

void State::Normalize() { this->m = normalize(this->m); }

void State::set_m(long int aptr) {
    void** a = (void**)aptr;
    m = *(new af::array(*a));
    normalize_inplace(this->m);
}

long int State::get_m_addr() const {
    af::array* a = new af::array(m);
    return (long int)a->get();
}

void State::set_Ms_field(long int aptr) {
    void** a = (void**)aptr;
    Ms_field = *(new af::array(*a)); // TODO rename Ms_field -> micro_Ms_field
}

long int State::wrapping_get_Ms_field() const {
    af::array* a = new af::array(Ms_field);
    return (long int)a->get();
}

void State::write_vti(std::string outputname) const { vti_writer_micro(m.as(f64), mesh, outputname); }
void State::_vti_writer_atom(std::string outputname) const { vti_writer_atom(m.as(f64), mesh, outputname); }
void State::_vti_reader(std::string inputname) { vti_reader(m, mesh, inputname); }

double State::meani(const int i) const {
    auto m = mean_m();
    return m[i];
}

} // namespace magnumafcpp
