#include "field_terms/micro/dmi_field.hpp"
#include "math.hpp"
#include "util/util.hpp"

namespace magnumafcpp {

void showdims(const af::array& a) {
    std::cout << "Exchange matrix: dims=" << a.dims(0) << "\t" << a.dims(1) << "\t" << a.dims(2) << "\t" << a.dims(3)
              << std::endl;
}
void apply_boundary_condition(af::array& hfield, const State& state);
void correct_edges(af::array& out, const af::array& in, Mesh mesh);

DmiField::DmiField(double D, std::array<double, 3> D_axis) : D(D), D_axis(D_axis) {}

DmiField::DmiField(af::array D_constants, std::array<double, 3> D_axis)
    : D_constants(D_constants.dims(3) == 1 ? af::tile(D_constants, 1, 1, 1, 3) : std::move(D_constants)),
      D_axis(D_axis) {
    if (this->D_constants.dims(3) == 3) {
        printf("\33[1;31mWarning:\33[0m DmiField: You are using legacy "
               "dimension [nx, ny, nz, 3] for D, please now use scalar field "
               "dimensions [nx, ny, nz, 1].\n");
    }
}

DmiField::DmiField(double D, double D_axis_x, double D_axis_y, double D_axis_z)
    : D(D), D_axis({D_axis_x, D_axis_y, D_axis_z}) {}

DmiField::DmiField(long int D_constants_ptr, double D_axis_x, double D_axis_y, double D_axis_z)
    : DmiField(util::pywrap::make_copy_form_py(D_constants_ptr), {D_axis_x, D_axis_y, D_axis_z}) {}

///
/// Bulk Dzyaloshinskii–Moriya interaction.
/// Calculates effective field [1, eq (7)]
/// \f[
///     \boldsymbol{H}_{DMI} = \frac{2 D}{\mu_0 M_s} [(\nabla \cdot
///     \boldsymbol{m}) \hat{n} - \nabla (\boldsymbol{m} \cdot \hat{n}) ]
/// \f]
/// , where \f$ \hat{n} \f$ defines the dmi unit-vector which is passed in the
/// consturctor.
///
/// [1] Rohart S and Thiaville A 2013 Skyrmion confinement in ultrathin film
/// nanostructures in the presence of Dzyaloshinskii–Moriya interaction Phys.
/// Rev. B 88 184422
///
af::array DmiField::impl_H_in_Apm(const State& state) const {
    // Normal vector
    double norm = sqrt(pow(D_axis[0], 2) + pow(D_axis[1], 2) + pow(D_axis[2], 2));
    af::array n = af::array(state.mesh.nx, state.mesh.ny, state.mesh.nz, 3, f64);
    n(af::span, af::span, af::span, 0) = D_axis[0] / norm;
    n(af::span, af::span, af::span, 1) = D_axis[1] / norm;
    n(af::span, af::span, af::span, 2) = D_axis[2] / norm;
    // print("n", n);

    // initialize finite difference first order derivative filter
    // Central finite difference in 2nd order
    // https://en.wikipedia.org/wiki/Finite_difference_coefficient
    af::array filtr_fd1 = af::constant(0.0, 3, 3, 3, 3, f64);
    // dmx/dx
    filtr_fd1(0, 1, 1, 0) = 1 / (2. * state.mesh.dx);
    filtr_fd1(2, 1, 1, 0) = -1 / (2. * state.mesh.dx);

    // dmy/dy
    filtr_fd1(1, 0, 1, 1) = 1 / (2. * state.mesh.dy);
    filtr_fd1(1, 2, 1, 1) = -1 / (2. * state.mesh.dy);

    // dmz/dz
    filtr_fd1(1, 1, 0, 2) = 1 / (2. * state.mesh.dz);
    filtr_fd1(1, 1, 2, 2) = -1 / (2. * state.mesh.dz);

    // First: n(div m)
    // Gradient and edges
    af::array first = convolve(state.m, filtr_fd1, AF_CONV_DEFAULT, AF_CONV_SPATIAL);

    // correct_edges(first, state.m);
    // TODO causes segfault//
    // apply_boundary_condition(first, state);
    // Make Divergence
    first = sum(first, 3);
    first = tile(first, 1, 1, 1, 3);
    first = n * first;

    // Second: fd1(n . m)
    // Dot product
    af::array second = sum(n * state.m, 3);
    // Expand for fd1 convolution
    second = tile(second, 1, 1, 1, 3);
    second = convolve(second, filtr_fd1, AF_CONV_DEFAULT, AF_CONV_SPATIAL);
    // correct_edges(second, state.m);
    // TODO causes segfault//
    // apply_boundary_condition(second, state);

    // if (state.Ms_field.isempty()){
    //  return 2.* material.D/(constants::mu0*state.Ms) * (first-second);//Note:
    //  Js=mu0*Ms
    //}
    // else{
    //  return 2.* material.D/(constants::mu0*state.Ms_field) *
    //  (first-second);//Note: Js=mu0*Ms
    //}

    if (state.Ms_field.isempty() && this->D_constants.isempty()) {
        return (2. * this->D) / (constants::mu0 * state.Ms) * (first - second);
    } else if (!state.Ms_field.isempty() && this->D_constants.isempty()) {
        af::array heff = (2. * this->D) / (constants::mu0 * state.Ms_field) * (first - second);
        af::replace(heff, state.Ms_field != 0, 0); // set all cells where Ms==0 to 0
        return heff;
    } else if (state.Ms_field.isempty() && !this->D_constants.isempty()) {
        return (2. * this->D_constants) / (constants::mu0 * state.Ms) * (first - second);
    } else {
        af::array heff = (2. * this->D_constants) / (constants::mu0 * state.Ms_field) * (first - second);
        af::replace(heff, state.Ms_field != 0, 0); // set all cells where Ms==0 to 0
        return heff;
    }
}

// EXPERIMENTAL
// This version assumes the same boundaries as in the exchange filed of 70Lines
// of Numpy On the edges, we assume (not existing) outer cells with the same
// values as the boundary cells
// in the middle ...|-1/2|  0 |1/2|...
// on the edges         ||  0 |1/2|...
// we want         |-1/2||  0 |1/2|...
// thus we take         ||-1/2|1/2|...
// so after the convolution we have to add the edges with -1/(2*dx)

void correct_edges(af::array& out, const af::array& in, Mesh mesh) {
    // Lower x edge:
    out(0, af::span, af::span, 0) += -0.5 * in(0, af::span, af::span, 0) / mesh.dx;
    // Upper x edge:
    out(-1, af::span, af::span, 0) -= -0.5 * in(-1, af::span, af::span, 0) / mesh.dx;

    // Lower y edge:
    out(af::span, 0, af::span, 1) += -0.5 * in(af::span, 0, af::span, 1) / mesh.dy;
    // Upper y edge:
    out(af::span, -1, af::span, 1) -= -0.5 * in(af::span, -1, af::span, 1) / mesh.dy;

    // z
    if (in.dims(2) == 1) {
        out(af::span, af::span, af::span, 2) = 0.;
    } else {
        // Lower z edge:
        out(af::span, af::span, 0, 2) += -0.5 * in(af::span, af::span, 0, 2) / mesh.dz;
        // Upper z edge:
        out(af::span, af::span, -1, 2) -= -0.5 * in(af::span, af::span, -1, 2) / mesh.dz;
    }
}

//--------------------------------------------------------------------------------------------------------------------------------
// EXPERIMENTAL
// TODO unit tests !
// Boundary Conditions according to DOI: 10.1103/PhysRevB.88.184422 Skyrmion
// confinement in ultrathin film nanostructures ... dm/dn = 1/xi (n_DM x
// n_surface) x m ; xi = 2 A/D here: dm/dn = 1/xi (D_axis x n_surface) x m ; xi
// = 2 A/D
void DmiField::apply_boundary_condition(af::array& hfield, const State& state) const {
    // DM Vector:
    const af::array n_DM(1, 1, 1, 3, this->D_axis.data());
    double A = 15e-12; // TODO set exchange A as class member and pass in constructor

    if (state.m.dims(0) == 1) {
        hfield(af::span, af::span, af::span, 0) = 0.;
    } else {
        // low x boundary:
        // n_surface=(-1, 0, 0)
        // dm/dn=dm/d(-x)= (m_-1 - m_1) / 2 * dx = 1/xi (D_axis x n_surface) x m
        // => m_-1 = m_1 + 2 * dx * 1/xi (D_axis x n_surface) x m_0
        const double c_n_x_surface_low[] = {-1, 0, 0};
        const af::array n_x_surface_low(1, 1, 1, 3, c_n_x_surface_low); // normal vector to the x-surface at
                                                                        // boundary with lower index i=0
        const af::array n_DMxn_x_surf_low = tile(math::cross4(n_DM, n_x_surface_low), 1, state.m.dims(1), state.m.dims(2), 1);
        const af::array x_minus_1 =
            state.m(1, af::span, af::span, 0) +
            2 * state.mesh.dx * (this->D / (2 * A)) *
                math::cross4(n_DMxn_x_surf_low, state.m(0, af::span, af::span, af::span))(af::span, af::span, af::span, 0);
        hfield(0, af::span, af::span, 0) += -0.5 * x_minus_1 / state.mesh.dx; // Minus due to: (m_{i+1} - m_{i-1})/(
                                                                              // 2*dx )  with m_{i-1} being replaced

        // high x boundary:
        // n_surface=(1, 0, 0)
        // dm/dn=dm/d(x)= (m_{n+1} - m_{n-1}) / 2 * dx = 1/xi (D_axis x
        // n_surface) x m_n
        // => m_{i+1} = m_{i-1} + 2 * dx * 1/xi (D_axis x n_surface) x m_i
        const double c_n_x_surface_high[] = {1, 0, 0};
        const af::array n_x_surface_high(1, 1, 1, 3,
                                         c_n_x_surface_high); // normal vector to the x-surface at boundary
                                                              // with higher index i=0
        const af::array n_DMxn_x_surf_high =
            tile(math::cross4(n_DM, n_x_surface_high), 1, state.m.dims(1), state.m.dims(2), 1);
        const af::array x_i_plus_1 =
            state.m(-1, af::span, af::span, 0) +
            2 * state.mesh.dx * (this->D / (2 * A)) *
                math::cross4(n_DMxn_x_surf_high, state.m(-1, af::span, af::span, af::span))(af::span, af::span, af::span, 0);
        hfield(-1, af::span, af::span, 0) += 0.5 * x_i_plus_1 / state.mesh.dx;
    }

    if (state.m.dims(1) == 1) {
        hfield(af::span, af::span, af::span, 1) = 0.;
    } else {
        // low y boundary:
        // n_surface=(0, -1, 0)
        // dm/dn=dm/d(-y)= (m_-1 - m_1) / 2 * dx = 1/xi (D_axis x n_surface) x
        // m_0
        // => m_-1 = m_1 + 2 * dx * 1/xi (D_axis x n_surface) x m_0
        const double c_n_y_surface_low[] = {0, -1, 0};
        const af::array n_y_surface_low(1, 1, 1, 3, c_n_y_surface_low); // normal vector to the y-surface at
                                                                        // boundary with lower indey i=0
        const af::array n_DMxn_y_surf_low = tile(math::cross4(n_DM, n_y_surface_low), state.m.dims(0), 1, state.m.dims(2), 1);
        const af::array y_minus_1 =
            state.m(af::span, 1, af::span, 1) +
            2 * state.mesh.dy * (this->D / (2 * A)) *
                math::cross4(n_DMxn_y_surf_low, state.m(af::span, 0, af::span, af::span))(af::span, af::span, af::span, 1);
        hfield(af::span, 0, af::span, 1) += -0.5 * y_minus_1 / state.mesh.dy; // Minus due to: (m_{i+1} - m_{i-1})/(
                                                                              // 2*dy )  with m_{i-1} being replaced

        // high y boundary:
        // n_surface=(0, 1, 0)
        const double c_n_y_surface_high[] = {0, 1, 0};
        const af::array n_y_surface_high(1, 1, 1, 3,
                                         c_n_y_surface_high); // normal vector to the y-surface at boundary
                                                              // with higher indey i=0
        const af::array n_DMxn_y_surf_high =
            tile(math::cross4(n_DM, n_y_surface_high), state.m.dims(0), 1, state.m.dims(2), 1);
        const af::array y_i_plus_1 =
            state.m(af::span, -1, af::span, 1) +
            2 * state.mesh.dy * (this->D / (2 * A)) *
                math::cross4(n_DMxn_y_surf_high, state.m(af::span, -1, af::span, af::span))(af::span, af::span, af::span, 1);
        hfield(af::span, -1, af::span, 1) += 0.5 * y_i_plus_1 / state.mesh.dy; // Minus due to: (m_{i+1} - m_{i-1})/(
                                                                               // 2*dy )  with m_{i-1} being replaced
    }

    // z
    if (state.m.dims(2) == 1) {
        hfield(af::span, af::span, af::span, 2) = 0.;
    } else {
        // low z boundary:
        // n_surface=(0, 0, -1)
        // dm/dn=dm/d(-z)= (m_-1 - m_1) / 2 * dz = 1/xi (D_axis x n_surface) x m
        // => m_-1 = m_1 + 2 * dz * 1/xi (D_axis x n_surface) x m_0
        const double c_n_z_surface_low[] = {0, 0, -1};
        const af::array n_z_surface_low(1, 1, 1, 3, c_n_z_surface_low); // normal vector to the z-surface at
                                                                        // boundary with lower index i=0
        const af::array n_DMxn_z_surf_low = tile(math::cross4(n_DM, n_z_surface_low), state.m.dims(0), state.m.dims(1), 1, 1);
        const af::array z_minus_1 =
            state.m(af::span, af::span, 1, 2) +
            2 * state.mesh.dz * (this->D / (2 * A)) *
                math::cross4(n_DMxn_z_surf_low, state.m(af::span, af::span, 0, af::span))(af::span, af::span, af::span, 2);
        hfield(af::span, af::span, 0, 2) += -0.5 * z_minus_1 / state.mesh.dz; // Minus due to: (m_{i+1} - m_{i-1})/(
                                                                              // 2*dz )  with m_{i-1} being replaced

        // high z boundary:
        // n_surface=(1, 0, 0)
        // dm/dn=dm/d(z)= (m_{n+1} - m_{n-1}) / 2 * dz = 1/xi (D_axis x
        // n_surface) x m_n
        // => m_{i+1} = m_{i-1} + 2 * dz * 1/xi (D_axis x n_surface) x m_i
        const double c_n_z_surface_high[] = {0, 0, 1};
        const af::array n_z_surface_high(1, 1, 1, 3,
                                         c_n_z_surface_high); // normal vector to the z-surface at boundary
                                                              // with higher index i=0
        const af::array n_DMxn_z_surf_high =
            tile(math::cross4(n_DM, n_z_surface_high), state.m.dims(0), state.m.dims(1), 1, 1);
        const af::array z_i_plus_1 =
            state.m(af::span, af::span, -1, 2) +
            2 * state.mesh.dz * (this->D / (2 * A)) *
                math::cross4(n_DMxn_z_surf_high, state.m(af::span, af::span, -1, af::span))(af::span, af::span, af::span, 2);
        hfield(af::span, af::span, -1, 2) += 0.5 * z_i_plus_1 / state.mesh.dz;
    }
}
} // namespace magnumafcpp
