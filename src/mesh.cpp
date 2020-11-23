#include "mesh.hpp"
#include "func.hpp"
#include <cmath>

namespace magnumafcpp {

Mesh::Mesh(unsigned nx, unsigned ny, unsigned nz, double dx, double dy, double dz)
    : n0(nx), n1(ny), n2(nz), dx(dx), dy(dy), dz(dz), n0_exp(2 * n0), n1_exp(2 * n1), n2_exp((n2 == 1) ? 1 : 2 * n2),
      dims(af::dim4(n0, n1, n2, 3)), dims_scalar(af::dim4(n0, n1, n2, 1)),
      dims_expanded(af::dim4(n0_exp, n1_exp, n2_exp, 3)) {}

std::ostream& operator<<(std::ostream& os, const Mesh& mesh) {
    os << "n0=" << mesh.n0 << " n1=" << mesh.n1 << " n2=" << mesh.n2 << " dx=" << mesh.dx << " dy=" << mesh.dy
       << " dz=" << mesh.dz << " n0_exp=" << mesh.n0_exp << " n1_exp=" << mesh.n1_exp << " n2_exp=" << mesh.n2_exp;
    return os;
}

af::array Mesh::skyrmconf(const bool point_up) {
    // Returns a initial configuration to be relaxed into a skyrmion
    // if point_up is true, skyrmion centers points in +z, if false in -z
    af::array m = af::constant(0.0, this->n0, this->n1, this->n2, 3, f64);
    if (point_up) {
        m(af::span, af::span, af::span, 2) = 1.;
    } else {
        m(af::span, af::span, af::span, 2) = -1.;
    }
    for (unsigned ix = 0; ix < this->n0; ix++) {
        for (unsigned iy = 0; iy < this->n1; iy++) {
            const double rx = double(ix) - this->n0 / 2.;
            const double ry = double(iy) - this->n1 / 2.;
            const double r = sqrt(pow(rx, 2) + pow(ry, 2));
            if (r > this->n0 / 4.) {
                if (point_up) {
                    m(ix, iy, af::span, 2) = -1.;
                } else {
                    m(ix, iy, af::span, 2) = 1.;
                }
            }
        }
    }
    return m;
}

/// Initializes a homogeneous magnetic field pointing in the direction of \param
/// vector within the largest cylinder which fits into the mesh.
// TODO should be rename to cylinder
af::array Mesh::ellipse(std::array<double, 3> vector, const bool verbose) {
    const double norm = std::sqrt(std::pow(vector[0], 2) + std::pow(vector[1], 2) + std::pow(vector[2], 2));
    vector[0] = vector[0] / norm;
    vector[1] = vector[1] / norm;
    vector[2] = vector[2] / norm;
    if (verbose)
        std::cout << "Mesh::ellipse: norm=" << norm << std::endl;
    if (verbose)
        std::cout << "Mesh::ellipse: normalized vector={" << vector[0] << ", " << vector[1] << ", " << vector[2] << "}"
                  << std::endl;

    long unsigned cells_within = 0;
    af::array m = af::constant(0.0, this->n0, this->n1, this->n2, 3, f64);
    for (unsigned ix = 0; ix < this->n0; ix++) {
        for (unsigned iy = 0; iy < this->n1; iy++) {
            const double a = (double)(this->n0 / 2);
            const double b = (double)(this->n1 / 2);
            const double rx = double(ix) - this->n0 / 2.;
            const double ry = double(iy) - this->n1 / 2.;
            const double r = pow(rx, 2) / pow(a, 2) + pow(ry, 2) / pow(b, 2);
            if (r < 1) {
                m(ix, iy, af::span, 0) = vector[0];
                m(ix, iy, af::span, 1) = vector[1];
                m(ix, iy, af::span, 2) = vector[2];
                cells_within++;
            }
        }
    }
    if (verbose)
        std::cout << "Info: Mesh::ellipse(): cells within cylinder = " << cells_within
                  << ", which should be approx a*b*M_PI*n2 = " << this->n0 / 2 * this->n1 / 2 * M_PI * this->n2
                  << std::endl;
    return m;
}

af::array Mesh::ellipse(const unsigned xyz, const bool positive_direction) {
    // Returns an initial elliptical magnetization
    // n_cells gives number of cells with non-zero Ms
    // xyz gives direction of initial magnetization direction,
    // positive_direction true points +, false in - direction
    af::array m = af::constant(0.0, this->n0, this->n1, this->n2, 3, f64);
    for (unsigned ix = 0; ix < this->n0; ix++) {
        for (unsigned iy = 0; iy < this->n1; iy++) {
            const double a = (double)(this->n0 / 2);
            const double b = (double)(this->n1 / 2);
            const double rx = double(ix) - this->n0 / 2.;
            const double ry = double(iy) - this->n1 / 2.;
            const double r = pow(rx, 2) / pow(a, 2) + pow(ry, 2) / pow(b, 2);
            if (r < 1) {
                for (unsigned iz = 0; iz < this->n2; iz++) {
                }
                if (positive_direction)
                    m(ix, iy, af::span, xyz) = 1;
                else
                    m(ix, iy, af::span, xyz) = -1;
            }
        }
    }
    std::cout << "Info: Mesh::ellipse(): n_cells should be approx a*b*M_PI*this->n2= "
              << this->n0 / 2 * this->n1 / 2 * M_PI * this->n2 << std::endl;
    return m;
}

af::array Mesh::init_vortex(const bool positive_direction) {
    // Returns an initial vortex magnetization
    // n_cells gives number of cells with non-zero Ms
    // positive_direction true, core points in +, false in - direction
    af::array m = af::constant(0.0, this->n0, this->n1, this->n2, 3, f64);
    for (unsigned ix = 0; ix < this->n0; ix++) {
        for (unsigned iy = 0; iy < this->n1; iy++) {
            const double rx = double(ix) - this->n0 / 2.;
            const double ry = double(iy) - this->n1 / 2.;
            const double r = sqrt(pow(rx, 2) + pow(ry, 2));
            if (r < this->n0 / 2.) {
                for (unsigned iz = 0; iz < this->n2; iz++) {
                }
                if (r == 0.) {
                    if (positive_direction)
                        m(ix, iy, af::span, 2) = 1;
                    else
                        m(ix, iy, af::span, 2) = -1;
                } else {
                    m(ix, iy, af::span, 0) = -ry / r;
                    m(ix, iy, af::span, 1) = rx / r;
                    if (positive_direction)
                        m(ix, iy, af::span, 2) = sqrt(this->n0) / r;
                    else
                        m(ix, iy, af::span, 2) = -sqrt(this->n0) / r;
                }
            }
        }
    }

    std::cout << "n_cells should be approx nx^2*M_PI/4.= " << pow(this->n0, 2) * M_PI / 4. << std::endl;
    m = normalize_handle_zero_vectors(m);
    return m;
}
af::array Mesh::init_sp4() {
    af::array m = af::constant(0.0, this->n0, this->n1, this->n2, 3, f64);
    m(af::seq(1, af::end - 1), af::span, af::span, 0) = 1;
    m(0, af::span, af::span, 1) = 1;
    m(-1, af::span, af::span, 1) = 1;
    return m;
}
} // namespace magnumafcpp
