#pragma once
#include "arrayfire.h"
#include "util/util.hpp" // for handle zero vals
#include "mesh.hpp"

namespace magnumafcpp::util {

inline af::array skyrmconf(const Mesh& mesh, const bool point_up = false) {
    // Returns a initial configuration to be relaxed into a skyrmion
    // if point_up is true, skyrmion centers points in +z, if false in -z
    af::array m = af::constant(0.0, mesh.nx, mesh.ny, mesh.nz, 3, f64);
    if (point_up) {
        m(af::span, af::span, af::span, 2) = 1.;
    } else {
        m(af::span, af::span, af::span, 2) = -1.;
    }
    for (unsigned ix = 0; ix < mesh.nx; ix++) {
        for (unsigned iy = 0; iy < mesh.ny; iy++) {
            const double rx = double(ix) - mesh.nx / 2.;
            const double ry = double(iy) - mesh.ny / 2.;
            const double r = sqrt(pow(rx, 2) + pow(ry, 2));
            if (r > mesh.nx / 4.) {
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
inline af::array ellipse(const Mesh& mesh, std::array<double, 3> vector, const bool verbose = true) {
    const double norm = std::sqrt(std::pow(vector[0], 2) + std::pow(vector[1], 2) + std::pow(vector[2], 2));
    vector[0] = vector[0] / norm;
    vector[1] = vector[1] / norm;
    vector[2] = vector[2] / norm;
    if (verbose)
        std::cout << "ellipse: norm=" << norm << std::endl;
    if (verbose)
        std::cout << "ellipse: normalized vector={" << vector[0] << ", " << vector[1] << ", " << vector[2] << "}"
                  << std::endl;

    long unsigned cells_within = 0;
    af::array m = af::constant(0.0, mesh.nx, mesh.ny, mesh.nz, 3, f64);
    for (unsigned ix = 0; ix < mesh.nx; ix++) {
        for (unsigned iy = 0; iy < mesh.ny; iy++) {
            const double a = (double)(mesh.nx / 2);
            const double b = (double)(mesh.ny / 2);
            const double rx = double(ix) - mesh.nx / 2.;
            const double ry = double(iy) - mesh.ny / 2.;
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
        std::cout << "Info: ellipse(): cells within cylinder = " << cells_within
                  << ", which should be approx a*b*M_PI*nz = " << mesh.nx / 2 * mesh.ny / 2 * M_PI * mesh.nz
                  << std::endl;
    return m;
}

inline af::array ellipse(const Mesh& mesh, const unsigned xyz = 0, const bool positive_direction = true) {
    // Returns an initial elliptical magnetization
    // n_cells gives number of cells with non-zero Ms
    // xyz gives direction of initial magnetization direction,
    // positive_direction true points +, false in - direction
    af::array m = af::constant(0.0, mesh.nx, mesh.ny, mesh.nz, 3, f64);
    for (unsigned ix = 0; ix < mesh.nx; ix++) {
        for (unsigned iy = 0; iy < mesh.ny; iy++) {
            const double a = (double)(mesh.nx / 2);
            const double b = (double)(mesh.ny / 2);
            const double rx = double(ix) - mesh.nx / 2.;
            const double ry = double(iy) - mesh.ny / 2.;
            const double r = pow(rx, 2) / pow(a, 2) + pow(ry, 2) / pow(b, 2);
            if (r < 1) {
                for (unsigned iz = 0; iz < mesh.nz; iz++) {
                }
                if (positive_direction)
                    m(ix, iy, af::span, xyz) = 1;
                else
                    m(ix, iy, af::span, xyz) = -1;
            }
        }
    }
    std::cout << "Info: ellipse(): n_cells should be approx a*b*M_PI*mesh.nz= "
              << mesh.nx / 2 * mesh.ny / 2 * M_PI * mesh.nz << std::endl;
    return m;
}

inline af::array init_vortex(const Mesh& mesh, const bool positive_direction = true) {
    // Returns an initial vortex magnetization
    // n_cells gives number of cells with non-zero Ms
    // positive_direction true, core points in +, false in - direction
    af::array m = af::constant(0.0, mesh.nx, mesh.ny, mesh.nz, 3, f64);
    for (unsigned ix = 0; ix < mesh.nx; ix++) {
        for (unsigned iy = 0; iy < mesh.ny; iy++) {
            const double rx = double(ix) - mesh.nx / 2.;
            const double ry = double(iy) - mesh.ny / 2.;
            const double r = sqrt(pow(rx, 2) + pow(ry, 2));
            if (r < mesh.nx / 2.) {
                for (unsigned iz = 0; iz < mesh.nz; iz++) {
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
                        m(ix, iy, af::span, 2) = sqrt(mesh.nx) / r;
                    else
                        m(ix, iy, af::span, 2) = -sqrt(mesh.nx) / r;
                }
            }
        }
    }

    std::cout << "n_cells should be approx nx^2*M_PI/4.= " << pow(mesh.nx, 2) * M_PI / 4. << std::endl;
    m = normalize_handle_zero_vectors(m);
    return m;
}
inline af::array init_sp4(const Mesh& mesh) {
    af::array m = af::constant(0.0, mesh.nx, mesh.ny, mesh.nz, 3, f64);
    m(af::seq(1, af::end - 1), af::span, af::span, 0) = 1;
    m(0, af::span, af::span, 1) = 1;
    m(-1, af::span, af::span, 1) = 1;
    return m;
}

} // namespace magnumafcpp::util
