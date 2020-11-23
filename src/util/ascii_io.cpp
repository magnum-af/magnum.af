#include "ascii_io.hpp"
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace magnumafcpp {

unsigned colum_major_stride(const unsigned i, const unsigned j, const unsigned k, const unsigned l, const unsigned ni,
                            const unsigned nj, const unsigned nk) {
    return i + ni * (j + nj * (k + nk * l));
}

void write_ascii(const af::array& a, const Mesh& mesh, std::string filename, bool verbose, int precision) {
    const unsigned nx = a.dims(0);
    const unsigned ny = a.dims(1);
    const unsigned nz = a.dims(2);
    const unsigned n_scalar = a.dims(3);

    std::ofstream stream(filename);
    if (stream.is_open()) {
        af::timer timer = af::timer::start();
        stream.precision(precision);
        int val_setw = precision + 4; // formatted at least for 0.*** values
        if (verbose) {
            printf("writing to ascii file %s with precision of %d\n", filename.c_str(), precision);
        }
        stream << "# " << a.dims(0) << " " << a.dims(1) << " " << a.dims(2) << " " << a.dims(3) << " " << mesh.dx << " "
               << mesh.dy << " " << mesh.dz << std::endl;
        stream << "# magnum.af ascii-encoded " << a.dims(3) << "d spacial field" << std::endl;
        stream << "# plot with gnuplot using: set view equal xyz; splot '" << filename
               << "' using 1:2:3:4:5:6 with vectors" << std::endl;
        stream << "# color code mx component: set view equal xyz; splot '" << filename
               << "' using ($1*1e9):($2*1e9):($3*1e9):4:5:6:4 with vectors "
                  "linecolor palette"
               << std::endl;
        // stream << "#                        : set view equal xyz; splot '" <<
        // filename << "' using ($1*1e9):($2*1e9):($3*1e9):4:5:6 with vectors"
        // << std::endl; stream << "#                        : set view equal
        // xyz; splot '" << filename << "' using
        // 1:2:3:4:5:6:(sqrt($4**2+$5**2+$6**2)) with vectors linecolor palette
        // z" << std::endl;
        stream << "#\n";
        stream << std::setw(val_setw) << std::left << "# x-position[m]"
               << "\t";
        stream << std::setw(val_setw) << std::left << "y-position[m]"
               << "\t";
        stream << std::setw(val_setw) << std::left << "z-position[m]"
               << "\t";
        if (n_scalar == 3) {
            stream << std::setw(val_setw) << std::left << "mx[a.u.]"
                   << "\t";
            stream << std::setw(val_setw) << std::left << "my[a.u.]"
                   << "\t";
            stream << std::setw(val_setw) << std::left << "mz[a.u.]"
                   << "\n";
        } else {
            for (unsigned i_scalar = 0; i_scalar < n_scalar; i_scalar++) {
                stream << std::setw(val_setw) << std::left << std::to_string(i_scalar) + "-component";
                if (i_scalar != (n_scalar - 1)) {
                    stream << "\t";
                }
            }
            stream << "\n";
        }

        // copying raw data to host is factor ~10 faster than accessing
        // af::array with indices
        double* host_a = a.host<double>();

        for (unsigned ix = 0; ix < nx; ix++) {
            for (unsigned iy = 0; iy < ny; iy++) {
                for (unsigned iz = 0; iz < nz; iz++) {
                    stream << std::setw(val_setw) << std::left << ix * mesh.dx << "\t";
                    stream << std::setw(val_setw) << std::left << iy * mesh.dy << "\t";
                    stream << std::setw(val_setw) << std::left << iz * mesh.dz << "\t";

                    for (unsigned i_scalar = 0; i_scalar < n_scalar; i_scalar++) {
                        stream << std::setw(val_setw) << std::left
                               << host_a[colum_major_stride(ix, iy, iz, i_scalar, nx, ny, nz)];
                        if (i_scalar != (n_scalar - 1)) {
                            stream << "\t";
                        }
                    }
                    stream << "\n";
                }
            }
        }
        af::freeHost(host_a);
        if (verbose) {
            printf("Wrote ascii file in %f [s]\n", timer.stop());
        }
    }
}
std::pair<af::array, Mesh> read_ascii(std::string filename, bool verbose) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        printf("read_ascii: Could not read file! Aborting...\n");
        std::exit(EXIT_FAILURE);
    }
    af::timer timer = af::timer::start();
    if (verbose) {
        printf("reading from ascii file %s\n", filename.c_str());
    }

    std::string line;

    // Getting mesh size and discretization
    unsigned nx, ny, nz, n_scalar;
    double dx, dy, dz;
    std::getline(infile, line);
    line.erase(0, 2);
    std::istringstream iss1(line);
    if (!(iss1 >> nx >> ny >> nz >> n_scalar >> dx >> dy >> dz)) {
        printf("Error while reading line!\n");
    }
    std::vector<double> raw_read_a(nx * ny * nz * n_scalar);

    // skipping all further lines starting with '#'
    while (std::getline(infile, line)) {
        if (line.at(0) == '#') {
            continue;
        }

        // reading in with loop
        for (unsigned ix = 0; ix < nx; ix++) {
            for (unsigned iy = 0; iy < ny; iy++) {
                for (unsigned iz = 0; iz < nz; iz++) {
                    std::istringstream iss(line);
                    double ddx, ddy, ddz;
                    if (!(iss >> ddx >> ddy >> ddz)) {
                        printf("Error while reading line!\n");
                        break;
                    }
                    for (unsigned i_scalar = 0; i_scalar < n_scalar; i_scalar++) {
                        double value;
                        if (!(iss >> value)) {
                            printf("Error while reading line!\n");
                            break;
                        }
                        raw_read_a[colum_major_stride(ix, iy, iz, i_scalar, nx, ny, nz)] = value;
                    }
                    std::getline(infile, line);
                }
            }
        }
    }
    if (verbose) {
        printf("Wrote ascii file in %f [s]\n", timer.stop());
    }
    return {af::array(nx, ny, nz, n_scalar, raw_read_a.data()), Mesh(nx, ny, nz, dx, dy, dz)};
}

} // namespace magnumafcpp
