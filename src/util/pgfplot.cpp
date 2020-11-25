#include "pgfplot.hpp"
#include "func.hpp"
#include <fstream>
#include <iostream>

namespace magnumafcpp {

void pgfplot_mz(af::array m, Mesh mesh, std::string outputfile, int mz_slice) {
    std::cout << "pgfplot_mz writing file " << outputfile << std::endl;

    std::ofstream stream;
    stream.precision(6);
    stream.open(outputfile.c_str());

    for (int y = 0; y < m.dims(1); y++) {
        for (int x = 0; x < m.dims(0); x++) {
            // if(x % 1 == 0 && y % 1 == 0){//NOTE: needed if file becomes too
            // big for latex
            stream << 1e9 * mesh.dx * x << " " << 1e9 * mesh.dy * y << " " << afvalue(m(x, y, mz_slice, 2))
                   << std::endl;
            std::cout << 1e9 * mesh.dx * x << " " << 1e9 * mesh.dy * y << " " << afvalue(m(x, y, mz_slice, 2))
                      << std::endl;
            //}
        }
        stream << std::endl;
        std::cout << std::endl;
    }
    stream.close();
}

void pgfplot_mz(af::array m, Mesh mesh, std::string outputfile, int mz_slice, int nmod) {
    std::cout << "pgfplot_mz writing file " << outputfile << std::endl;

    std::ofstream stream;
    stream.precision(6);
    stream.open(outputfile.c_str());

    for (int y = 0; y < m.dims(1); y++) {
        for (int x = 0; x < m.dims(0); x++) {
            if (x % nmod == nmod / 2 && y % nmod == nmod / 2) { // NOTE: needed if file becomes too big for latex
                stream << 1e9 * mesh.dx * x << " " << 1e9 * mesh.dy * y << " " << afvalue(m(x, y, mz_slice, 2))
                       << std::endl;
                std::cout << 1e9 * mesh.dx * x << " " << 1e9 * mesh.dy * y << " " << afvalue(m(x, y, mz_slice, 2))
                          << std::endl;
            }
        }
        stream << std::endl;
        std::cout << std::endl;
    }
    stream.close();
}

void pgfplot_mi(af::array m, Mesh mesh, std::string outputfile, int i, int mz_slice, int nmod) {
    std::cout << "pgfplot_mz writing file " << outputfile << std::endl;

    std::ofstream stream;
    stream.precision(6);
    stream.open(outputfile.c_str());

    for (int y = 0; y < m.dims(1); y++) {
        for (int x = 0; x < m.dims(0); x++) {
            if (x % nmod == nmod / 2 && y % nmod == nmod / 2) { // NOTE: needed if file becomes too big for latex
                stream << 1e9 * mesh.dx * x << " " << 1e9 * mesh.dy * y << " " << afvalue(m(x, y, mz_slice, i))
                       << std::endl;
                std::cout << 1e9 * mesh.dx * x << " " << 1e9 * mesh.dy * y << " " << afvalue(m(x, y, mz_slice, i))
                          << std::endl;
            }
        }
        stream << std::endl;
        std::cout << std::endl;
    }
    stream.close();
}

void pgfplot_nz_quiver(af::array m, Mesh mesh, std::string outputfile, int mz_slice, int nmod) {
    std::cout << "pgfplot_nz_quiver writing file " << outputfile << std::endl;
    std::ofstream stream;
    stream.precision(6);
    stream.open(outputfile.c_str());

    stream << "x"
           << " "
           << "y"
           << " "
           << "mx"
           << " "
           << "my" << std::endl;
    for (int y = 0; y < m.dims(1); y++) {
        for (int x = 0; x < m.dims(0); x++) {
            if (x % nmod == nmod / 2 && y % nmod == nmod / 2) { // NOTE: needed if file becomes too big for latex
                stream << 1e9 * mesh.dx * x << " " << 1e9 * mesh.dy * y << " " << afvalue(m(x, y, mz_slice, 0)) << " "
                       << afvalue(m(x, y, mz_slice, 1)) << std::endl;
                std::cout << 1e9 * mesh.dx * x << " " << 1e9 * mesh.dy * y << " " << afvalue(m(x, y, mz_slice, 0))
                          << " " << afvalue(m(x, y, mz_slice, 1)) << std::endl;
            }
        }
        stream << std::endl;
        std::cout << std::endl;
    }
    stream.close();
}

} // namespace magnumafcpp