// usage: ./convert_vti_to_ascii /path/to/vtifile.vti /path/to/newascii.txt

#include "arrayfire.h"
#include "magnum_af.hpp"
#include <string>

using namespace magnumafcpp;

int main(int argc, char** argv) {
    // Checking input variables and setting GPU Device
    for (int i = 0; i < argc; i++) {
        std::cout << "argv[" << i << "]: " << argv[i] << std::endl;
    }

    Mesh mesh(0, 0, 0, 0, 0, 0);
    af::array m;

    vti_reader(m, mesh, argv[1]);
    // vti_writer_micro(m(af::span, af::span, 0, af::span), Mesh(mesh.n0,
    // mesh.n1, 1, mesh.dx, mesh.dy, mesh.dz), filepath + "toplayer_" +
    // argv[3]); write_ascii(m(af::span, af::span, 0, af::span), Mesh(mesh.n0,
    // mesh.n1, 1, mesh.dx, mesh.dy, mesh.dz), filepath + "toplayer_" +
    // argv[4]);
    write_ascii(m, mesh, argv[2]);
    // read_ascii(m, mesh, filepath + argv[4]);

    return 0;
}
