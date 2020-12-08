// usage: ./convert_vti_to_ascii /path/to/vtifile.vti /path/to/newascii.txt

#include "util/ascii_io.hpp"
#include "util/vtk_IO.hpp"

using namespace magnumafcpp;

int main(int argc, char** argv) {
    // Checking input variables and setting GPU Device
    for (int i = 0; i < argc; i++) {
        std::cout << "argv[" << i << "]: " << argv[i] << std::endl;
    }

    Mesh mesh(0, 0, 0, 0, 0, 0);
    af::array m;

    vti_reader(m, mesh, argv[1]);
    // vti_writer_micro(m(af::span, af::span, 0, af::span), Mesh(mesh.nx,
    // mesh.ny, 1, mesh.dx, mesh.dy, mesh.dz), outdir / "toplayer_" +
    // posargs[0]); write_ascii(m(af::span, af::span, 0, af::span), Mesh(mesh.nx,
    // mesh.ny, 1, mesh.dx, mesh.dy, mesh.dz), outdir / "toplayer_" +
    // posargs[1]);
    write_ascii(m, mesh, argv[2]);
    // auto [arr, mesh] = read_ascii(filename, false);

    return 0;
}
