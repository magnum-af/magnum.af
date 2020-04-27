#include "magnum_af.hpp"
#include "arrayfire.h"
#include <string>

using namespace magnumafcpp;

int main(int argc, char **argv)
{
    // Checking input variables and setting GPU Device
    for (int i = 0; i < argc; i++)
    {
        std::cout << "Parameter " << i << " was " << argv[i] << std::endl;
    }
    std::string filepath(argc > 1 ? argv[1] : "output_magnum.af/");
    af::setDevice(argc > 2 ? std::stoi(argv[2]) : 0);
    af::info();

    Mesh mesh;
    af::array m;
    
    vti_reader(m, mesh, filepath + argv[3]);
    //vti_writer_micro(m(af::span, af::span, 0, af::span), Mesh(mesh.n0, mesh.n1, 1, mesh.dx, mesh.dy, mesh.dz), filepath + "toplayer_" + argv[3]);
    //write_ascii(m(af::span, af::span, 0, af::span), Mesh(mesh.n0, mesh.n1, 1, mesh.dx, mesh.dy, mesh.dz), filepath + "toplayer_" + argv[4]);
    write_ascii(m, mesh, filepath + argv[4]);
    //read_ascii(m, mesh, filepath + argv[4]);

    return 0;
}
