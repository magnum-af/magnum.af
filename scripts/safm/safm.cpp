#include "magnum_af.hpp"

std::string filepath;
std::ofstream stream;

double calc_hz(double dz){
    // Parameter initialization
    const double x=40.e-9, y=40e-9, z=4.e-9;
    const int nx = 16, ny=16 , nz=4;
    const double Ms = 8e5;

    //Generating Objects
    std::vector<double> z_spacing = {z/nz, dz, z/nz, z/nz};
    NonequispacedMesh mesh(nx, ny, x/nx, y/ny, z_spacing);
    //mesh.print(std::cout);

    // Initial magnetic field
    af::array m = af::constant(0, mesh.dims, f64);
    m(af::span, af::span, 0, 2) = 1;
    m(af::span, af::span, 1, 2) = -1;
    State state(mesh, Ms, m, false);
    vtr_writer(state.m, Mesh(nx, ny, nz, x/nx, y/ny, 0) , z_spacing, filepath + "minit");

    LlgTerms llgterm;
    NonEquiDemagField demag(mesh, false, false, 0);

    af::array h = demag.h(state);
    vtr_writer(h, Mesh(nx, ny, nz, x/nx, y/ny, 0) , z_spacing, filepath + "h");
    //af::print("h slice", h(nx/2, ny/2, af::span, af::span));
    //af::print("h softmagnetic", h(nx/2, ny/2, 3, af::span));
    //mesh.print(stream);
    stream << z_spacing[1] << ", " << afvalue(h(nx/2, ny/2, 3, 0)) << ", " << afvalue(h(nx/2, ny/2, 3, 1)) << ", " << afvalue(h(nx/2, ny/2, 3, 2)) << std::endl;
    //std::cout << z_spacing[1] << ", " << afvalue(h(nx/2, ny/2, 3, 0)) << ", " << afvalue(h(nx/2, ny/2, 3, 1)) << ", " << afvalue(h(nx/2, ny/2, 3, 2)) << std::endl;
    return afvalue(h(nx/2, ny/2, 3, 2));
}


double return_x(double x){return x;}

int main(int argc, char** argv)
{
    // Checking input variables and setting GPU Device
    af::timer total_time = af::timer::start();
    for (int i=0; i<argc; i++){cout << "Parameter " << i << " was " << argv[i] << std::endl;}
    filepath = std::string(argc>1? argv[1]: "output_magnum.af/");
    af::setDevice(argc>2? std::stoi(argv[2]):0);
    af::info();

    std::cout.precision(18);
    stream.precision(12);
    stream.open (filepath + "m.dat");

    //std::cout << zc.f(0.975e-9) << std::endl;
    //ZeroCrossing zc(return_x, 1e-12, 10, -1, 1.01, 101);
    //calc(filepath, stream);
    //ZeroCrossing zc(0, 1, 100);
    ZeroCrossing zc(calc_hz, 1e-6, 10, 0.9e-9, 1.0e-9);
    auto result = zc.calc_x_and_f();
    std::cout << "x = " << result.first << ", f(x) = " << result.second << std::endl;

    stream.close();
    std::cout<<"total [af-s]: "<< af::timer::stop(total_time) <<std::endl;
    return 0;
}
