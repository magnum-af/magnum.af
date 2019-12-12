#include "magnum_af.hpp"

using namespace magnumafcpp;

int main(int argc, char **argv)
{
    // Checking input variables and setting GPU Device
    af::timer total_time = af::timer::start();
    for (int i = 0; i < argc; i++)
    {
        std::cout << "Parameter " << i << " was " << argv[i] << std::endl;
    }
    std::string filepath = std::string(argc > 1 ? argv[1] : "output_magnum.af/");
    af::setDevice(argc > 2 ? std::stoi(argv[2]) : 0);
    af::info();
    const int nx(argc > 3 ? std::stoi(argv[3]) : 128);
    const int ny(argc > 4 ? std::stoi(argv[4]) : 32);
    const double x(argc > 5 ? std::stod(argv[5]) : 100e-9);
    const double y(argc > 6 ? std::stod(argv[6]) : 25e-9);
    std::cout << "read in values: " << nx << "\t" << ny << "\t" << x << "\t" << y << "\t" << std::endl;

    int i_callback = 0;
    auto calc_hz = [nx, ny, x, y, filepath, &i_callback](double dz) -> double {
        // Parameter initialization
        double Ms = 1.393e6; //Vortex val//[J/T/m^3] == [Joule/Tesla/meter^3] = 1.75 T/mu_0
        //TODO//double A     = 1.5e-11;//Vortex val//[J/m]

        //Generating Objects
        const double dz_fl = 3e-9;
        std::vector<double> z_spacing = {5e-9, 1e-9, dz, 1e-9, dz_fl};
        NonequispacedMesh mesh(nx, ny, x / nx, y / ny, z_spacing);
        //mesh.print(std::cout);

        // Initial magnetic field
        af::array m = af::constant(0, mesh.dims, f64);

        for (int ix = 0; ix < nx; ix++)
        {
            for (int iy = 0; iy < ny; iy++)
            {
                const double a = (double)(nx / 2);
                const double b = (double)(ny / 2);
                const double rx = double(ix) - nx / 2.;
                const double ry = double(iy) - ny / 2.;
                const double r = pow(rx, 2) / pow(a, 2) + pow(ry, 2) / pow(b, 2);
                if (r < 1)
                {
                    for (int iz = 0; iz < mesh.nz; iz++)
                    {
                    }
                    //if(positive_direction) m(ix, iy, af::span, xyz)=1;
                    //else m(ix, iy, af::span, xyz)=-1;
                    m(ix, iy, 0, 2) = 1;  //SAFM Layer 0 in z
                    m(ix, iy, 2, 2) = -1; //SAFM Layer 1 in -z
                    //not here//m(ix, iy, 4, 0) = 1;// Free Layer in x
                }
            }
        }

        //m(af::span, af::span, 0, 2) = 1;
        //m(af::span, af::span, 2, 2) = -1;
        State state(mesh, Ms, m, false);
        state.vtr_writer(filepath + "minit");

        NonEquiDemagField demag(mesh, false, false, 0);

        af::array h = demag.h(state);
        vtr_writer(h, mesh, filepath + "h");
        vtr_writer(h, mesh, filepath + "h_" + std::to_string(i_callback));
        //af::print("h slice", h(nx/2, ny/2, af::span, af::span));
        //af::print("h softmagnetic", h(nx/2, ny/2, 3, af::span));
        //mesh.print(stream);
        const int index_free_layer = 4;
        Mesh mesh_fl(nx, ny, 1, x / nx, y / ny, dz_fl);
        vti_writer_micro(h(af::span, af::span, index_free_layer, af::span), mesh_fl, filepath + "h_free_layer");
        vti_writer_micro(h(af::span, af::span, index_free_layer, af::span), mesh_fl, filepath + "h_free_layer_it_" + std::to_string(i_callback));

        std::ofstream stream;
        stream.precision(12);
        stream.open(filepath + "m.dat");

        stream << z_spacing[1] << ", " << afvalue(h(nx / 2, ny / 2, index_free_layer, 0)) << ", " << afvalue(h(nx / 2, ny / 2, index_free_layer, 1)) << ", " << afvalue(h(nx / 2, ny / 2, index_free_layer, 2)) << std::endl;
        //std::cout << z_spacing[1] << ", " << afvalue(h(nx/2, ny/2, 3, 0)) << ", " << afvalue(h(nx/2, ny/2, 3, 1)) << ", " << afvalue(h(nx/2, ny/2, 3, 2)) << std::endl;
        //NOTE: previously used//return afvalue(h(nx/2, ny/2, index_free_layer, 2));
        //std::cout << "test dims: " << af::mean(af::mean(h(af::span, af::span, index_free_layer, 2), 1), 0).dims() << std::endl;
        stream.close();
        i_callback++;
        return afvalue(af::mean(af::mean(h(af::span, af::span, index_free_layer, 2), 1), 0));
    };

    std::cout.precision(18);
    //magnumaf::ZeroCrossing zc(calc_hz, 1e-6, 10, 9.9e-9, 10e-9, 10, 3);
    //auto result = zc.calc_x_and_f();
    NewtonIteration ni(calc_hz);
    auto result = ni.run(X0(5e-9), Precision(1e-8), EpsilonFactor(1e-10), Imax(100));
    std::cout << "result: x = " << result.first << ", f(x) = " << result.second << std::endl;

    std::cout << "total [af-s]: " << af::timer::stop(total_time) << std::endl;
    return 0;
}
