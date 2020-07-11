#include "magnum_af.hpp"

using namespace magnumafcpp;

int main(int argc, char** argv) {
    // Checking input variables and setting GPU Device
    af::timer total_time = af::timer::start();
    for (int i = 0; i < argc; i++) {
        std::cout << "Parameter " << i << " was " << argv[i] << std::endl;
    }
    std::string filepath = std::string(argc > 1 ? argv[1] : "output_magnum.af/");
    af::setDevice(argc > 2 ? std::stoi(argv[2]) : 0);
    af::info();
    const int nx(argc > 3 ? std::stoi(argv[3]) : 128);
    const int ny(argc > 4 ? std::stoi(argv[4]) : 128);
    const double x(argc > 5 ? std::stod(argv[5]) : 1000e-9);
    const double y(argc > 6 ? std::stod(argv[6]) : 1000e-9);
    const int idz_vortex(argc > 7 ? std::stoi(argv[7]) : 8); // 1/3 vortex freelayer thickness
    std::cout << "read in values: " << nx << "\t" << ny << "\t" << x << "\t" << y << "\t" << std::endl;

    int i_callback = 0;
    auto calc_hz = [nx, ny, x, y, filepath, &i_callback, idz_vortex](double dz) -> double {
        // Parameter initialization
        double Ms = 0.5 / constants::mu0; // Vortex val//[J/T/m^3] == [Joule/Tesla/meter^3]

        const double vortex_z_height = 75e-9;              // vortex freelayer thickness
        const double dz_fl = vortex_z_height / idz_vortex; // cell-wise spacing in vortex freelayer
        std::vector<double> z_spacing = {dz, 0.5e-9, 10e-9, 1e-9};
        for (int i = 0; i < idz_vortex; i++) {
            z_spacing.push_back(dz_fl);
        }
        NonequispacedMesh mesh(nx, ny, x / nx, y / ny, z_spacing);
        // mesh.print(std::cout);

        // Initial magnetic field
        af::array m_free_geom_for_meandemag = af::constant(0, nx, ny, 1, 1, f64);
        af::array m = af::constant(0, mesh.dims, f64);

        for (int ix = 0; ix < nx; ix++) {
            for (int iy = 0; iy < ny; iy++) {
                const double a = (double)(nx / 2);
                const double b = (double)(ny / 2);
                const double rx = double(ix) - nx / 2.;
                const double ry = double(iy) - ny / 2.;
                const double r = pow(rx, 2) / pow(a, 2) + pow(ry, 2) / pow(b, 2);
                if (r < 1) {
                    // if(positive_direction) m(ix, iy, af::span, xyz)=1;
                    // else m(ix, iy, af::span, xyz)=-1;
                    m(ix, iy, 0, 0) = 1;  // SAFM Layer 0 in x
                    m(ix, iy, 2, 0) = -1; // SAFM Layer 1 in -x
                    // not here//m(ix, iy, 4, 0) = 1;// Free Layer in x
                    m_free_geom_for_meandemag(ix, iy, 0, 0) = 1;
                }
            }
        }

        // m(af::span, af::span, 0, 2) = 1;
        // m(af::span, af::span, 2, 2) = -1;
        State state(mesh, Ms, m, false);
        state.vtr_writer(filepath + "minit");

        NonEquiDemagField demag(mesh, false, false, 0);

        af::array h = demag.h(state);
        vtr_writer(h, mesh, filepath + "h");
        vtr_writer(h, mesh, filepath + "h_" + std::to_string(i_callback));
        // af::print("h slice", h(nx/2, ny/2, af::span, af::span));
        // af::print("h softmagnetic", h(nx/2, ny/2, 3, af::span));
        // mesh.print(stream);
        // const int index_free_layer = 4;
        Mesh mesh_fl(nx, ny, idz_vortex, x / nx, y / ny, dz_fl);
        auto fl_seq = af::seq(4, 4 + idz_vortex - 1);
        vti_writer_micro(h(af::span, af::span, fl_seq, af::span) *
                             af::tile(m_free_geom_for_meandemag, 1, 1, idz_vortex, 3),
                         mesh_fl, filepath + "h_free_layer");
        vti_writer_micro(h(af::span, af::span, fl_seq, af::span) *
                             af::tile(m_free_geom_for_meandemag, 1, 1, idz_vortex, 3),
                         mesh_fl, filepath + "h_free_layer_it_" + std::to_string(i_callback));

        std::ofstream stream;
        stream.precision(12);
        stream.open(filepath + "m.dat");

        stream << z_spacing[1] << ", " << afvalue(h(nx / 2, ny / 2, 5, 0)) << ", " << afvalue(h(nx / 2, ny / 2, 5, 1))
               << ", " << afvalue(h(nx / 2, ny / 2, 5, 2)) << std::endl;
        // std::cout << z_spacing[1] << ", " << afvalue(h(nx/2, ny/2, 3, 0)) <<
        // ", " << afvalue(h(nx/2, ny/2, 3, 1)) << ", " << afvalue(h(nx/2, ny/2,
        // 3, 2)) << std::endl; NOTE: previously used//return afvalue(h(nx/2,
        // ny/2, index_free_layer, 2)); std::cout << "test dims: " <<
        // af::mean(af::mean(h(af::span, af::span, index_free_layer, 2), 1),
        // 0).dims() << std::endl;
        stream.close();
        i_callback++;
        // return afvalue(
        //    af::mean(af::mean(h(af::span, af::span, index_free_layer, 0) *
        //                          m_free_geom_for_meandemag,
        //                      1),
        //             0));
        return afvalue(af::mean(af::mean(af::mean(af::mean(h(af::span, af::span, fl_seq, af::span) *
                                                               af::tile(m_free_geom_for_meandemag, 1, 1, idz_vortex, 3),
                                                           3),
                                                  2),
                                         1),
                                0));
    };

    std::cout.precision(18);
    // magnumaf::ZeroCrossing zc(calc_hz, 1e-6, 10, 9.9e-9, 10e-9, 10, 3);
    // auto result = zc.calc_x_and_f();
    NewtonIteration ni(calc_hz);
    auto result = ni.run(X0(10.001e-9), Precision(1e-8), EpsilonFactor(1e-10), Imax(100));
    std::cout << "result: x = " << result.first << ", f(x) = " << result.second << std::endl;

    std::cout << "total [af-s]: " << af::timer::stop(total_time) << std::endl;
    return 0;
}
