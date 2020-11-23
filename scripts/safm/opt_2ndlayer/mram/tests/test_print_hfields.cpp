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
    const int nx(argc > 3 ? std::stoi(argv[3]) : 64);
    const int ny(argc > 4 ? std::stoi(argv[4]) : 64);
    const double x(argc > 5 ? std::stod(argv[5]) : 60e-9);
    const double y(argc > 6 ? std::stod(argv[6]) : 60e-9);
    std::cout << "read in values: " << nx << "\t" << ny << "\t" << x << "\t" << y << "\t" << std::endl;

    // Parameter initialization
    int i_callback = 0;
    bool bottom_over_middle = true; // TODO use switch
    double dz = 5e-9;
    double Ms = 1.393e6; // Vortex val//[J/T/m^3] == [Joule/Tesla/meter^3] = 1.75 T/mu_0

    // Generating Objects
    const double dz_fl = 3e-9;
    std::vector<double> z_spacing = {5e-9, 1e-9, dz, 1e-9, dz_fl};
    NonequispacedMesh mesh(nx, ny, x / nx, y / ny, z_spacing);
    // std::cout << mesh << std::endl;

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
                if (bottom_over_middle) {
                    m(ix, iy, 0, 2) = 1; // SAFM Layer 0 in z
                } else {
                    m(ix, iy, 2, 2) = -1; // SAFM Layer 1 in -z
                }
                m_free_geom_for_meandemag(ix, iy, 0, 0) = 1;
            }
        }
    }

    State state(mesh, Ms, m, false);
    state.vtr_writer(filepath + "minit");

    NonEquiDemagField demag(mesh, false, false, 0);

    af::array h = demag.h(state);
    vtr_writer(h, mesh, filepath + "h");
    vtr_writer(h, mesh, filepath + "h_" + std::to_string(i_callback));
    const int index_free_layer = 4;
    Mesh mesh_fl(nx, ny, 1, x / nx, y / ny, dz_fl);
    vti_writer_micro(h(af::span, af::span, index_free_layer, af::span) *
                         af::tile(m_free_geom_for_meandemag, 1, 1, 1, 3),
                     mesh_fl, filepath + "h_free_layer");
    vti_writer_micro(h(af::span, af::span, index_free_layer, af::span) *
                         af::tile(m_free_geom_for_meandemag, 1, 1, 1, 3),
                     mesh_fl, filepath + "h_free_layer_it_" + std::to_string(i_callback));

    std::ofstream stream;
    stream.precision(12);
    stream.open(filepath + "m.dat");

    stream << z_spacing[1] << ", " << afvalue(h(nx / 2, ny / 2, index_free_layer, 0)) << ", "
           << afvalue(h(nx / 2, ny / 2, index_free_layer, 1)) << ", " << afvalue(h(nx / 2, ny / 2, index_free_layer, 2))
           << std::endl;
    stream.close();
    return 0;
}
