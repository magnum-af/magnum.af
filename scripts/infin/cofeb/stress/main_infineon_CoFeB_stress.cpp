#include "arrayfire.h"
#include "magnum_af.hpp"

int main(int argc, char** argv) {
    std::string filepath(argc > 1 ? argv[1] : "../Data/Testing");
    if (argc > 1) {
        filepath.append("/");
    }
    af::setDevice(argc > 2 ? std::stoi(argv[2]) : 0);
    const double A = double(argc > 3 ? std::stod(argv[3]) * 1e-3 / (4e-7 * M_PI)
                                     : (double)(0.05 / (4e-7 * M_PI)));     // Input a in mT, argv[3]=25 mT is
                                                                            // converted to 0.025 T and divided by mu0
    const double B = double(argc > 4 ? std::stod(argv[4]) / 100 : 1.0) * A; // Input a in percent, B=1.0 == 100%
    std::cout << "Writing into path " << filepath.c_str() << std::endl;
    std::cout << "A=" << A << "B= " << B << std::endl;
    std::cout.precision(24);
    af::info();

    // Defining lamdas
    auto zee_func = [A, B](State state) -> af::array {
        double phi = 2 * M_PI * (state.t);
        af::dim4 dim = af::dim4(state.mesh.nx, state.mesh.ny, state.mesh.nz, 1);
        af::array zee = af::array(dims_vector(state.mesh), f64);
        zee(af::span, af::span, af::span, 0) = constant(A * std::cos(phi), dim, f64);
        zee(af::span, af::span, af::span, 1) = constant(A * std::sin(phi), dim, f64);
        zee(af::span, af::span, af::span, 2) = constant(0.0, dim, f64);
        // zee(af::span, af::span, af::span, 2)=constant( A * std::sin(phi) ,
        // state.mesh.nx, state.mesh.ny, state.mesh.nz, 1, f64);
        return zee;
    };

    // Parameter initialization
    const double x = 800e-9, y = 800e-9,
                 z = 1.3e-3 / 1.056e6; //[m] // z for 100mT lin range t_CoFeB
                                       //= 1.3e-3/1.056e6
    const int nx = 250, ny = 250, nz = 1;

    // Generating Objects
    Mesh mesh(nx, ny, nz, x / nx, y / ny, z / nz);
    std::cout << mesh << std::endl;
    Material material = Material();
    state.Ms = 1.58 / constants::mu0; // [J/T/m^3] = Ms = Js/mu0 = 1.58 Tesla
                                      // /mu_0 // Js = 1.58 Tesla
    double A_ex = 15e-12;             // [J/m]
    double Ku1 = 1.3e-3 / z;          // [J/m^3] // Ku1 = K_total - K_shape = Hk*Js/2/mu0 + Js^2/2/mu0 = |
                                      // [Hk and Js in Tesla] | = ((0.1*1.58)/2/(4*pi*1e-7) +
                                      // (1.58)^2/(2)/(4*pi*1e-7)) = 1.056e6

    double Ku1_stress = 1400; // TODO

    State state(mesh, material, util::ellipse(mesh, 2));
    std::cout << "ncells= " << state.get_n_cells_() << std::endl;
    std::cout << "Mean i for check: " << state.meani(0) << "\t" << state.meani(1) << "\t" << state.meani(2)
              << std::endl;

    vti_writer_micro(state.Ms, mesh, (filepath + "Ms").c_str());
    vti_writer_micro(state.m, mesh, (filepath + "minit").c_str());
    std::cout << mesh << std::endl;

    af::timer timer_llgterms = af::timer::start();
    LBFGS_Minimizer minimizer;
    // LBFGS_Minimizer minimizer = LBFGS_Minimizer();//Fails on GTO, maybe due
    // to gcc verions < 7.2
    minimizer.llgterms_.push_back(uptr_Fieldterm(new DemagField(mesh)));
    minimizer.llgterms_.push_back(uptr_Fieldterm(new ExchangeField(A_ex)));
    minimizer.llgterms_.push_back(uptr_Fieldterm(new UniaxialAnisotropyField(Ku1)));
    minimizer.llgterms_.push_back(uptr_Fieldterm(new UniaxialAnisotropyField(Ku1_stress, {1, 0, 0})));
    minimizer.llgterms_.push_back(uptr_Fieldterm(new ExternalField(zee_func)));
    std::cout << "Llgterms assembled in " << af::timer::stop(timer_llgterms) << std::endl;

    // Checking aniso
    vti_writer_micro(minimizer.llgterms_.end()[-3]->h(state), mesh, filepath + "check_h_ani_z");
    vti_writer_micro(minimizer.llgterms_.end()[-2]->h(state), mesh,
                     filepath + "check_h_ani_stress"); // TODO this looks strange
    vti_writer_micro(minimizer.llgterms_.end()[-1]->h(state), mesh, filepath + "check_h_zee");

    // obtaining relaxed magnetization
    af::timer t = af::timer::start();
    minimizer.Minimize(state);
    std::cout << "timerelax [af-s]: " << af::timer::stop(t) << std::endl;
    vti_writer_micro(state.m, mesh, filepath + "mrelax");

    // Starting elliptical field loop
    std::ofstream stream;
    stream.precision(12);
    stream.open((filepath + "m.dat").c_str());
    stream << "# id-val  	<mx>    <my>    <mz>    H_x    H_y    H_z" << std::endl;
    af::timer t_hys = af::timer::start();
    const int max_i = 100;
    for (int i = 0; i <= max_i; i++) {
        state.t = (double)i / (double)max_i;
        state.steps++;
        minimizer.Minimize(state);
        state.calc_mean_m(stream, minimizer.llgterms_.end()[-1]->h(state)(0, 0, 0, af::span));
        vti_writer_micro(state.m, mesh, filepath + "m_" + std::to_string(state.steps));
        vti_writer_micro(minimizer.llgterms_.end()[-2]->h(state), mesh,
                         filepath + "check_h_ani_stress" +
                             std::to_string(state.steps)); // TODO this looks interesting, value drops
                                                           // at the boundaries of the disc
    }
    stream.close();
    std::cout << "time full hysteresis [af-s]: " << af::timer::stop(t_hys) << std::endl;
    return 0;
}
