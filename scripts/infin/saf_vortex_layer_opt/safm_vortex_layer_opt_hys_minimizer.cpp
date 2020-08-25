#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace magnumafcpp;

int main(int argc, char** argv) {
    std::cout << "argc" << argc << std::endl;
    for (int i = 0; i < argc; i++)
        std::cout << "Parameter " << i << " was " << argv[i] << "\n";
    std::string filepath(argc > 1 ? argv[1] : "../Data/Testing");
    if (argc > 0)
        filepath.append("/");
    std::cout << "Writing into path " << filepath << std::endl;
    af::setDevice(argc > 2 ? std::stoi(argv[2]) : 0);
    const double hzee_max = (argc > 3 ? std::stod(argv[3]) : 0.12); //[Tesla]
    const unsigned int steps_full_hysteresis = (argc > 4 ? std::stoi(argv[4]) : 200);
    std::string path_h_fl(argc > 5 ? argv[5] : filepath + "h_free_layer.vti"); // path to freelayer h vti

    af::info();
    std::cout.precision(24);

    af::array h_demag_safm;
    Mesh mesh(0, 0, 0, 0, 0, 0);
    vti_reader(h_demag_safm, mesh, path_h_fl);
    // Defining H_zee via lamdas
    const int im_free_layer = 0; // 0 == x, 2 == z
    auto zee_func = [h_demag_safm, hzee_max, steps_full_hysteresis](State state) -> af::array {
        double field_Tesla;
        if (state.steps < steps_full_hysteresis / 4) {
            field_Tesla = hzee_max * 4. * state.steps / steps_full_hysteresis;
        } else if (state.steps < 3 * steps_full_hysteresis / 4) {
            field_Tesla = -hzee_max * 4. * state.steps / steps_full_hysteresis + 2 * hzee_max;
        } else if (state.steps < steps_full_hysteresis * 5. / 4.) {
            field_Tesla = hzee_max * 4. * state.steps / steps_full_hysteresis - 4 * hzee_max;
        } else {
            field_Tesla = 0;
            std::cout << "WARNING ZEE time out of range, setting external "
                         "field to zero"
                      << std::endl;
        }
        // std::cout << "fild= "<< field_Tesla << std::endl;
        //{0.08715574, 0, 0.996194698};//field tile 5 degree sin\cos(5/360 * 2
        // pi)
        af::array zee = af::constant(0.0, state.mesh.n0, state.mesh.n1, state.mesh.n2, 3, f64);
        // const double hext_x_factor = std::sin(hext_angle / 360 * 2 * M_PI);
        // const double hext_z_factor = std::cos(hext_angle / 360 * 2 * M_PI);
        // zee(af::span, af::span, af::span, 0) =
        //    af::constant(hext_x_factor * field_Tesla / constants::mu0,
        //                 state.mesh.n0, state.mesh.n1, state.mesh.n2, 1, f64);
        zee(af::span, af::span, af::span, 0) =
            af::constant(field_Tesla / constants::mu0, state.mesh.n0, state.mesh.n1, state.mesh.n2, 1, f64);
        // zee(af::span, af::span, af::span, im_free_layer) =
        //    af::constant(hext_z_factor * field_Tesla / constants::mu0,
        //                 state.mesh.n0, state.mesh.n1, state.mesh.n2, 1, f64);
        // std::cout << "dims= " << zee.dims() << ", "<< h_demag_safm.dims() <<
        // std::endl;
        return zee + h_demag_safm;
    };

    // Parameter initialization
    ////Generating Objects
    const double Ms = 0.5 / constants::mu0; //[J/T/m^3] == [Joule/Tesla/meter^3] = 1.75 T/mu_0
    const double A = 1.5e-11;               //[J/m]
    const double Ku1 = 0.2e6;               // J/m^3

    State state(mesh, Ms, mesh.init_vortex(true));
    vti_writer_micro(state.Ms_field, mesh, filepath + "2nd_Ms");
    std::cout << "ncells= " << state.get_n_cells_() << std::endl;

    vti_writer_micro(state.m, mesh, (filepath + "2nd_minit_nonnormalized").c_str());
    state.m = normalize_handle_zero_vectors(state.m);
    vti_writer_micro(state.m, mesh, (filepath + "2nd_minit_renorm").c_str());

    af::timer timer_llgterms = af::timer::start();
    LBFGS_Minimizer minimizer = LBFGS_Minimizer(1e-6, 1000, 0);
    minimizer.of_convergence.open(filepath + "minimizer_convergence.dat");
    minimizer.llgterms_.push_back(LlgTerm(new DemagField(mesh)));
    minimizer.llgterms_.push_back(LlgTerm(new ExchangeField(A)));
    minimizer.llgterms_.push_back(LlgTerm(new UniaxialAnisotropyField(Ku1, {1, 0, 0})));
    minimizer.llgterms_.push_back(LlgTerm(new ExternalField(zee_func)));
    std::cout << "Llgterms assembled in " << af::timer::stop(timer_llgterms) << std::endl;

    std::ofstream stream;
    stream.precision(24);
    stream.open((filepath + "m.dat").c_str());
    stream << "# t	<mx>    <my>    <mz>    hzee" << std::endl;
    af::timer t_hys = af::timer::start();
    for (unsigned i = 0; i < steps_full_hysteresis * 5. / 4.; i++) {
        minimizer.Minimize(state);
        // state.calc_mean_m_steps(stream, constants::mu0 *
        // afvalue(minimizer.llgterms_[minimizer.llgterms_.size()-1]->h(state)(0,
        // 0, 0, im_free_layer)));
        state.calc_mean_m_steps(
            stream, constants::mu0 * minimizer.llgterms_[minimizer.llgterms_.size() - 1]->h(state)(0, 0, 0, af::span));
        if (state.steps % 10 == 0) {
            vti_writer_micro(state.m, mesh, (filepath + "m_hysteresis_" + std::to_string(state.steps)).c_str());
        }
        state.steps++;
        std::cout << "i=" << i << ", mean=" << state.meani(im_free_layer)
                  << ", h=" << constants::mu0 * afvalue(minimizer.llgterms_.back()->h(state)(0, 0, 0, im_free_layer))
                  << std::endl;
    }
    stream.close();
    std::cout << "time full hysteresis [af-s]: " << af::timer::stop(t_hys) << std::endl;
    return 0;
}
