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
    const double hzee_max = (argc > 3 ? std::stod(argv[3]) : 0.12);            //[Tesla]
    const double integr_time = (argc > 4 ? std::stod(argv[4]) : 4 * 50e-9);    //[s]
    std::string path_h_fl(argc > 5 ? argv[5] : filepath + "h_free_layer.vti"); // path to freelayer h vti

    af::info();
    // std::cout.precision(24);

    af::array h_demag_safm;
    Mesh mesh;
    vti_reader(h_demag_safm, mesh, path_h_fl);
    // Defining H_zee via lamdas
    double integr_time_per_quater = integr_time / 4.;
    double rate = hzee_max / integr_time_per_quater; //[T/s]
    std::cout << "hzee_max= " << hzee_max << ", rate=" << rate << ", integr_time_per_quater=" << integr_time_per_quater
              << std::endl;
    auto zee_func = [h_demag_safm, hzee_max, rate](State state) -> af::array {
        double field_Tesla = 0;
        if (state.t < hzee_max / rate)
            field_Tesla = rate * state.t;
        else if (state.t < 3 * hzee_max / rate)
            field_Tesla = -rate * state.t + 2 * hzee_max;
        else if (state.t < 5 * hzee_max / rate)
            field_Tesla = rate * state.t - 4 * hzee_max;
        else {
            field_Tesla = rate * state.t - 4 * hzee_max;
            std::cout << "NOTE: zee time out of range" << std::endl;
        }
        // double field_Tesla;
        // if(state.steps < steps_full_hysteresis/4){
        //    field_Tesla = hzee_max * 4. * state.steps/steps_full_hysteresis;
        //}
        // else if (state.steps < 3*steps_full_hysteresis/4){
        //    field_Tesla = - hzee_max * 4. * state.steps/steps_full_hysteresis
        //    + 2*hzee_max;
        //}
        // else if(state.steps < steps_full_hysteresis){
        //    field_Tesla = hzee_max * 4. *state.steps/steps_full_hysteresis -
        //    4*hzee_max;
        //}
        // else {
        //    field_Tesla = 0; std::cout << "WARNING ZEE time out of range,
        //    setting external field to zero" << std::endl;
        //}
        // std::cout << "fild= "<< field_Tesla << std::endl;
        af::array zee = af::constant(0.0, state.mesh.nx, state.mesh.ny, state.mesh.nz, 3, f64);
        zee(af::span, af::span, af::span, 0) =
            af::constant(field_Tesla / constants::mu0, state.mesh.nx, state.mesh.ny, state.mesh.nz, 1, f64);
        // std::cout << "dims= " << zee.dims() << ", "<< h_demag_safm.dims() <<
        // std::endl;
        return zee + h_demag_safm;
    };

    // Parameter initialization
    ////Generating Objects
    double Ms = 1.393e6; //[J/T/m^3] == [Joule/Tesla/meter^3] = 1.75 T/mu_0
    double A = 1.5e-11;  //[J/m]

    State state(mesh, Ms, util::ellipse(mesh, 1));
    vti_writer_micro(state.Ms_field, mesh, filepath + "2nd_Ms");
    std::cout << "ncells= " << state.get_n_cells_() << std::endl;

    vti_writer_micro(state.m, mesh, (filepath + "2nd_minit_nonnormalized").c_str());
    state.m = normalize_handle_zero_vectors(state.m);
    vti_writer_micro(state.m, mesh, (filepath + "2nd_minit_renorm").c_str());

    af::timer timer_llgterms = af::timer::start();
    auto demag = uptr_Fieldterm(new DemagField(mesh));
    auto exch = uptr_Fieldterm(new ExchangeField(A));
    auto zee = uptr_Fieldterm(new ExternalField(zee_func));
    vec_uptr_Fieldterm llgterms = {demag, exch, zee};
    LLGIntegrator llg(1, llgterms);
    // LBFGS_Minimizer minimizer = LBFGS_Minimizer(1e-6, 1000, 0);
    // minimizer.of_convergence.open(filepath + "minimizer_convergence.dat");
    // minimizer.llgterms_.push_back( uptr_Fieldterm (new DemagField(mesh)));
    // minimizer.llgterms_.push_back( uptr_Fieldterm (new ExchangeField(A)));
    // minimizer.llgterms_.push_back( uptr_Fieldterm (new ExternalField(zee_func)));
    std::cout << "Llgterms assembled in " << af::timer::stop(timer_llgterms) << std::endl;

    std::ofstream stream;
    stream.precision(24);
    stream.open((filepath + "m.dat").c_str());
    stream << "# t	<mx>    <my>    <mz>    hzee" << std::endl;
    af::timer t_hys = af::timer::start();
    while (state.t < 5 * hzee_max / rate) {
        if (state.steps % 100 == 0)
            std::cout << "%=" << state.t / integr_time * 4. / 5. << ", i=" << state.steps << ", t=" << state.t
                      << ", <m>=" << state.meani(0)
                      << ", hzee=" << constants::mu0 * afvalue(llg.llgterms.back()->h(state)(0, 0, 0, 0)) << std::endl;
        llg.step(state);
        state.calc_mean_m_steps(stream, constants::mu0 * afvalue(llg.llgterms.back()->h(state)(0, 0, 0, 0)));
        if (state.steps % 1000 == 0) {
            vti_writer_micro(state.m, mesh, (filepath + "m_hysteresis_" + std::to_string(state.steps)).c_str());
        }
    }
    // for (unsigned i = 0; i < steps_full_hysteresis; i++){
    //    minimizer.Minimize(state);
    //    state.calc_mean_m_steps(stream,
    //    afvalue(minimizer.llgterms_[minimizer.llgterms_.size()-1]->h(state)(0,
    //    0, 0, 0))); if( state.steps % 10 == 0){
    //        vti_writer_micro(state.m, mesh , (filepath +
    //        "m_hysteresis_"+std::to_string(state.steps)).c_str());
    //    }
    //    state.steps++;
    //    std::cout << "i=" << i << std::endl;
    stream.close();
    std::cout << "time full hysteresis [af-s]: " << af::timer::stop(t_hys) << std::endl;
    return 0;
}
