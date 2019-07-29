#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace magnumaf;


int main(int argc, char** argv)
{
    std::cout<<"argc"<<argc<<std::endl;
    for (int i=0; i<argc; i++) std::cout << "Parameter " << i << " was " << argv[i] << "\n";
    std::string filepath(argc>1? argv[1]: "../Data/Testing");
    if(argc>0)filepath.append("/");
    std::cout<<"Writing into path "<< filepath <<std::endl;
    af::setDevice(argc>2? std::stoi(argv[2]):0);
    std::string path_h_fl(argc>3? argv[3]: filepath + "h_free_layer.vti"); // path to freelayer h vti
    const double hzee_max = (argc > 4 ? std::stod(argv[4]): 0.12); //[Tesla]
    const unsigned int steps_full_hysteresis =(argc > 5 ? std::stoi(argv[5]) : 200);

    af::info();
    std::cout.precision(24);

    af::array h_demag_safm;
    Mesh mesh;
    vti_reader(h_demag_safm, mesh, path_h_fl);
    // Defining H_zee via lamdas
    auto zee_func= [h_demag_safm, hzee_max, steps_full_hysteresis ] ( State state ) -> af::array {
        double field_Tesla;
        if(state.steps < steps_full_hysteresis/4){
            field_Tesla = hzee_max * 4. * state.steps/steps_full_hysteresis;
        }
        else if (state.steps < 3*steps_full_hysteresis/4){
            field_Tesla = - hzee_max * 4. * state.steps/steps_full_hysteresis + 2*hzee_max;
        }
        else if(state.steps < steps_full_hysteresis){
            field_Tesla = hzee_max * 4. *state.steps/steps_full_hysteresis - 4*hzee_max;
        }
        else {
            field_Tesla = 0; std::cout << "WARNING ZEE time out of range, setting external field to zero" << std::endl;
        }
        //std::cout << "fild= "<< field_Tesla << std::endl;
        af:: array zee = af::constant(0.0, state.mesh.n0, state.mesh.n1, state.mesh.n2, 3, f64);
        zee(af::span, af::span, af::span, 0) = af::constant(field_Tesla/constants::mu0 , state.mesh.n0, state.mesh.n1, state.mesh.n2, 1, f64);
        //std::cout << "dims= " << zee.dims() << ", "<< h_demag_safm.dims() << std::endl;
        return  zee + h_demag_safm;
    };

    // Parameter initialization
    ////Generating Objects
    double Ms    = 1.393e6;//[J/T/m^3] == [Joule/Tesla/meter^3] = 1.75 T/mu_0
    double A     = 1.5e-11;//[J/m]

    State state(mesh, Ms, mesh.ellipse());
    vti_writer_micro(state.Ms_field, mesh, filepath + "2nd_Ms");
    std::cout << "ncells= "<< state.get_n_cells_() << std::endl;

    vti_writer_micro(state.m, mesh , (filepath + "2nd_minit_nonnormalized").c_str());
    state.m=renormalize_handle_zero_values(state.m);
    vti_writer_micro(state.m, mesh , (filepath + "2nd_minit_renorm").c_str());

    af::timer timer_llgterms = af::timer::start();
    LBFGS_Minimizer minimizer = LBFGS_Minimizer(1e-6, 1000, 0);
    minimizer.of_convergence.open(filepath + "minimizer_convergence.dat");
    minimizer.llgterms_.push_back( LlgTerm (new DemagField(mesh)));
    minimizer.llgterms_.push_back( LlgTerm (new ExchangeField(A)));
    minimizer.llgterms_.push_back( LlgTerm (new ExternalField(zee_func)));
    std::cout<<"Llgterms assembled in "<< af::timer::stop(timer_llgterms) <<std::endl;

    std::ofstream stream;
    stream.precision(24);
    stream.open ((filepath + "m.dat").c_str());
    stream << "# t	<mx>    <my>    <mz>    hzee" << std::endl;
    af::timer t_hys = af::timer::start();
    for (unsigned i = 0; i < steps_full_hysteresis; i++){
        minimizer.Minimize(state);
        state.calc_mean_m_steps(stream, afvalue(minimizer.llgterms_[minimizer.llgterms_.size()-1]->h(state)(0, 0, 0, 0)));
        if( state.steps % 10 == 0){
            vti_writer_micro(state.m, mesh , (filepath + "m_hysteresis_"+std::to_string(state.steps)).c_str());
        }
        state.steps++;
        std::cout << "i=" << i << std::endl;
    }
    stream.close();
    std::cout<<"time full hysteresis [af-s]: "<< af::timer::stop(t_hys) <<std::endl;
    return 0;
}
