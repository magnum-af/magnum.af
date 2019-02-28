#include "arrayfire.h"
#include "magnum_af.hpp"

int main(int argc, char** argv)
{
    std::cout<<"argc"<<argc<<std::endl;
    for (int i=0; i<argc; i++) cout << "Parameter " << i << " was " << argv[i] << "\n";
    std::string filepath(argc>1? argv[1]: "../Data/Testing");
    if(argc>0)filepath.append("/");
    std::cout<<"Writing into path "<< filepath <<std::endl;
    setDevice(argc>2? std::stoi(argv[2]):0);
    const double hzee_max = (argc > 3 ? std::stod(argv[3]): 0.12); //[Tesla]
    const int steps_full_hysteresis =(argc > 4 ? std::stoi(argv[4]) : 200); 

    af::info();
    std::cout.precision(24);

    // Defining H_zee via lamdas
    auto zee_func= [ hzee_max, steps_full_hysteresis ] ( State state ) -> af::array {
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
        std::cout << "fild= "<< field_Tesla << std::endl;
        array zee = constant(0.0,state.mesh.n0,state.mesh.n1,state.mesh.n2,3,f64);
        zee(span,span,span,0)=constant(field_Tesla/state.material.mu0 ,state.mesh.n0,state.mesh.n1,state.mesh.n2,1,f64);
        return  zee;
    };

    // Parameter initialization
    const int nx = 250, ny=250 ,nz=1;
    const double x=1600e-9, y=1600e-9, z=65e-9;//[m] // Physical dimensions
  
    //Generating Objects
    Mesh mesh(nx,ny,nz,x/nx,y/ny,z/nz);
    Material material = Material();
    material.ms    = 1.393e6;//[J/T/m^3] == [Joule/Tesla/meter^3] = 1.75 T/mu_0
    material.A     = 1.5e-11;//[J/m]

    State state(mesh, material, mesh.init_vortex());
    vti_writer_micro(state.Ms, mesh ,(filepath + "Ms").c_str());
    std::cout << "ncells= "<< state.get_n_cells_() << std::endl;

    vti_writer_micro(state.m, mesh ,(filepath + "minit_nonnormalized").c_str());
    state.m=renormalize_handle_zero_values(state.m);
    vti_writer_micro(state.m, mesh ,(filepath + "minit_renorm").c_str());

    af::timer timer_llgterms = af::timer::start();
    LBFGS_Minimizer minimizer = LBFGS_Minimizer(1e-6, 1000, 0);
    minimizer.of_convergence.open(filepath + "minimizer_convergence.dat");
    minimizer.llgterms_.push_back( LlgTerm (new DemagField(mesh,material)));
    minimizer.llgterms_.push_back( LlgTerm (new ExchangeField(mesh,material)));
    minimizer.llgterms_.push_back( LlgTerm (new Zee(zee_func)));
    std::cout<<"Llgterms assembled in "<< af::timer::stop(timer_llgterms) <<std::endl;

    std::ofstream stream;
    stream.precision(24);
    stream.open ((filepath + "m.dat").c_str());
    stream << "# t	<mx>    <my>    <mz>    hzee" << std::endl;
    timer t_hys = af::timer::start();
    for (int i = 0; i < steps_full_hysteresis; i++){
        minimizer.Minimize(state);
        state.calc_mean_m_steps(stream, afvalue(minimizer.llgterms_[minimizer.llgterms_.size()-1]->h(state)(0,0,0,0)));
        if( state.steps % 10 == 0){
            vti_writer_micro(state.m, mesh ,(filepath + "m_hysteresis_"+std::to_string(state.steps)).c_str());
        }
        state.steps++;
    }
    stream.close();
    std::cout<<"time full hysteresis [af-s]: "<< af::timer::stop(t_hys) <<std::endl;
    return 0;
}
