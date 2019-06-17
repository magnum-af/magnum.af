//on GTO
//make -j && ../bin/magnum.af-opencl /home/paul/git/magnum.af/Data/Testing 0 0.01296 7200000
#include "arrayfire.h"
#include "magnum_af.hpp"
using namespace af; typedef std::shared_ptr<LLGTerm> llgt_ptr;
int main(int argc, char** argv)
{
    if(argc>1) setDevice(std::stoi(argv[2]));
    info();

    std::cout<<"argc"<<argc<<std::endl;
     for (int i=0; i<argc; i++)
          cout << "Parameter " << i << " was " << argv[i] << "\n";
    // Parameter initialization
    const int nx = 112, ny=112 ,nz=1;//nz=5 -> lz=(5-1)*dx
    const double dx=2.715e-10;

    double n_interp = 60;
    double string_dt=5e-14;
    const int string_steps = 10000;
    double string_abort_rel_diff = 1e-12;
    double string_abort_abs_diff = 1e-27;

    std::string filepath(argc>0? argv[1]: "../Data/Testing/");
    if(argc>0)filepath.append("/");
    std::cout<<"Writing into path "<<filepath.c_str()<<std::endl;


    //Generating Objects
    Mesh mesh(nx,ny,nz,dx,dx,dx);
    Material material = Material();
    state.Ms    = 1.1e6;
    material.A     = 1.6e-11;
    material.alpha = 1;
    state.material.afsync  = false;

    //material.D=2*5.76e-3;
    material.D=std::stod(argv[3]);
    std::cout<<"D="<<material.D<<std::endl;
    //material.D_axis[2]=-1;

    //material.Ku1=6.4e6;
    material.Ku1=std::stod(argv[4]);
    std::cout<<"Ku1="<<material.Ku1<<std::endl;
    //material.Ku1_axis[2]=1;

    material.J_atom=2.*material.A*dx;
    material.D_atom= material.D * pow(dx,2);
    //old values in prev versions, this is now wrong in pthmag:
    //const double wrong_J_atom=4.*material.A*dx;
    const double wrong_D_atom= 2.* material.D * pow(dx,2);
    std::cout<<"D_atom="<<material.D_atom<<std::endl;
    material.K_atom=material.Ku1*pow(dx,3);
    std::cout<<"Ku1_atom="<<material.K_atom<<std::endl;
    material.p=state.Ms*pow(dx,3);//Compensate nz=1 instead of nz=4


    array m;
    Mesh testmesh(nx,ny,nz,dx,dx,dx);

    State state(mesh,material, m);
    //vti_writer_atom(state.m, mesh ,(filepath + "minit").c_str());

    std::vector<llgt_ptr> llgterm;
    llgterm.push_back( llgt_ptr (new AtomisticDipoleDipoleField(mesh)));
    llgterm.push_back( llgt_ptr (new AtomisticExchangeField(mesh)));
    llgterm.push_back( llgt_ptr (new AtomisticDmiField(mesh,material)));
    llgterm.push_back( llgt_ptr (new AtomisticUniaxialAnisotropyField(mesh,material)));

    LLG Llg(state,llgterm);

    std::vector<State> inputimages;
    for(int i=0; i < n_interp; i++){
        vti_reader(m, testmesh, std::string("/home/paul/git/magnum.af/Data/phasespace/1n/hotfix_sign/d")+argv[3]+"/k"+argv[4]+"/skyrm_image"+std::to_string(i)+".vti");
        inputimages.push_back(State(mesh, material, m));
    }

    String string(state,inputimages, n_interp, string_dt ,llgterm);
    string.write_vti(filepath+"init_string");
    string.calc_E();
    for(unsigned int i=0; i<string.images.size(); i++){
       std::cout << i << "   " << string.E[i] << std::endl;
    }
    auto max = std::max_element(std::begin(string.E), std::end(string.E));
    std::cout << *max << std::endl;
    int idx_max_element;
    for(unsigned int i=0; i<string.images.size(); i++){
       if (*max == string.E[i]) idx_max_element = i;
    }
    std::cout << idx_max_element << "   " << std::endl;
    af::array tau_s=string.images[idx_max_element].m - string.images[idx_max_element-1].m;
    tau_s=renormalize(tau_s);
    //TODO add detRK4 of stoch_integrator with fdmdt= -gradV + 2 (-gradV, tau_s) tau_s, -gradV is  detfdmdt
    return 0;
}
