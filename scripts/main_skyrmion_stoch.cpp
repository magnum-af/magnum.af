#include "arrayfire.h"
#include "pth_mag.hpp"

using namespace af; 
typedef std::shared_ptr<LLGTerm> llgt_ptr; 

void calcm(State state, std::ostream& myfile){
    myfile << std::setw(12) << state.t << "\t" <<meani(state.m,0)<< "\t" <<meani(state.m,1)<< "\t" <<meani(state.m,2)<< "\t" << std::endl;
}

int main(int argc, char** argv)
{

    std::cout<<"argc = "<<argc<<std::endl;
     for (int i=0; i<argc; i++)
          cout << "Parameter " << i << " was " << argv[i] << "\n";
    
    std::string filepath(argc>1? argv[1]: "../Data/skyrmion_stoch");
    if(argc>0)filepath.append("/");
    std::cout<<"Writing into path "<<filepath.c_str()<<std::endl;

    setDevice(argc>2? std::stoi(argv[2]):0);
    //if(argc>1) setDevice(std::stoi(argv[2]));
    info();

    // Parameter initialization
    const int nx = 30, ny=30 ,nz=1;
    const double dx=1e-9;
  
    double n_interp = 60;
    double string_dt=1e-13;
    const int string_steps = 10000;
  
    const double dt = 2e-14;//Integration step
  
    //Generating Objects
    Mesh mesh(nx,ny,nz,dx,dx,dx);
    Param param = Param();
    param.gamma = 2.211e5;
    param.ms    = 1.1e6;
    param.A     = 1.6e-11;
    param.alpha = 1;
    param.D=2*5.76e-3;
    param.Ku1=6.4e6;
    param.T = 400;
  
    param.J_atom=2.*param.A*dx;
    param.D_atom= param.D * pow(dx,2)/2.;
    param.K_atom=param.Ku1*pow(dx,3);
    param.p=param.ms*pow(dx,3);//Compensate nz=1 instead of nz=4

    array m; 
    Mesh testmesh(nx,ny,nz,dx,dx,dx);
    vti_reader(m, testmesh, filepath+"E_barrier/relax.vti");
    //vti_writer_atom(m, mesh ,(filepath + "test_readin").c_str());
    
    //Reading energy barrier from calculation performed with main_n30.cpp
    double e_barrier;
    ifstream stream(filepath+"E_barrier/E_barrier.dat");
    if (stream.is_open()){
        stream >> e_barrier;
    }
    stream.close();
    std::cout.precision(12);
    std::cout <<"E_barrier from prev calculation = \n"<< e_barrier <<std::endl;

    double A = pow(10,9);
    double k = pow(10,7);
    double T = e_barrier/(param.kb*(log(A)-log(k)));
    std::cout << "propsed T for decay in 1e-7 sec is \n"<< T << std::endl;
    
    
    std::ofstream stream_m((filepath + "m.dat").c_str());
    stream_m.precision(12);
    //stream_m.open ((filepath + "m.dat").c_str());
 
    //Assemble Stochastic Integrator and State object
    State state(mesh,param, m);
    std::vector<llgt_ptr> llgterm;
    llgterm.push_back( llgt_ptr (new ATOMISTIC_EXCHANGE(mesh)));
    llgterm.push_back( llgt_ptr (new ATOMISTIC_DMI(mesh,param)));
    llgterm.push_back( llgt_ptr (new ATOMISTIC_ANISOTROPY(mesh,param)));
    Stochastic_LLG Stoch(state, llgterm, dt , "Heun");


    unsigned long int steps = 1e5;
    int image=0;
    for (unsigned long int i = 0; i < steps; i++){
        Stoch.step(state,dt); 
        calcm(state,stream_m);
        if(i < 100)  vti_writer_atom(state.m, state.mesh, filepath+"dense_skyrm"+std::to_string(i));//state.t*pow(10,9)
        if(i % 100 == 0) vti_writer_atom(state.m, state.mesh, filepath+"skyrm"+std::to_string(i));//state.t*pow(10,9)
    }
    stream_m.close();
    std::cout<<"state.t = "<< state.t << std::endl;
  
    return 0;
}
