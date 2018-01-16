#include <list>
#include "arrayfire.h"
#include "pth_mag.hpp"

using namespace af; 
typedef std::shared_ptr<LLGTerm> llgt_ptr; 

void calcm(State state, std::ostream& myfile){
    myfile << std::setw(12) << state.t << "\t" <<meani(state.m,0)<< "\t" <<meani(state.m,1)<< "\t" <<meani(state.m,2) << "\t" ;
    myfile <<fabs(meani(state.m,0))<< "\t" <<fabs(meani(state.m,1))<< "\t" <<fabs(meani(state.m,2))<< std::endl;
}

class Detector{
    public:
        Detector(unsigned int length_in, double threshold_in): length( length_in ), threshold( threshold_in ){}
        bool gt {true};// avg greater than (gt) threshold:  if true, avg > threshold, else avg < threshold
        unsigned int length;// = 3000;
        double threshold; // = 0.75;
        bool event {false};// The event is the annihilation/decay of the skyrmion
        void add_data(double value){
            if(data.size()<length){data.push_back(value);}
            else{
                data.pop_front(); 
                data.push_back(value); 
                avg=0;
                for( double values : data){
                    avg+=values;
                }
                avg/=(double)length;
                if(gt == true  && avg > threshold) {event = true;}
                if(gt == false && avg < threshold) {event = true;}
            }
        }
        double get_avg(){return avg;}


        std::list<double> data; //The data to take the average on
        double avg{NaN};
    private:
    
};

int main(int argc, char** argv)
{
//    Detector detector = Detector ();
//    int i = -1;
//    detector.threshold = 3000;
//    while (detector.event == false ){
//        i++;
//    //for (int i=0; i < 4000; i++){
//        detector.add_data(i);
//        std::cout << i <<", "<< detector.data.size()<<", "<< detector.data.front()<<" , "<< detector.data.back()<<" , "<< detector.get_avg()<<std::endl;
//    }

    std::cout<<"argc = "<<argc<<std::endl;
     for (int i=0; i<argc; i++)
          cout << "Parameter " << i << " was " << argv[i] << "\n";
    
    std::string filepath(argc>1? argv[1]: "../Data/skyrmion_stoch/Testing");
    if(argc>0)filepath.append("/");
    std::cout<<"Writing into path "<<filepath.c_str()<<std::endl;

    setDevice(argc>2? std::stoi(argv[2]):0);
    //if(argc>1) setDevice(std::stoi(argv[2]));
    info();

    // Parameter initialization
    const int nx = 30, ny=30 ,nz=1;
    const double dx=1e-9;
  
    const double dt = 1e-13;//Integration step
    const double maxtime = 1e-5;
    std::cout << "Simulating "<<maxtime<<" [s] with " << maxtime/dt << " steps, estimated computation time is " << maxtime/dt*3e-3 << " [s] " << std::endl;
  
    //Generating Objects
    Mesh mesh(nx,ny,nz,dx,dx,dx);
    Param param = Param();
    param.gamma = 2.211e5;
    param.ms    = 1.1e6;
    param.A     = 1.6e-11;
    param.alpha = 1;
    param.D=2*5.76e-3;
    param.Ku1=6.4e6;
    param.T = 450;
  
    param.J_atom=2.*param.A*dx;
    param.D_atom= param.D * pow(dx,2)/2.;
    param.K_atom=param.Ku1*pow(dx,3);
    param.p=param.ms*pow(dx,3);//Compensate nz=1 instead of nz=4

    array m; 
    Mesh testmesh(nx,ny,nz,dx,dx,dx);
    vti_reader(m, testmesh, "../E_barrier/relax.vti");
    //vti_reader(m, testmesh, filepath+"E_barrier/relax.vti");
    //vti_writer_atom(m, mesh ,(filepath + "test_readin").c_str());
    
    //Reading energy barrier from calculation performed with main_n30.cpp
    double e_barrier;
    //ifstream stream(filepath+"E_barrier/E_barrier.dat");
    ifstream stream("../E_barrier/E_barrier.dat");
    if (stream.is_open()){
        stream >> e_barrier;
    }
    stream.close();
    std::cout.precision(12);
    std::cout <<"E_barrier from prev calculation = "<< e_barrier <<std::endl;

    double A = pow(10,9);
    double k = pow(10,8);
    double T = e_barrier/(param.kb*(log(A)-log(k)));
    std::cout << "Calculated T for decay in 1e-8 sec = "<< T << std::endl;
    //param.T = T;
    
    
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


    Detector detector = Detector(3000, 0.75);
    unsigned long int i=0;
    while (detector.event == false){
        if(state.t > maxtime){
             std::cout << "WARNING: no event detected within maxtime seconds, aborting" << std::endl;
             break;
        }
    //unsigned long int steps = 1e6;
    //for (unsigned long int i = 0; i < steps; i++){
        Stoch.step(state,dt); 
        calcm(state,stream_m);
        detector.add_data(meani(state.m,2));
        if(i < 100)  vti_writer_atom(state.m, state.mesh, filepath+"/skyrm/dense_skyrm"+std::to_string(i));//state.t*pow(10,9)
        if(i % 100 == 0) vti_writer_atom(state.m, state.mesh, filepath+"/skyrm/skyrm"+std::to_string(i));//state.t*pow(10,9)
        if(i % 1000 == 0) std::cout << i << ", mz= "<< meani(state.m,2) << std::endl;
        if(i % 1000 == 0) std::cout << i << ", avg_mz= "<< detector.get_avg() << std::endl;
        i++;
    }
    
    double detect_time = state.t;
    Detector detector2 = Detector(100, 0.75);
    detector2.gt=false;
    int reverse=0;
    for (std::list<double>::reverse_iterator rit=detector.data.rbegin(); rit!=detector.data.rend(); ++rit){
        detector2.add_data(*rit);
        reverse++;
        if(detector2.event == true) break;
    }
    detect_time-=reverse*dt;
    std::cout<< "Preliminiary annihilation time at " << detect_time << "[s]" << std::endl;
    

    stream_m.close();
    std::cout<<"state.t = "<< state.t << std::endl;
  
    return 0;
}
