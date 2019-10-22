#include <list>
#include <sstream>
#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace magnumafcpp;


using namespace af;
typedef std::shared_ptr<LLGTerm> llgt_ptr;

void calcm(State state, std::ostream& myfile, double get_avg){
    myfile << std::setw(12) << state.t << "\t" <<meani(state.m, 2) << "\t" << get_avg << std::endl;
    //myfile << std::setw(12) << state.t << "\t" <<meani(state.m, 0)<< "\t" <<meani(state.m, 1)<< "\t" <<meani(state.m, 2) << "\t" << get_avg << std::endl;
    //myfile <<fabs(meani(state.m, 0))<< "\t" <<fabs(meani(state.m, 1))<< "\t" <<fabs(meani(state.m, 2))<< std::endl;
}

void set_boundary_mz(array& m);
class Detector{
    public:
        Detector(unsigned int length_in, double threshold_in): length( length_in ), threshold( threshold_in ){}
        bool gt {true};// avg greater than (gt) threshold:  if true, avg > threshold, else avg < threshold
        unsigned int length;// = 3000;
        double threshold; // = 0.75;
        bool event {false};// The event is the annihilation/decay of the skyrmion
        double get_avg(){return avg;}

        std::list<double> data; //The data to take the average on
        std::list<double> data_avg; //The calculated avarages of the data
        void add_data(double value);
        double avg{NaN};
    private:

};

void Detector::add_data(double value){
    if(data.size()<length){
        data.push_back(value);
    }
    else{
        data.pop_front();
        data.push_back(value);
        avg=0;
        for( double values : data){
            avg+=values;
        }
        avg/=(double)length;
        if(data_avg.size()<length){
            data_avg.push_back(avg);
        }
        else {
            data_avg.pop_front();
            data_avg.push_back(avg);
        }
        if(gt == true  && avg > threshold) {event = true;}
        if(gt == false && avg < threshold) {event = true;}
    }
}

int main(int argc, char** argv)
{
    std::cout.precision(32);
    std::cout<<time<<std::endl;
    // .magnum.af /path/to/write path/to/read dt T ID GPU
    // .magnum.af $PWD  ~/git/magnum.af/Data/skyrmion_stoch/E_barrier 1e-9 300 0 0
    std::cout<<"argc = "<<argc<<std::endl;
    for (int i=0; i<argc; i++){
        cout << "Parameter " << i << " was " << argv[i] << "\n";
    }

    std::string filepath(argc>1? argv[1]: "../Data/skyrmion_stoch/Testing");
    filepath.append("/");
    std::cout<<"Writing into path "<<filepath<<std::endl;

    std::string indatapath(argc>2? argv[2]: "/home/pth/git/magnum.af/Data/skyrmion_stoch/E_barrier");
    indatapath.append("/");

    //Timestep
    double dt;
    if(argc>3){
        std::istringstream os(argv[3]); // Reading double in scientific notation
        os >> dt;
    }
    else{
        dt = 1e-13;
    }
    std::cout << "dt = " << dt << std::endl;

    //Temperature
    double T_in(argc>4? std::stod(argv[4]):300);

    int ID (argc>5? std::stoi(argv[5]) : -1);
    //const double dt = 1e-13;//Integration step

    //GPU device
    setDevice(argc>6? std::stoi(argv[6]):0);
    //if(argc>1) setDevice(std::stoi(argv[2]));
    info();

    // Parameter initialization
    const int nx = 30, ny=30 , nz=1;
    const double dx=1e-9;

    const double threshold = 0.85;
    const int samples = 1;
    const double maxtime = 1e-5;

    std::list<af::array> images; //The data to take the average on
    unsigned int images_length=25;

    //std::cout << "Simulating "<<maxtime<<" [s] with " << maxtime/dt << " steps, estimated computation time is " << maxtime/dt*3e-3 << " [s] " << std::endl;

    //Generating Objects
    Mesh mesh(nx, ny, nz, dx, dx, dx);
    Material material = Material();
    state.Ms    = 1.1e6;
    material.A     = 1.6e-11;
    material.alpha = 1;
    material.D=2*5.76e-3;
    material.Ku1=6.4e6;
    material.T = T_in;
    //material.T = argc>3? std::stod(argv[3]):300;
    std::cout<< "material.T = " << material.T << std::endl;

    material.J_atom=2.*material.A*dx;
    material.D_atom= material.D * pow(dx, 2);
    material.K_atom=material.Ku1*pow(dx, 3);
    material.p=state.Ms*pow(dx, 3);//Compensate nz=1 instead of nz=4

    array m;
    Mesh testmesh(nx, ny, nz, dx, dx, dx);
    vti_reader(m, testmesh, indatapath + "relax.vti");
    //Reading energy barrier from calculation performed with main_n30.cpp
    double e_barrier;
    //ifstream stream(filepath+"E_barrier/E_barrier.dat");
    ifstream stream(indatapath + "E_barrier.dat");
    if (stream.is_open()){
        stream >> e_barrier;
    }
    stream.close();
    std::cout.precision(12);
    std::cout <<"E_barrier from prev calculation = "<< e_barrier <<std::endl;

    double A = pow(10, 9);
    double k = pow(10, 8);
    double T = e_barrier/(constants::kb*(log(A)-log(k)));
    std::cout << "Calculated T for decay in 1e-8 sec = "<< T << std::endl;
    //material.T = T;


    //Assemble Stochastic Integrator and State object
    State state(mesh, material, m);
    std::vector<llgt_ptr> llgterm;
    llgterm.push_back( llgt_ptr (new AtomisticExchangeField(mesh)));
    llgterm.push_back( llgt_ptr (new AtomisticDmiField(mesh, material)));
    llgterm.push_back( llgt_ptr (new AtomisticUniaxialAnisotropyField(mesh, material)));
    Stochastic_LLG Stoch(state, llgterm, dt , "Heun");

    //ofs_antime<<"#detect_time << \t << state.t << \t << i <<\t << reverse " << std::endl;
  //  std::vector<double> antimes;// appends annihilation time for each iteration to calculate mean
    // Iterations to obtain mean annihilaiton time
    state = State(mesh, material, m);
    std::cout<<"TEST (should be 0): state.t = "<< state.t << std::endl;

    Detector detector = Detector( (int)(3e-10/dt), threshold);
    unsigned long int i=0;
    while (detector.event == false){
        if(state.t > maxtime){
             std::cout << "WARNING: no event detected within maxtime seconds, aborting" << std::endl;
             break;
        }
        Stoch.step(state, dt);
        //set_boundary_mz(state.m);
        detector.add_data(meani(state.m, 2));
        //if( ID < id_outp ) calcm(state, stream_m, detector.get_avg());
        //if(i < 100)  vti_writer_atom(state.m, state.mesh, filepath+"/vti/dense_skyrm"+std::to_string(ID)+"_"+std::to_string(i));//state.t*pow(10, 9)
        //if(i < 100)  vti_writer_atom(state.m, state.mesh, filepath+"/skyrm/dense_skyrm"+std::to_string(j)+"_"+std::to_string(i));//state.t*pow(10, 9)
        //if(i % (int)(1e-10/dt) == 0) vti_writer_atom(state.m, state.mesh, filepath+"/skyrm/skyrm"+std::to_string(j)+"_"+std::to_string(i));//state.t*pow(10, 9)
        //if( ID < id_outp && i % (int)(1e-10/dt) == 0) vti_writer_atom(state.m, state.mesh, filepath+"/vti/skyrm"+std::to_string(ID)+"_"+std::to_string(i));//state.t*pow(10, 9)
        //if( ID < id_outp && i % (int)(1e-10/dt) == 0) std::cout << state.t << " , i= " << i << ", mz    = "<< meani(state.m, 2) << std::endl;
        //if( ID < id_outp && i % (int)(1e-10/dt) == 0) std::cout << state.t << " , i= " << i << ", avg_mz= "<< detector.get_avg() << std::endl;
        if( images.size() < images_length &&  i % (int)(1e-11/dt) == 0) images.push_back(state.m);
        if( images.size() > images_length &&  i % (int)(1e-11/dt) == 0) {images.pop_front(); images.push_back(state.m);}
        i++;
    }

    //Writing mz and mz_avag into file mID.dat
    std::ofstream stream_m;
    stream_m.open(filepath + "m"+std::to_string(ID)+".dat");
    stream_m.precision(12);
    for (auto const& i : detector.data) {
        stream_m << i <<std::endl;
    }
    stream_m.close();

    if(detector.data_avg.size()>0){
        stream_m.open(filepath + "m_avg"+std::to_string(ID)+".dat");
        for (auto const& i : detector.data_avg) {
            stream_m << i <<std::endl;
        }
        stream_m.close();
    }

    vti_writer_atom(state.m, state.mesh, filepath+"skyrm"+std::to_string(ID));//state.t*pow(10, 9)

    //last images:
    int count1=0;
    for (auto const& i : images) {
        vti_writer_atom(i, state.mesh, filepath+"/vti/skyrm"+std::to_string(ID)+"_"+std::to_string(count1));//state.t*pow(10, 9)
        count1++;
    }
    //vti_writer_atom(state.m, state.mesh, filepath+"skyrm"+std::to_string(j));//state.t*pow(10, 9)

    double detect_time = state.t;
    Detector detector2 = Detector((int)(5e-12/dt), threshold);
    detector2.gt=false;
    int reverse=0;
    for (std::list<double>::reverse_iterator rit=detector.data.rbegin(); rit!=detector.data.rend(); ++rit){
        detector2.add_data(*rit);
        reverse++;
        if(detector2.event == true) break;
    }
    detect_time-=reverse*dt;
    std::cout<< "Preliminiary annihilation time at " << detect_time << "[s]" << std::endl;

    std::ofstream ofs_antime(filepath + "anihilationtime.dat", std::ios_base::app);
    ofs_antime.precision(12);
    ofs_antime<< detect_time << "\t" << state.t << "\t" << i <<"\t" << reverse << "\t" << ID << std::endl;
    ofs_antime.close();


    std::cout<<"state.t = "<< state.t << std::endl;
    return 0;
}
void set_boundary_mz(array& m){
    m( 0, span, span, 0)=0.;
    m(-1, span, span, 0)=0.;
    m(span,  0, span, 0)=0.;
    m(span, -1, span, 0)=0.;

    m( 0, span, span, 1)=0.;
    m(-1, span, span, 1)=0.;
    m(span,  0, span, 1)=0.;
    m(span, -1, span, 1)=0.;

    m( 0, span, span, 2)=1.;
    m(-1, span, span, 2)=1.;
    m(span,  0, span, 2)=1.;
    m(span, -1, span, 2)=1.;
}

    //vti_reader(m, testmesh, "../../E_barrier/relax.vti");
    //set_boundary_mz(m);
    //vti_reader(m, testmesh, filepath+"E_barrier/relax.vti");
    //vti_writer_atom(m, mesh , (filepath + "test_readin").c_str());

  //  double mean_antime = 0;
  //  for (double n : antimes){
  //      mean_antime+=n;
  //  }
  //  mean_antime/= (double) antimes.size();

  //  double unbiased_sample_variance = 0; // s^2= 1/(n-1) sum(y_i - y_mean)^2 from i = 1 to n
  //  for (double n : antimes){
  //      unbiased_sample_variance += pow( n - mean_antime , 2);
  //  }
  //  unbiased_sample_variance/=( (double)antimes.size() -1 );
  //  double unbiased_sample_sigma = sqrt(unbiased_sample_variance);

  //  std::cout<<"antimes.size() = "<< antimes.size() << std::endl;
  //  std::cout<<"mean_antime = "<< mean_antime << std::endl;
  //  std::cout<<"unbiased_sample_variance = "<< unbiased_sample_variance << std::endl;
  //  std::cout<<"unbiased_sample_sigma = "<< unbiased_sample_sigma << std::endl;
  //  ofs_antime.open(filepath + "mean_antime"+std::to_string(ID)+".dat"); ofs_antime.precision(12);
  //  ofs_antime<<"#mean_antime << unbiased_sample_variance << unbiased_sample_sigma<< antimes.size() "<< std::endl;
  //  ofs_antime<< mean_antime <<  "\t" << unbiased_sample_variance <<"\t" << unbiased_sample_sigma<< "\t" << antimes.size() <<  std::endl;
  //  ofs_antime.close();

    //  antimes.push_back(detect_time);

//    Detector detector = Detector ();
//    int i = -1;
//    detector.threshold = 3000;
//    while (detector.event == false ){
//        i++;
//    //for (int i=0; i < 4000; i++){
//        detector.add_data(i);
//        std::cout << i <<", "<< detector.data.size()<<", "<< detector.data.front()<<" , "<< detector.data.back()<<" , "<< detector.get_avg()<<std::endl;
//    }

