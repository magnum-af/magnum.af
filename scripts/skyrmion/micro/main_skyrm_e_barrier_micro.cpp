#include "magnum.af.hpp"

using namespace magnumafcpp;

//TODO this file was used to see wether the skyrmion changes by employing boundary conditions as stated in Paper in implementation. This currently lacks unit testing
using namespace af; typedef std::shared_ptr<LLGTerm> llgt_ptr;
int main(int argc, char** argv)
{
  std::cout<<"argc"<<argc<<std::endl;
   for (int i=0; i<argc; i++)
        cout << "Parameter " << i << " was " << argv[i] << "\n";
  if(argc>1) setDevice(std::stoi(argv[2]));
  info();

  //char* charptr;
  //std::cout<<"argv"<<std::strtod(argv[1], &charptr)<<std::endl;
  //std::cout<<"argv"<<std::strtod(argv[2], &charptr)<<std::endl;

  // Parameter initialization
  const int nx = 111, ny=111 , nz=1;

  const double dx=2.715e-10;//dx max = sqrt(A/Ku1)=1.58114e-9
  const double x=111*dx, y=111*dx, z=4*dx;


  //Simulation Parameters
  //double hmax = 3.5e-10;
  //double hmin = 1.0e-15;

  //double atol = 1e-6;
  //
  //double rtol = atol;

  double n_interp = 60;
  double string_dt=5e-13;
  const int string_steps = 300;
  std::string filepath(argc>0? argv[1]: "../Data/Testing/");
  if(argc>0)filepath.append("/");
  //else std::string filepath("../Data/Testing/");
  std::cout<<"Writing into path "<<filepath.c_str()<<std::endl;


  //Generating Objects
  Mesh mesh(nx, ny, nz, x/nx, y/ny, z/nz);
  Material material = Material();
  material.gamma = 2.211e5;
  state.Ms    = 1.1e6;
  material.A     = 1.6e-11;
  material.alpha = 1;
  state.material.afsync  = false;

  material.D=2*5.76e-3;
  //material.D_axis[2]=-1;

  material.Ku1=6.4e6;
  //material.Ku1_axis[2]=1;

  material.J_atom=4.*material.A*dx;
  material.D_atom= 2.* material.D * pow(dx, 2);
  material.K_atom=4.* material.Ku1/state.Ms;
  material.p=state.Ms*pow(dx, 3);//Compensate nz=1 instead of nz=4


   // Initial magnetic field
   array m = constant(0.0, mesh.n0, mesh.n1, mesh.n2, 3, f64);
   m(span, span, span, 2) = -1;
   for(int ix=0;ix<mesh.n0;ix++){
     for(int iy=0;iy<mesh.n1;iy++){
       const double rx=double(ix)-mesh.n0/2.;
       const double ry=double(iy)-mesh.n1/2.;
       const double r = sqrt(pow(rx, 2)+pow(ry, 2));
       if(r>nx/4.) m(ix, iy, span, 2)=1.;
     }
   }

  State state(mesh, material, m);
  std::cout << "test" << std::endl;
  vti_writer_micro(state.m, mesh , (filepath + "minit").c_str());

  std::cout << "test" << std::endl;
  std::vector<llgt_ptr> llgterm;
  llgterm.push_back( llgt_ptr (new DemagField(mesh, material)));
  std::cout << "test" << std::endl;
  llgterm.push_back( llgt_ptr (new ExchangeField(mesh, material)));
  llgterm.push_back( llgt_ptr (new DmiField(mesh, material)));
  llgterm.push_back( llgt_ptr (new UniaxialAnisotropyField(mesh, material)));

  //llgterm.push_back( llgt_ptr (new AtomisticDipoleDipoleField(mesh)));
  //llgterm.push_back( llgt_ptr (new AtomisticExchangeField(mesh)));
  //llgterm.push_back( llgt_ptr (new AtomisticDmiField(mesh, material)));
  //llgterm.push_back( llgt_ptr (new AtomisticUniaxialAnisotropyField(mesh, material)));

  LLG Llg(state, llgterm);

  timer t = af::timer::start();
  while (state.t < 2.e-10){
    state.m=Llg.step(state);
  }
  double timerelax= af::timer::stop(t);
  vti_writer_micro(state.m, mesh , (filepath + "relax").c_str());

  std::cout<<"timerelax [af-s]: "<< timerelax << " for "<<Llg.counter_accepted+Llg.counter_reject<<" steps, thereof "<< Llg.counter_accepted << " Steps accepted, "<< Llg.counter_reject<< " Steps rejected" << std::endl;




  array last   = constant( 0, mesh.dims, f64);
  last(span, span, span, 2)=1;

  std::vector<State> inputimages;
  inputimages.push_back(state);
  inputimages.push_back(State(mesh, material, last));

  String string(state, inputimages, n_interp, string_dt , llgterm);
  //String* string = new String(state, inputimages, n_interp , llgterm);
  std::cout.precision(12);

  std::ofstream stream_E_barrier;
  stream_E_barrier.precision(12);
  stream_E_barrier.open ((filepath + "E_barrier.dat").c_str());

  std::ofstream stream_steps;
  stream_steps.precision(12);
  stream_steps.open ((filepath + "steps.dat").c_str());

  std::ofstream stream_E_curves;
  stream_E_curves.precision(12);
  stream_E_curves.open ((filepath + "E_curves.dat").c_str());

  double max_lowest=1e100;
  double max_prev_step=1e100;
  int i_max_lowest=-1;
  std::vector<State> images_max_lowest;
  std::vector<double> E_max_lowest;
  for(int i=0; i<string_steps;i++){
    string.step();
    string.calc_E();
    auto max = std::max_element(std::begin(string.E), std::end(string.E));
    if(*max-string.E[0]<max_lowest) {
      max_lowest=*max-string.E[0];
      i_max_lowest=i;
      images_max_lowest=string.images;
      E_max_lowest=string.E;
    }
    else if(i>50){
      std::cout   << "Exiting loop: Energy barrier after 50step relaxation becomes bigger "<<std::endl;
      stream_steps<<"#Exiting loop: Energy barrier after 50step relaxation becomes bigger "<<std::endl;
      break;
    }
    //std::cout<<"Test: fabs= "<<fabs(2*(*max-string.E[0]-max_prev_step)/(*max-string.E[0]+max_prev_step))<<std::endl;
    if(i>25 && fabs(2*(*max-string.E[0]-max_prev_step)/(*max-string.E[0]+max_prev_step))<1e-6){
      std::cout   <<      "Exiting loop: Energy barrier relative difference smaller than 1e-6"<<std::endl;
      stream_steps<<     "#Exiting loop: Energy barrier relative difference smaller than 1e-6"<<std::endl;
      break;
    }
    std::cout   <<i<<"\t"<<*max-string.E[0]<<"\t"<<*max-string.E[-1]<< "\t"<<*max<<"\t"<<fabs(2*(*max-string.E[0]-max_prev_step)/(*max-string.E[0]+max_prev_step))<<std::endl;
    stream_steps<<i<<"\t"<<*max-string.E[0]<<"\t"<<*max-string.E[-1]<< "\t"<<*max<<"\t"<<fabs(2*(*max-string.E[0]-max_prev_step)/(*max-string.E[0]+max_prev_step))<<std::endl;
    for(unsigned j=0;j<string.E.size();++j)
    {
      stream_E_curves<<i<<" "<<j<<" "<<string.E[j]-string.E[0]<<" "<<string.E[j]-string.E[-1]<<" "<<string.E[j]<<std::endl;
    }
    stream_E_curves<<i<<"\n\n"<<std::endl;
    max_prev_step=*max-string.E[0];
  }
    std::cout   <<"#i , lowest overall:   max-[0], max-[-1] max [J]: "<<i_max_lowest<<"\t"<<max_lowest<<"\t"<<max_lowest+E_max_lowest[0]-E_max_lowest[-1]<<"\t"<<max_lowest+E_max_lowest[0]<< std::endl;
    stream_steps<<"#i , lowest overall:   max-[0], max-[-1] max [J]: "<<i_max_lowest<<"\t"<<max_lowest<<"\t"<<max_lowest+E_max_lowest[0]-E_max_lowest[-1]<<"\t"<<max_lowest+E_max_lowest[0]<< std::endl;
    stream_E_barrier<<max_lowest<<"\t"<<std::endl;

  std::ofstream myfileE;
  myfileE.precision(12);
  myfileE.open ((filepath + "E_last_step.dat").c_str());

  std::ofstream stream_max_lowest;
  stream_max_lowest.precision(12);
  stream_max_lowest.open ((filepath + "E_max_lowest.dat").c_str());

  std::cout<<string.E.size()<<"\t"<<string.images.size()<< "\t" <<std::endl;
  for(unsigned i = 0; i < string.images.size(); i++){
    std::cout<<"i="<<i<< "\t" << "E= "<<string.E[i]<<std::endl;
    myfileE<<i<< "\t" << string.E[i]<< "\t" << string.E[i]-string.E[0]<< "\t" << string.E[i]-string.E[-1]<<std::endl;
    std::string name = filepath;
    name.append("skyrm_image");
    name.append(std::to_string(i));
    vti_writer_micro(string.images[i].m, mesh , name.c_str());
    //af_to_vtk(string.images_interp[i], name.c_str());
    stream_max_lowest<<i<< "\t" << E_max_lowest[i]<<"\t" << E_max_lowest[i]-E_max_lowest[0]<<"\t" << E_max_lowest[i]-E_max_lowest[-1]<<std::endl;
    name = filepath;
    name.append("skyrm_image_max_lowest");
    name.append(std::to_string(i));
    vti_writer_micro(images_max_lowest[i].m, mesh , name.c_str());
  }

  for(unsigned i=0;i<Llg.Fieldterms.size();++i){
    std::cout<<"get_cpu_time()"<<std::endl;
    std::cout<<i<<"\t"<<Llg.Fieldterms[i]->get_cpu_time()<<std::endl;
    stream_steps<<"#"<<"get_cpu_time()"<<std::endl;
    stream_steps<<"#"<<i<<"\t"<<Llg.Fieldterms[i]->get_cpu_time()<<std::endl;
  }


  myfileE.close();
  stream_E_barrier.close();
  stream_steps.close();
  stream_E_curves.close();
  stream_max_lowest.close();
  //delete[] string;

  return 0;
}
