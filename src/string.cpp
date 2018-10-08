#include "string.hpp"
#include "func.hpp"

using namespace af;

String::String(State statein, std::vector<State> inputimages, int n_interp_in, double dt_in, std::vector<std::shared_ptr<LLGTerm> > Fieldterms_in):
  state(statein), Llg("RKF45", Controller(), true), n_interp(n_interp_in), dt(dt_in){

  //If set true, this only uses the energy dissipation term (i.e. Mx(MxH)) in the LLG
  // in init  Llg.fdmdt_dissipation_term_only=true;
  Llg.llgterms = Fieldterms_in;

    calc_x(inputimages);

    for(int i=0;i<n_interp;i++){
        images.push_back(State(state.mesh, state.param, array(state.mesh.dims,f64)));
    }
    for(int i=0;i<n_interp;i++){
        int j=0;
        while(x[j]<x_interp[i] && j<n_interp) j++;
        if(j>0) j--;
        if(j<n_interp-1){
            images[i].m=inputimages[j].m+(x_interp[i]-x[j])*(inputimages[j+1].m-inputimages[j].m)/(x[j+1]-x[j]);
        }
        else{
            std::cout<<"Warning: x_interp[j="<<j<<"] exceedes x range, y_interp[j] value set to y[j]"<<std::endl;
            images[i]=inputimages[j];
        }
    }
    vec_renormalize();
}

void String::run(const std::string filepath, const double string_abort_rel_diff, const double string_abort_abs_diff, const int string_steps){

    this->write_vti(filepath+"init_string");
    std::cout.precision(12);
  
    std::ofstream stream_E_barrier;
    stream_E_barrier.precision(12);
  
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
        af::timer t = af::timer::start();
        //af::printMemInfo();
        this->step();
        this->calc_E();
        // Update maximal values
        auto max = std::max_element(std::begin(this->E), std::end(this->E));
        if(*max-this->E[0]<max_lowest) {
            max_lowest=*max-this->E[0];
            i_max_lowest=i;
            images_max_lowest=this->images;
            E_max_lowest=this->E;
        }    
        // Check for convergence
        const double abs_diff = fabs(*max-this->E[0]-max_prev_step);
        const double rel_diff = fabs((2.* abs_diff)/(*max-this->E[0]+max_prev_step));
        if(i > 25 && rel_diff < string_abort_rel_diff){
            std::cout <<      "String run: relative difference of two consecutive E_barriers is: rel_diff= " << rel_diff << " This is smaller than " << string_abort_rel_diff <<std::endl;
            stream_steps <<  "#String run: relative difference of two consecutive E_barriers is: rel_diff= " << rel_diff << " This is smaller than " << string_abort_rel_diff <<std::endl;
            break;
        }
        if(i > 25 && abs_diff < string_abort_abs_diff){
            std::cout   <<  "String run: absolute difference of two consecutive E_barriers is: abs_diff= " << abs_diff << " is smaller than " << string_abort_abs_diff <<std::endl;
            stream_steps<< "#String run: absolute difference of two consecutive E_barriers is: abs_diff= " << abs_diff << " is smaller than " << string_abort_abs_diff <<std::endl;
            break;
        }
        stream_E_barrier.open ((filepath + "E_barrier.dat").c_str());
        stream_E_barrier<<max_lowest<<"\t"<<this->images[0].mesh.n0 <<"\t"<<this->images[0].mesh.dx<<"\t"<<this->images[0].param.D<<"\t"<<this->images[0].param.Ku1<<"\t"<<this->images[0].param.K_atom<<"\t"<<this->images[0].param.D_atom<<std::endl;
        stream_E_barrier.close();
        for(unsigned j=0;j<this->E.size();++j)
        {
            stream_E_curves<<i<<" "<<j<<" "<<this->E[j]-this->E[0]<<" "<<this->E[j]-this->E[-1]<<" "<<this->E[j]<<std::endl;
        }
        stream_E_curves<<i<<"\n\n"<<std::endl;
        max_prev_step=*max-this->E[0];
        if(i%50==1){ 
            std::cout<<"Writing current skyrm images for iteration"<<i<<std::endl;
            for(unsigned j = 0; j < this->images.size(); j++){
                std::string name = filepath;
                name.append("current_skyrm_image");
                name.append(std::to_string(j));
                vti_writer_micro(this->images[j].m, this->images[0].mesh ,name.c_str());
            }
        }
        std::cout   << i << "\t" << *max-this->E[0 ]<< "\t rate= "<< 1./(af::timer::stop(t)) << " [1/s]" << std::endl;
        stream_steps<< i << "\t" << *max-this->E[0 ]<< "\t"<< std::endl;
    }
    std::cout   <<"#i ,lowest overall:   max-[0], max-[-1] max [J]: "<<i_max_lowest<<"\t"<<max_lowest<<"\t"<<max_lowest+E_max_lowest[0]-E_max_lowest[-1]<<"\t"<<max_lowest+E_max_lowest[0]<< std::endl;
    stream_steps<<"#i ,lowest overall:   max-[0], max-[-1] max [J]: "<<i_max_lowest<<"\t"<<max_lowest<<"\t"<<max_lowest+E_max_lowest[0]-E_max_lowest[-1]<<"\t"<<max_lowest+E_max_lowest[0]<< std::endl;
  
    std::ofstream myfileE;
    myfileE.precision(12);
    myfileE.open ((filepath + "E_last_step.dat").c_str());
  
    std::ofstream stream_max_lowest;
    stream_max_lowest.precision(12);
    stream_max_lowest.open ((filepath + "E_max_lowest.dat").c_str());
  
    std::cout<<this->E.size()<<"\t"<<this->images.size()<< "\t" <<std::endl;
    for(unsigned i = 0; i < this->images.size(); i++){
      std::cout<<"i="<<i<< "\t" << "E= "<<this->E[i]<<std::endl;
      myfileE<<i<< "\t" << this->E[i]<< "\t" << this->E[i]-this->E[0]<< "\t" << this->E[i]-this->E[-1]<<std::endl;
      std::string name = filepath;
      name.append("skyrm_image");
      name.append(std::to_string(i));
      vti_writer_micro(this->images[i].m, this->images[0].mesh ,name.c_str());
      stream_max_lowest<<i<< "\t" << E_max_lowest[i]<<"\t" << E_max_lowest[i]-E_max_lowest[0]<<"\t" << E_max_lowest[i]-E_max_lowest[-1]<<std::endl;
      name = filepath;
      name.append("skyrm_image_max_lowest");
      name.append(std::to_string(i));
      vti_writer_micro(images_max_lowest[i].m, this->images[0].mesh ,name.c_str());
    }
  
    //for(unsigned i=0;i<Llg.llgterms.size();++i){
    //  std::cout<<"get_cpu_time()"<<std::endl;
    //  std::cout<<i<<"\t"<<Llg.cpu_time()<<std::endl;
    //  stream_steps<<"#"<<"get_cpu_time()"<<std::endl;
    //  stream_steps<<"#"<<i<<"\t"<<Llg.cpu_time()<<std::endl;
    //}
  
    myfileE.close();
    stream_steps.close();
    stream_E_curves.close();
    stream_max_lowest.close();
};



void String::calc_E(){
    if(E.empty()==false) E.clear();
    for(int i=0; i<n_interp; i++) E.push_back(Llg.E(images[i]));
}

void String::calc_x(){
    x.clear();
    x.push_back(0.);
    for(unsigned int i=1; i<images.size(); i++){
        x.push_back(x[i-1] + FrobeniusNorm(images[i].m-images[i-1].m));
    }
    x_interp.clear();
    for(int i=0;i<n_interp;i++){
        x_interp.push_back((double)i/(double)(n_interp-1)*x.back());
    }
}

void String::calc_x(std::vector<State> inputimages){
    x.clear();
    x.push_back(0.);
    for(unsigned int i=1; i<inputimages.size(); i++){
        x.push_back(x[i-1] + FrobeniusNorm(inputimages[i].m-inputimages[i-1].m));
    }
    x_interp.clear();
    for(int i=0;i<n_interp;i++){
        x_interp.push_back((double)i/(double)(n_interp-1)*x.back());
    }
}

void String::lin_interpolate(){
    std::vector<State>images_temp=images;
    for(int i=0;i<n_interp;i++){
        int j=0;
        while(x[j]<x_interp[i] && j<n_interp) j++;
        if(j>0) j--;
        if(j<n_interp-1){
            images[i].m=images_temp[j].m+(x_interp[i]-x[j])*(images_temp[j+1].m-images_temp[j].m)/(x[j+1]-x[j]);
        }
        else{
            std::cout<<"Warning: x_interp[j="<<j<<"] exceedes x range, y_interp[j] value set to y[j]"<<std::endl;
            images[i]=images_temp[j];
        }
        //af::eval(images[i].m);//If memory error occurs, uncomment this, check performance, maybe only evaluate for i%5=0 or so
    }
    vec_renormalize();
}

void String::integrate(){
    for(unsigned int i=0;i<images.size();i++){
        double imagtime=images[i].t;
        while (images[i].t < imagtime + dt){
            Llg.step(images[i]);
        }
        // Now skipping step backwards
        // double h=imagtime+dt-images[i].t;
        // double dummy_err;
        // images[i].m += Llg.RKF45(images[i],h,dummy_err);
        
        // NOTE:
        // af::eval(images[i].m);//If memory error occurs, uncomment this 
    }
}

void String::step(){
    integrate();
    calc_x();
    lin_interpolate();
}

void String::vec_renormalize(){
    for(unsigned int i=0; i<images.size();i++){
        images[i].m= renormalize(images[i].m);
        //af::eval avoids JIT crash here!
        af::eval(images[i].m);
    }
}

void String::write_vti(std::string file){
    for(unsigned j = 0; j < images.size(); j++){
        vti_writer_atom(images[j].m, state.mesh , file + std::to_string(j));
    }
}


//#include "string.hpp"
//#include "func.hpp"
//
//using namespace af;
//
//String::String(State statein, std::vector<State> inputimages, int n_interp_in, std::vector<std::shared_ptr<LLGTerm> > Fieldterms_in):
// state(statein), Llg(state,1e-6,1e-6,3.5e-10,1.0e-15, Fieldterms_in), n_interp(n_interp_in), images(inputimages){
//
//  
//  Llg.fdmdt_dissipation_term_only=true;//To only use the first term in llg equ
//
//  calc_x();
//  for(int i=0;i<n_interp;i++){
//    images_interp.push_back(State(state.mesh, state.param, array(state.mesh.dims,f64)));
//  }
//  lin_interpolate();
//}
//
//void String::calc_E(){
//  if(E.empty()==false) E.clear();
//  for(int i=0; i<n_interp; i++) E.push_back(Llg.E(images_interp[i]));
//  //for(int i=0; i<n_interp; i++) std::cout<<"E= "<<E[i]<<"\t["<<i<<"]"<<std::endl;
//}
//
//void String::calc_x(){
//  x.clear();
//  x.push_back(0.);
//  for(unsigned int i=1; i<images.size(); i++){
//      x.push_back(x[i-1] + FrobeniusNorm(images[i].m-images[i-1].m));
//  }
//
//  x_interp.clear();
//  for(int i=0;i<n_interp;i++){
//      x_interp.push_back((double)i/(double)(n_interp-1)*x.back());
//  }
////  std::cout<<"test1"<<x.size()<<std::endl;
////  for(std::vector<double>::size_type i = 0; i != x.size(); i++){std::cout<< x[i]<<std::endl;}
////  for(std::vector<double>::size_type i = 0; i != x_interp.size(); i++) std::cout<< x_interp[i]<<std::endl;
//}
//
//void String::lin_interpolate(){
//  for(int i=0;i<n_interp;i++){
//    int j=0;
//    while(x[j]<x_interp[i] && j<n_interp) j++;
//    if(j>0) j--;
//    if(j<n_interp-1){
//      images_interp[i].m=images[j].m+(x_interp[i]-x[j])*(images[j+1].m-images[j].m)/(x[j+1]-x[j]);
//    }
//    else{
//      std::cout<<"Warning: x_interp[j="<<j<<"] exceedes x range, y_interp[j] value set to y[j]"<<std::endl;
//      images_interp[i]=images[j];
//    }
//    vec_renormalize();
//  }
//  //for(int n=0; n<images.size();n++) print("images",images[n]);
//  //for(int i=0; i<n_interp; i++) std::cout<< x[i]<<std::endl;
//  //for(int i=0; i<n_interp; i++) std::cout<< x_interp[i]<<std::endl;
//  //for(int n=0; n<n_interp;n++) print("images_interp",images_interp[n]);
//
//  //for(int i=0; i<n_interp; i++){
//  //  std::cout<< "i="<<i<<" x[i]"<< x[i]<<" x_interp[i] "<<x_interp[i]<<std::endl;
//  //  print("images",images[i]);
//  //  print("images_interp",images_interp[i]);
//  //}
//}
//
//void String::integrate(){
//
//  images.clear();
//  for(unsigned int i=0;i<images_interp.size();i++){
//    images.push_back(images_interp[i]);
//    images[i].t=0;//TODO handle time
//    while (images[i].t < dt){
//      images[i].m=Llg.step(images[i]);
//    }
//    double h=dt-images[i].t;
//    double dummy_err;
//    images[i].m += Llg.RKF45(images[i].m,h,dummy_err);
//    //TODO time+=images.t;
//  }
//}
//
//void String::step(){
//  integrate();
//  calc_x();
//  lin_interpolate();
//}
//
//void String::vec_renormalize(){
//  for(unsigned int i=0; i<images.size();i++) images[i].m= renormalize(images[i].m);
//  for(unsigned int i=0; i<images_interp.size();i++) images_interp[i].m= renormalize(images_interp[i].m);
//}

