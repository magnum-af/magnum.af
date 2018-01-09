#include "string.hpp"
#include "func.hpp"

using namespace af;

String::String(State statein, std::vector<State> inputimages, int n_interp_in, double dt_in, std::vector<std::shared_ptr<LLGTerm> > Fieldterms_in):
  state(statein), Llg(state, Fieldterms_in), n_interp(n_interp_in), dt(dt_in){

  //If set true, this only uses the energy dissipation term (i.e. Mx(MxH)) in the LLG
    Llg.fdmdt_dissipation_term_only=true;

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
            images[i].m=Llg.llgstep(images[i]);
        }
        double h=imagtime+dt-images[i].t;
        double dummy_err;
        images[i].m += Llg.RKF45(images[i].m,h,dummy_err);
        //af::eval(images[i].m);//If memory error occurs, uncomment this 
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
//      images[i].m=Llg.llgstep(images[i]);
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

