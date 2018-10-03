#include "atomistic_demag.hpp"
using namespace af;

//Energy calculation
//Edemag=-mu0/2 integral(M . Hdemag) dx
double ATOMISTIC_DEMAG::E(const State& state){
  return -state.param.mu0/2 * state.param.p * afvalue(sum(sum(sum(sum(h(state)*state.m,0),1),2),3));
  //return -state.param.p/2 * afvalue(sum(sum(sum(sum(h(state)*state.m,0),1),2),3));
}

array N_atomistic(int n0_exp, int n1_exp, int n2_exp, double dx, double dy, double dz);

ATOMISTIC_DEMAG::ATOMISTIC_DEMAG (Mesh mesh){
  Nfft=N_atomistic(mesh.n0_exp,mesh.n1_exp,mesh.n2_exp,mesh.dx,mesh.dy,mesh.dz);
  //h   =array (mesh.n0_exp    ,mesh.n1_exp,mesh.n2_exp,3,f64);
}


array ATOMISTIC_DEMAG::h(const State& state){
    timer_demagsolve = timer::start();
  // FFT with zero-padding of the m field
  array mfft;
  if (state.mesh.n2_exp == 1){
    mfft=fftR2C<2>(state.m,dim4(state.mesh.n0_exp,state.mesh.n1_exp));
  }
  else {
    mfft=fftR2C<3>(state.m,dim4(state.mesh.n0_exp,state.mesh.n1_exp,state.mesh.n2_exp));
  }

  array hfft=array (state.mesh.n0_exp/2+1,state.mesh.n1_exp,state.mesh.n2_exp,3,c64);
  // Pointwise product
  hfft(span,span,span,0)= Nfft(span,span,span,0) * mfft(span,span,span,0) 
                                 + Nfft(span,span,span,1) * mfft(span,span,span,1)
                                 + Nfft(span,span,span,2) * mfft(span,span,span,2);
  hfft(span,span,span,1)= Nfft(span,span,span,1) * mfft(span,span,span,0) 
                                 + Nfft(span,span,span,3) * mfft(span,span,span,1)
                                 + Nfft(span,span,span,4) * mfft(span,span,span,2);
  hfft(span,span,span,2)= Nfft(span,span,span,2) * mfft(span,span,span,0)
                                 + Nfft(span,span,span,4) * mfft(span,span,span,1)
                                 + Nfft(span,span,span,5) * mfft(span,span,span,2);

  // IFFT reversing padding
  array h_field;
  if (state.mesh.n2_exp == 1){
    h_field=fftC2R<2>(hfft);
    //af::print("h_dip",state.param.p * h_field(seq(0,state.mesh.n0_exp/2-1),seq(0,state.mesh.n1_exp/2-1)));//TODO hack
    if(state.param.afsync) sync();
    cpu_time += timer::stop(timer_demagsolve);
    return state.param.p * h_field(seq(0,state.mesh.n0_exp/2-1),seq(0,state.mesh.n1_exp/2-1));//TODO consider p density, then we have to multip at m before fft
  }
  else {
    h_field=fftC2R<3>(hfft);
    if(state.param.afsync) sync();
    cpu_time += timer::stop(timer_demagsolve);
    return state.param.p * h_field(seq(0,state.mesh.n0_exp/2-1),seq(0,state.mesh.n1_exp/2-1),seq(0,state.mesh.n2_exp/2-1),span);
  }
}

array N_atomistic(int n0_exp, int n1_exp, int n2_exp, double dx, double dy, double dz){
  double* N = NULL;
  N = new double[n0_exp*n1_exp*n2_exp*6];
  //Experimental
  for (int i0 = 0; i0 < n0_exp; i0++){
    const int j0 = (i0 + n0_exp/2) % n0_exp - n0_exp/2;
    for (int i1 = 0; i1 < n1_exp; i1++){
      const int j1 = (i1 + n1_exp/2) % n1_exp - n1_exp/2;
      for (int i2 = 0; i2 < n2_exp; i2++ ){
        const int j2 = (i2 + n2_exp/2) % n2_exp - n2_exp/2;
        const int idx = 6*(i2+n2_exp*(i1+n1_exp*i0));
        const double rx=j0*dx;
        const double ry=j1*dy;
        const double rz=j2*dz;
        const double r = sqrt(pow(rx,2)+pow(ry,2)+pow(rz,2));
        if(r==0.){ //TODO repsace with if (j0 == 0 && j1 == 0 && j2 == 0)
          //std::cout<<"In ATOMISTIC_DEMAG::N_atomistic: r==0"<<std::endl;
          //std::cout<<"In ATOMISTIC_DEMAG::setting n to 1/3."<<std::endl;
          //Accounting for self-interaction (would be inf, when approximated with sphere -1/3 in diag
          //TODO check
          N[idx+0] = 0.;
          N[idx+1] = 0.;
          N[idx+2] = 0.;
          N[idx+3] = 0.;
          N[idx+4] = 0.;
          N[idx+5] = 0.;
          
          //N[idx+0] = -1./3.;
          //N[idx+1] = 0.;
          //N[idx+2] = 0.;
          //N[idx+3] = -1./3.;
          //N[idx+4] = 0.;
          //N[idx+5] = -1./3.;
          
        }
        else{
          //N[idx+0] = 1./(4.*M_PI)*(3.*rx*rx/pow(r,5) - 1./pow(r,3));
          //N[idx+1] = 1./(4.*M_PI)*(3.*rx*ry/pow(r,5)              );
          //N[idx+2] = 1./(4.*M_PI)*(3.*rx*rz/pow(r,5)              );
          //N[idx+3] = 1./(4.*M_PI)*(3.*ry*ry/pow(r,5) - 1./pow(r,3));
          //N[idx+4] = 1./(4.*M_PI)*(3.*ry*rz/pow(r,5)              );
          //N[idx+5] = 1./(4.*M_PI)*(3.*rz*rz/pow(r,5) - 1./pow(r,3));
          N[idx+0] = 3.*rx*rx/pow(r,5) - 1./pow(r,3);
          N[idx+1] = 3.*rx*ry/pow(r,5)              ;
          N[idx+2] = 3.*rx*rz/pow(r,5)              ;
          N[idx+3] = 3.*ry*ry/pow(r,5) - 1./pow(r,3);
          N[idx+4] = 3.*ry*rz/pow(r,5)              ;
          N[idx+5] = 3.*rz*rz/pow(r,5) - 1./pow(r,3);
        }
      }
    }
  }
  array Naf(6,n2_exp,n1_exp,n0_exp,N);
  Naf=reorder(Naf,3,2,1,0);
  Naf *= 1./(4.*M_PI);
  //print("ATOMISTIC_DEMAG::N_atomistic: Naf",Naf);
  //print("Demag:Naf", Naf(0,0,0,span));
  delete [] N;
  N = NULL;
  if (n2_exp == 1){
    Naf = fftR2C<2>(Naf);
  }
  else {
    Naf = fftR2C<3>(Naf);
  }
  return Naf;
}
