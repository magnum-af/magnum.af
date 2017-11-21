#ifndef PARAM_H
#define PARAM_H
#include<math.h>
class Param{
  public:
    ~Param(){}; 
    double mu0=4e-7 * M_PI;
    double gamma{0},ms{0},A{0},alpha{0};
    double p{0}; //ATOMISTIC_DEMAG dipole p
    bool afsync{false};
    int mode{6};  // 0 -> rk4, 1 -> rk4_3o8, 2 -> rkf, 3 -> explicitEuler
    double D{0.};//DMI //Note Js=mu0*Ms
    //Uniaxial anisotropy
    double Ku1{0};//Anisotropy 
    double D_axis[3]={0,0,-1};//DMI axis
    //double eu_x{0.},eu_y{0.},eu_z{1.};
    double Ku1_axis[3]={0,0,1};//Anisotropy axis

    //Atomistic
    double J_atom{0.};//Exchange
    double D_atom{0.};//Atomistic
    double D_atom_axis[3]={0.,0.,-1.};//Atomistic
    double K_atom{0.};//Atomistic
    double K_atom_axis[3]={0.,0.,1.};//Atomistic

    //Param (double gamma, double ms, double A, double alpha, bool afsync, int mode);
    //Param (double gamma, double ms, double A, double alpha, bool afsync, int mode, double D, double*  D_axis, double Ku1, double* Ku1_axis);
};
#endif
