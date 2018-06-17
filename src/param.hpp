#ifndef PARAM_H
#define PARAM_H
#include<math.h>
struct Param{
    double mu0=4e-7 * M_PI;		// [H/m] 	// vaccum permeability	//= pi/2500000 [N/A^2] //TODO const leads to error in string.cpp
    double kb = 1.38064852e-23;		// [J/K]  	// Boltzmann Constant
    double gamma{221276.1488637255};	// [m/A/s]	// gyromagnetic ratio	// gamma = gamma_e * mu0; electron gyromagnetic ratio multiplied by magnetic constant
    double alpha{0};			// []		// Unitless Damping constant
    double T{0};			// [K]  	// Temperature 				//in stochastic.cpp: Temperature in Stochastic Integrator

    // Micromagneitc 
    double ms{0};			// [J/T/m^3] 
    double A{0};			// [J/m]
    double D{0.};			// [D/m^2] 	// DM interaction constant  
    double Ku1{0};			// [J/m^3] 	// Uniaxial Anisotropy 
    double D_axis[3]={0,0,-1};		//		// DMI axis
    double Ku1_axis[3]={0,0,1};		//		// Anisotropy axis

    // Atomistic
    double p{0};			// [J/T]  	// Atomistic magnetic moment
    double J_atom{0.};			// [J]   	// Atomistic exchange
    double D_atom{0.};			// [J]   	// Atomistic DMI
    double K_atom{0.};			// [J]   	// Atomistic anisotropy
    double D_atom_axis[3]={0.,0.,1.};		 	// Atomistic DMI axis
    double K_atom_axis[3]={0.,0.,1.};		 	// Atomistic anisotropy axis
    bool  hexagonal_close_packed{false};                // Selects hexagonal close packed atomistic structure

    // non-physical-parameters
    int mode{6};  			// 		// Inegration method	//TODO skipp inparam 0 -> rk4, 1 -> rk4_3o8, 2 -> rkf, 3 -> explicitEuler
    bool afsync{false};			//		// activate af::sync for timings //TODO skipp in param

};
#endif

//Note Js=mu0*Ms
