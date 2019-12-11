#pragma once
namespace magnumafcpp{

///Struct holding all simulation parameters.
struct Material{
    // Atomistic
    void set_atomistic_from_micromagnetic(double dx);
    void set_atomistic_from_micromagnetic(double dx, double ms, double A, double D, double Ku1);
    double p{0};			//!< [J/T]  	// Atomistic magnetic moment
    double J_atom{0.};			//!< [J]   	// Atomistic exchange
    double D_atom{0.};			//!< [J]   	// Atomistic DMI
    double K_atom{0.};			//!< [J]   	// Atomistic anisotropy
    double D_atom_axis[3]={0., 0., 1.};		 	//!< Atomistic DMI axis
    double K_atom_axis[3]={0., 0., 1.};		 	//!< Atomistic anisotropy axis
    bool  hexagonal_close_packed{false};                //!< Selects hexagonal close packed atomistic structure

    // non-physical-parameters
    bool afsync{false};			//!< activate af::sync for timings //TODO skipp in material

    // Default constructor
    Material(){};
    // For wrapping only
    Material(double p, double J_atom, double D_atom, double K_atom, double D_atom_axis_x , double D_atom_axis_y, double D_atom_axis_z, double K_atom_axis_x, double K_atom_axis_y, double K_atom_axis_z, bool hexagonal_close_packed);
};

//Note Js=mu0*Ms
}// namespace magnumafcpp
