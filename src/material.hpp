#pragma once
namespace magnumaf{

///Struct holding all simulation parameters.
struct Material{
    // Micromagneitc
    float D{0.};			//!< [D/m^2] 	// DM interaction constant
    float D_axis[3]={0, 0, -1};		//!<		// DMI axis

    // Atomistic
    void set_atomistic_from_micromagnetic(float dx);
    void set_atomistic_from_micromagnetic(float dx, float ms, float A, float D, float Ku1);
    float p{0};			//!< [J/T]  	// Atomistic magnetic moment
    float J_atom{0.};			//!< [J]   	// Atomistic exchange
    float D_atom{0.};			//!< [J]   	// Atomistic DMI
    float K_atom{0.};			//!< [J]   	// Atomistic anisotropy
    float D_atom_axis[3]={0., 0., 1.};		 	//!< Atomistic DMI axis
    float K_atom_axis[3]={0., 0., 1.};		 	//!< Atomistic anisotropy axis
    bool  hexagonal_close_packed{false};                //!< Selects hexagonal close packed atomistic structure

    // non-physical-parameters
    bool afsync{false};			//!< activate af::sync for timings //TODO skipp in material

    // Default constructor
    Material(){};
    // For wrapping only
    Material(float D, float D_axis_x, float D_axis_y, float D_axis_z, float p, float J_atom, float D_atom, float K_atom, float D_atom_axis_x , float D_atom_axis_y, float D_atom_axis_z, float K_atom_axis_x, float K_atom_axis_y, float K_atom_axis_z, bool hexagonal_close_packed);
};

//Note Js=mu0*Ms
}// namespace magnumaf
