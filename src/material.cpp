#include "material.hpp"

Material::Material(double alpha, double T, double ms, double D, double Ku1, double D_axis_x, double D_axis_y, double D_axis_z, double Ku1_axis_x, double Ku1_axis_y, double Ku1_axis_z, double p, double J_atom, double D_atom, double K_atom, double D_atom_axis_x , double D_atom_axis_y, double D_atom_axis_z, double K_atom_axis_x, double K_atom_axis_y, double K_atom_axis_z, bool hexagonal_close_packed, int mode, bool afsync )
: alpha(alpha), T(T), ms(ms), D(D), Ku1(Ku1), D_axis{D_axis_x, D_axis_y, D_axis_z}, Ku1_axis{Ku1_axis_x, Ku1_axis_y, Ku1_axis_z}, p(p), J_atom(J_atom), D_atom(D_atom), K_atom(K_atom), D_atom_axis{D_atom_axis_x, D_atom_axis_y, D_atom_axis_z}, K_atom_axis{K_atom_axis_x, K_atom_axis_y, K_atom_axis_z}, hexagonal_close_packed(hexagonal_close_packed), mode(mode), afsync(afsync)
{
}

void Material::set_atomistic_from_micromagnetic(double dx){
    this->p = this->ms * pow(dx,3);
    //this->J_atom =2. * this->A * dx;
    this->D_atom = this->D * pow(dx,2);
    this->K_atom = this->Ku1 * pow(dx,3);
}

void Material::set_atomistic_from_micromagnetic(double dx, double ms, double A, double D, double Ku1){
    this->p = this->ms * pow(dx,3);
    //this->J_atom =2. * this->A * dx;
    this->D_atom = this->D * pow(dx,2);
    this->K_atom = this->Ku1 * pow(dx,3);
}
