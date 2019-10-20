#include "material.hpp"

namespace magnumaf{


Material::Material(float D, float D_axis_x, float D_axis_y, float D_axis_z, float p, float J_atom, float D_atom, float K_atom, float D_atom_axis_x , float D_atom_axis_y, float D_atom_axis_z, float K_atom_axis_x, float K_atom_axis_y, float K_atom_axis_z, bool hexagonal_close_packed)
: D(D), D_axis{D_axis_x, D_axis_y, D_axis_z}, p(p), J_atom(J_atom), D_atom(D_atom), K_atom(K_atom), D_atom_axis{D_atom_axis_x, D_atom_axis_y, D_atom_axis_z}, K_atom_axis{K_atom_axis_x, K_atom_axis_y, K_atom_axis_z}, hexagonal_close_packed(hexagonal_close_packed)
{
}

//void Material::set_atomistic_from_micromagnetic(float dx){
    //this->p = this->ms * pow(dx, 3);
    //this->J_atom =2. * this->A * dx;
    //this->D_atom = this->D * pow(dx, 2);
    //this->K_atom = this->Ku1 * pow(dx, 3);
//}

//void Material::set_atomistic_from_micromagnetic(float dx, float ms, float A, float D, float Ku1){
    //this->p = this->ms * pow(dx, 3);
    //this->J_atom =2. * this->A * dx;
    //this->D_atom = this->D * pow(dx, 2);
    //this->K_atom = this->Ku1 * pow(dx, 3);
//}
}// namespace magnumaf
