#include "param.hpp"
void Param::set_atomistic_from_micromagnetic(const double dx){
    this->p = this->ms * pow(dx,3);
    this->J_atom =2. * this->A * dx;
    this->D_atom = this->D * pow(dx,2);
    this->K_atom = this->Ku1 * pow(dx,3);
}
