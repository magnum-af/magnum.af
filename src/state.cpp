#include "state.hpp"
#include "func.hpp"

void State::set_Ms_if_m_minvalnorm_is_zero(const af::array& m, af::array& Ms){
    // Initializes Ms if any entry of initial m has zero norm
    if(minval(vecnorm(m)) == 0){
        if(this->verbose_ > 0) {std::cout << "Info: in state.cpp: initial m has values with zero norm, building Ms array" << std::endl;}
        af::array nzero = !af::iszero(vecnorm(m));
        n_cells_ = afvalue_u32(af::sum(af::sum(af::sum(nzero,0), 1), 2));
        Ms = af::constant(this->param.ms, nzero.dims(), f64);
        Ms *= nzero;
        Ms = af::tile(Ms,1,1,1,3);
    }
}

void State::check_discretization(){
    if ( this->param.A != 0 && this->param.Ku1 != 0) { // TODO implement better way of checking
        double max_allowed_cellsize = sqrt(this->param.A/this->param.Ku1);
        if (this->verbose_ > 0 && (this->mesh.dx > max_allowed_cellsize || this->mesh.dy > max_allowed_cellsize || this->mesh.dz > max_allowed_cellsize )){
            std::cout << "Warning: State::check_discretization: cell size is too large (greater than sqrt(A/Ku1) " << std::endl;
        }
    }
}
//long int State::get_m_addr(){
//    u_out = this->m.copy();
//    return (long int) m_out.get();
//}
//
State::State (Mesh mesh_in, Param param_in, af::array m_in):
              mesh(mesh_in),param(param_in), m(m_in)
{
    set_Ms_if_m_minvalnorm_is_zero( this->m, this->Ms);
    check_discretization();
}

State::State (Mesh mesh_in, Param param_in, long int aptr):
              mesh(mesh_in),param(param_in)
{
    void **a = (void **)aptr;
    m = *( new af::array( *a ));
    m.lock();
    set_Ms_if_m_minvalnorm_is_zero( this->m, this->Ms);
    check_discretization();
}

void State::_vti_writer_micro(std::string outputname){
    vti_writer_micro(m, mesh, outputname); 
}
void State::_vti_writer_atom (std::string outputname){
    vti_writer_atom(m, mesh, outputname); 
}
void State::_vti_reader(std::string inputname){
    vti_reader(m, mesh, inputname);
}


void State::_vtr_writer(std::string outputname){
    vtr_writer(m, mesh, outputname); 
}
void State::_vtr_reader(std::string inputname){
    vtr_reader(m, mesh, inputname); 
}

double State::meani(const int i){
    double *norm_host=NULL;
    if(Ms.isempty()){
        norm_host = af::mean(af::mean(af::mean(m(af::span,af::span,af::span,i),0),1),2).host<double>();
    }
    else{
        norm_host = (af::sum(af::sum(af::sum(m(af::span,af::span,af::span,i),0),1),2)/n_cells_).host<double>();
    }
  double norm = norm_host[0];
  af::freeHost(norm_host);
  return norm;
}

///< Writing to filestrean: state.t, <mx>,  <my>,  <mz>
void State::calc_mean_m(std::ostream& myfile ){
    if(Ms.isempty()){
        af::array mean_dim3 = af::mean(af::mean(af::mean(this->m,0),1),2);
        myfile << std::setw(12) << this->t << "\t" << afvalue(mean_dim3(af::span,af::span,af::span,0)) << "\t" << afvalue(mean_dim3(af::span,af::span,af::span,1))<< "\t" << afvalue(mean_dim3(af::span,af::span,af::span,2)) << std::endl;
    }
    else{
        af::array sum_dim3 = af::sum(af::sum(af::sum(this->m,0),1),2);
        myfile << std::setw(12) << this->t << "\t" << afvalue(sum_dim3(af::span,af::span,af::span,0))/n_cells_ << "\t" << afvalue(sum_dim3(af::span,af::span,af::span,1))/n_cells_<< "\t" << afvalue(sum_dim3(af::span,af::span,af::span,2))/n_cells_ << std::endl;
    }

}

///< Writing to filestrean: state.t, <mx>,  <my>,  <mz>, hzee
void State::calc_mean_m( std::ostream& myfile, double hzee){
    if(Ms.isempty()){
        af::array mean_dim3 = af::mean(af::mean(af::mean(this->m,0),1),2);
        myfile << std::setw(12) << this->t << "\t" << afvalue(mean_dim3(af::span,af::span,af::span,0)) << "\t" << afvalue(mean_dim3(af::span,af::span,af::span,1))<< "\t" << afvalue(mean_dim3(af::span,af::span,af::span,2)) <<  "\t" << hzee << std::endl;
    }
    else{
        af::array sum_dim3 = af::sum(af::sum(af::sum(this->m,0),1),2);
        myfile << std::setw(12) << this->t << "\t" << afvalue(sum_dim3(af::span,af::span,af::span,0))/n_cells_ << "\t" << afvalue(sum_dim3(af::span,af::span,af::span,1))/n_cells_<< "\t" << afvalue(sum_dim3(af::span,af::span,af::span,2))/n_cells_ <<  "\t" << hzee << std::endl;
    }
}

///< Writing to filestrean: state.t, <mx>,  <my>,  <mz>, hzee_x, hzee_y, hzee_z
void State::calc_mean_m( std::ostream& myfile, const af::array& hzee){
    af::array sum_dim3 = sum(sum(sum(this->m,0),1),2);
    if(Ms.isempty()){
        af::array mean_dim3 = af::mean(af::mean(af::mean(this->m,0),1),2);
        myfile << std::setw(12) << this->t << "\t" << afvalue(mean_dim3(af::span,af::span,af::span,0)) << "\t" << afvalue(mean_dim3(af::span,af::span,af::span,1))<< "\t" << afvalue(mean_dim3(af::span,af::span,af::span,2)) << "\t" << afvalue(hzee(0,0,0,0)) << "\t" << afvalue(hzee(0,0,0,1)) << "\t" << afvalue(hzee(0,0,0,2)) << std::endl;
    }
    else{
        af::array sum_dim3 = af::sum(af::sum(af::sum(this->m,0),1),2);
        myfile << std::setw(12) << this->t << "\t" << afvalue(sum_dim3(af::span,af::span,af::span,0))/n_cells_ << "\t" << afvalue(sum_dim3(af::span,af::span,af::span,1))/n_cells_<< "\t" << afvalue(sum_dim3(af::span,af::span,af::span,2))/n_cells_ << "\t" << afvalue(hzee(0,0,0,0)) << "\t" << afvalue(hzee(0,0,0,1)) << "\t" << afvalue(hzee(0,0,0,2)) << std::endl;
    }
}

///< Writing to filestrean: state.steps, <mx>,  <my>,  <mz>, hzee
void State::calc_mean_m_steps( std::ostream& myfile, double hzee){
    if(Ms.isempty()){
        af::array mean_dim3 = af::mean(af::mean(af::mean(this->m,0),1),2);
        myfile << std::setw(12) << this->steps << "\t" << afvalue(mean_dim3(af::span,af::span,af::span,0)) << "\t" << afvalue(mean_dim3(af::span,af::span,af::span,1))<< "\t" << afvalue(mean_dim3(af::span,af::span,af::span,2)) <<  "\t" << hzee << std::endl;
    }
    else{
        af::array sum_dim3 = af::sum(af::sum(af::sum(this->m,0),1),2);
        myfile << std::setw(12) << this->steps << "\t" << afvalue(sum_dim3(af::span,af::span,af::span,0))/n_cells_ << "\t" << afvalue(sum_dim3(af::span,af::span,af::span,1))/n_cells_<< "\t" << afvalue(sum_dim3(af::span,af::span,af::span,2))/n_cells_ <<  "\t" << hzee << std::endl;
    }
}

///< Writing to filestrean: state.steps, <mx>,  <my>,  <mz>, hzee_x, hzee_y, hzee_z
void State::calc_mean_m_steps( std::ostream& myfile, const af::array& hzee){
    af::array sum_dim3 = sum(sum(sum(this->m,0),1),2);
    if(Ms.isempty()){
        af::array mean_dim3 = af::mean(af::mean(af::mean(this->m,0),1),2);
        myfile << std::setw(12) << this->steps << "\t" << afvalue(mean_dim3(af::span,af::span,af::span,0)) << "\t" << afvalue(mean_dim3(af::span,af::span,af::span,1))<< "\t" << afvalue(mean_dim3(af::span,af::span,af::span,2)) << "\t" << afvalue(hzee(0,0,0,0)) << "\t" << afvalue(hzee(0,0,0,1)) << "\t" << afvalue(hzee(0,0,0,2)) << std::endl;
    }
    else{
        af::array sum_dim3 = af::sum(af::sum(af::sum(this->m,0),1),2);
        myfile << std::setw(12) << this->steps << "\t" << afvalue(sum_dim3(af::span,af::span,af::span,0))/n_cells_ << "\t" << afvalue(sum_dim3(af::span,af::span,af::span,1))/n_cells_<< "\t" << afvalue(sum_dim3(af::span,af::span,af::span,2))/n_cells_ << "\t" << afvalue(hzee(0,0,0,0)) << "\t" << afvalue(hzee(0,0,0,1)) << "\t" << afvalue(hzee(0,0,0,2)) << std::endl;
    }
}
