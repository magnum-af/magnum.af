#include "state.hpp"
#include "func.hpp"

void State::set_Ms_if_m_minvalnorm_is_zero(const af::array& m, af::array& Ms){
    // Initializes Ms if any entry of initial m has zero norm
    if(minval(vecnorm(m)) == 0){
        std::cout << "Info: in state.cpp: initial m has values with zero norm, building Ms array" << std::endl;
        af::array nzero = !af::iszero(vecnorm(m));
        Ms = af::constant(this->param.ms, nzero.dims(), f64);
        Ms *= nzero;
        Ms = af::tile(Ms,1,1,1,3);
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
}

State::State (Mesh mesh_in, Param param_in, long int aptr):
              mesh(mesh_in),param(param_in)
{
  void **a = (void **)aptr;
  m = *( new af::array( *a ));
  m.lock();
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
  norm_host = af::mean(af::mean(af::mean(m(af::span,af::span,af::span,i),0),1),2).host<double>();
  double norm = norm_host[0];
  af::freeHost(norm_host);
  return norm;
}

void State::calc_mean_m(std::ostream& myfile ){
    af::array mean_dim3 = af::mean(af::mean(af::mean(this->m,0),1),2);
    myfile << std::setw(12) << this->t << "\t" << afvalue(mean_dim3(af::span,af::span,af::span,0)) << "\t" << afvalue(mean_dim3(af::span,af::span,af::span,1))<< "\t" << afvalue(mean_dim3(af::span,af::span,af::span,2)) << std::endl;
}

void State::calc_mean_m( std::ostream& myfile, const long int n_cells){
    af::array sum_dim3 = af::sum(af::sum(af::sum(this->m,0),1),2);
    myfile << std::setw(12) << this->t << "\t" << afvalue(sum_dim3(af::span,af::span,af::span,0))/n_cells << "\t" << afvalue(sum_dim3(af::span,af::span,af::span,1))/n_cells<< "\t" << afvalue(sum_dim3(af::span,af::span,af::span,2))/n_cells << std::endl;
}

void State::calc_mean_m( std::ostream& myfile, const long int n_cells, double hzee){
    af::array sum_dim3 = af::sum(af::sum(af::sum(this->m,0),1),2);
    myfile << std::setw(12) << this->t << "\t" << afvalue(sum_dim3(af::span,af::span,af::span,0))/n_cells << "\t" << afvalue(sum_dim3(af::span,af::span,af::span,1))/n_cells<< "\t" << afvalue(sum_dim3(af::span,af::span,af::span,2))/n_cells <<  "\t" << hzee << std::endl;
}

void State::calc_mean_m( std::ostream& myfile, const long int n_cells, const af::array& hzee){
    af::array sum_dim3 = sum(sum(sum(this->m,0),1),2);
    myfile << std::setw(12) << this->t << "\t" << afvalue(sum_dim3(af::span,af::span,af::span,0))/n_cells << "\t" << afvalue(sum_dim3(af::span,af::span,af::span,1))/n_cells<< "\t" << afvalue(sum_dim3(af::span,af::span,af::span,2))/n_cells << "\t" << afvalue(hzee(0,0,0,0)) << "\t" << afvalue(hzee(0,0,0,1)) << "\t" << afvalue(hzee(0,0,0,2)) << std::endl;
}
