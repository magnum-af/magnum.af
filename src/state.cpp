#include "state.hpp"
#include "func.hpp"
#include "misc.hpp"

void State::set_Ms_if_m_minvalnorm_is_zero(const af::array& m, af::array& Ms){
    // Initializes Ms if any entry of initial m has zero norm
    if(minval(vecnorm(m)) == 0){
        if(verbose_) {std::cout << "Info: in state.cpp: initial m has values with zero norm, building Ms array" << std::endl;}
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
        if (verbose_ && (this->mesh.dx > max_allowed_cellsize || this->mesh.dy > max_allowed_cellsize || this->mesh.dz > max_allowed_cellsize )){
            std::cout << red("Warning: State::check_discretization: cell size is too large (greater than sqrt(A/Ku1)") << std::endl;
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

///< State method taking additional boolean array for specific mean evaluation where this array is true (==1)
State::State (Mesh mesh_in, Param param_in, af::array m_in, af::array evaluate_mean):
              mesh(mesh_in),param(param_in), m(m_in), evaluate_mean_(evaluate_mean)
{
    set_Ms_if_m_minvalnorm_is_zero( this->m, this->Ms);
    check_discretization();
    evaluate_mean_is_1_ = afvalue_u32(af::sum(af::sum(af::sum(evaluate_mean_,0), 1), 2));
    evaluate_mean_ = af::tile(evaluate_mean_, 1, 1, 1, 3);// expanding to 3 vector dimensions, now calculating evaluate_mean_is_1_ would be 3 times too high
    if (verbose_) std::cout << "Info: state.cpp: evaluate_mean_is_1_= " << evaluate_mean_is_1_ << std::endl;
}

State::State (Mesh mesh_in, Param param_in, long int aptr): mesh(mesh_in), param(param_in), verbose_(false)
{
    void **a = (void **)aptr;
    m = *( new af::array( *a ));
    m.lock();
    set_Ms_if_m_minvalnorm_is_zero( this->m, this->Ms);
    check_discretization();
}

///< For wrapping: State method taking additional boolean array for specific mean evaluation where this array is true (==1)
State::State (Mesh mesh_in, Param param_in, long int aptr, long int evaluate_mean_ptr): mesh(mesh_in), param(param_in), verbose_(false)
{
    void **a = (void **)aptr;
    m = *( new af::array( *a ));
    m.lock();

    void **b = (void **)evaluate_mean_ptr;
    evaluate_mean_ = *( new af::array( *b ));
    evaluate_mean_.lock();

    set_Ms_if_m_minvalnorm_is_zero( this->m, this->Ms);
    check_discretization();
    evaluate_mean_is_1_ = afvalue_u32(af::sum(af::sum(af::sum(evaluate_mean_,0), 1), 2));
    evaluate_mean_ = af::tile(evaluate_mean_, 1, 1, 1, 3);// expanding to 3 vector dimensions, now calculating evaluate_mean_is_1_ would be 3 times too high
}

void State::_vti_writer_micro(std::string outputname){
    vti_writer_micro(m, mesh, outputname); 
}
void State::_vti_writer_micro_boolean(std::string outputname){
    vti_writer_micro(evaluate_mean_(af::span, af::span, af::span, 0).as(f64), mesh, outputname); //NOTE: as evaluate_mean_ is tiles to 3 dims, taking only first
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
    if (!evaluate_mean_.isempty()){
        //af::print ("temp", temp);
        //std::cout << "tem type = "<< temp.type() << std::endl;
        ///< Calculates the mean values for the specified values given in evaluate_mean_
        norm_host = (af::sum(af::sum(af::sum((m * evaluate_mean_)(af::span,af::span,af::span,i),0),1),2)/evaluate_mean_is_1_).host<double>();
    }
    else if(!Ms.isempty()){
        norm_host = (af::sum(af::sum(af::sum(m(af::span,af::span,af::span,i),0),1),2)/n_cells_).host<double>();
    }
    else{
        norm_host = af::mean(af::mean(af::mean(m(af::span,af::span,af::span,i),0),1),2).host<double>();
    }
    double norm = norm_host[0];
    //std::cout << "norm_host = "<< norm << std::endl;
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
