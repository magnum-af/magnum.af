#include "state.hpp"
#include "vtk_IO.hpp"
#include "func.hpp"
#include "misc.hpp"
#include <iomanip>

namespace magnumafcpp{


/// Overloaded '+' operator adds an af::array to af::array this->m
State State::operator+(const af::array& a) const{
    State result = *this;
    result.m += a;
    return result;
}

void State::set_Ms_field_if_m_minvalnorm_is_zero(const af::array& m, af::array& Ms_field){
    // Initializes Ms_field if any entry of initial m has zero norm
    if(minval(vecnorm(m)) == 0){
        if(verbose) {printf("%s in state.cpp: initial m has values with zero norm, building Ms_field array\n", Info());}
        af::array nzero = !af::iszero(vecnorm(m));
        n_cells_ = afvalue_u32(af::sum(af::sum(af::sum(nzero, 0), 1), 2));
        if(Ms == 0) printf("Wraning: State::set_Ms_field: State.Ms is used but set to zero. It appears that you are using a legacy constuctor. Please pass Ms in constructor!\n");
        Ms_field = af::constant(this->Ms, nzero.dims(), f64);//TODO this yields probem as Ms is not set in constuctor!
        Ms_field *= nzero;
        Ms_field = af::tile(Ms_field, 1, 1, 1, 3);
    }
}

//void State::check_discretization(){
//    if ( this->material.A != 0 && this->material.Ku1 != 0) { // TODO implement better way of checking
//        double max_allowed_cellsize = sqrt(this->material.A/this->material.Ku1);
//        if (verbose && (this->mesh.dx > max_allowed_cellsize || this->mesh.dy > max_allowed_cellsize || this->mesh.dz > max_allowed_cellsize )){
//            if( ! mute_warning) printf("%s State::check_discretization: cell size is too large (greater than sqrt(A/Ku1)\n", Warning());
//        }
//    }
//}

//void State::check_nonequispaced_discretization(){
//    if ( this->material.A != 0 && this->material.Ku1 != 0) { // TODO implement better way of checking
//        const double max_allowed_cellsize = sqrt(this->material.A/this->material.Ku1);
//        const double max_dz = *std::max_element(nonequimesh.z_spacing.begin(), nonequimesh.z_spacing.end());
//        if (verbose && (this->mesh.dx > max_allowed_cellsize || this->mesh.dy > max_allowed_cellsize || max_dz > max_allowed_cellsize )){
//            if( ! mute_warning) printf("%s State::check_discretization: cell size is too large (greater than sqrt(A/Ku1)\n", Warning());
//        }
//    }
//}

void State::check_m_norm(double tol){//allowed norm is 1 or 0 (for no Ms_field)
    af::array one_when_value_is_zero = af::iszero(vecnorm(m));
    double meannorm = afvalue(af::mean(af::mean(af::mean(af::mean(vecnorm(m)+1.*one_when_value_is_zero, 0), 1), 2), 3));
    if ( (fabs(meannorm - 1.) > tol) && ( this->mute_warning == false )) {
        printf("%s State::check_m_norm: non-zero parts of the magnetization are not normalized to 1! Results won't be physically meaningfull.\n", Warning());
    }
}
//long int State::get_m_addr(){
//    u_out = this->m.copy();
//    return (long int) m_out.get();
//}
//

State::State (Mesh mesh, double Ms, af::array m, bool verbose, bool mute_warning): mesh(mesh), Ms(Ms), m(m), verbose(verbose), mute_warning(mute_warning)
{
    check_m_norm();
    set_Ms_field_if_m_minvalnorm_is_zero( this->m, this->Ms_field);
    //check_discretization();
}


State::State (Mesh mesh, double Ms, long int m, bool verbose, bool mute_warning): mesh(mesh), Ms(Ms), m(*(new af::array( *((void **) m)))), verbose(verbose), mute_warning(mute_warning)
{
    check_m_norm();
    set_Ms_field_if_m_minvalnorm_is_zero( this->m, this->Ms_field);
    //check_discretization();
}


State::State (Mesh mesh, af::array Ms_field, af::array m, bool verbose, bool mute_warning): mesh(mesh), Ms_field(Ms_field), m(m), verbose(verbose), mute_warning(mute_warning)
{
    check_m_norm();
}


State::State (Mesh mesh, long int Ms_field_ptr, long int m, bool verbose, bool mute_warning): mesh(mesh), Ms_field(*(new af::array( *((void **) Ms_field_ptr)))), m(*(new af::array( *((void **) m)))), verbose(verbose), mute_warning(mute_warning)
{
    check_m_norm();
}


State::State (Mesh mesh, Material param, af::array m, bool verbose, bool mute_warning): mesh(mesh), material(param), m(m), verbose(verbose), mute_warning(mute_warning)
{
    check_m_norm();
    set_Ms_field_if_m_minvalnorm_is_zero( this->m, this->Ms_field);
    //check_discretization();
}


State::State (NonequispacedMesh nonequimesh, double Ms, af::array m, bool verbose, bool mute_warning):
              nonequimesh(nonequimesh), Ms(Ms), m(m), verbose(verbose), mute_warning(mute_warning)
{
    check_m_norm();
    set_Ms_field_if_m_minvalnorm_is_zero( this->m, this->Ms_field);
}


State::State (NonequispacedMesh nonequimesh, af::array Ms_field, af::array m, bool verbose, bool mute_warning):
              nonequimesh(nonequimesh), Ms_field(Ms_field), m(m), verbose(verbose), mute_warning(mute_warning)
{
    check_m_norm();
    set_Ms_field_if_m_minvalnorm_is_zero( this->m, this->Ms_field);
}


///< State method taking additional boolean array for specific mean evaluation where this array is true (==1)
State::State (Mesh mesh_in, Material param_in, af::array m_in, af::array evaluate_mean):
              mesh(mesh_in), material(param_in), m(m_in), evaluate_mean_(evaluate_mean)
{
    set_Ms_field_if_m_minvalnorm_is_zero( this->m, this->Ms_field);
    //check_discretization();
    check_m_norm();
    evaluate_mean_is_1_ = afvalue_u32(af::sum(af::sum(af::sum(evaluate_mean_, 0), 1), 2));
    evaluate_mean_ = af::tile(evaluate_mean_, 1, 1, 1, 3);// expanding to 3 vector dimensions, now calculating evaluate_mean_is_1_ would be 3 times too high
    if (verbose) printf("%s state.cpp: evaluate_mean_is_1_= %u\n", Info(), evaluate_mean_is_1_);
}

State::State (Mesh mesh_in, Material param_in, long int aptr): mesh(mesh_in), material(param_in), verbose(false)
{
    void **a = (void **)aptr;
    m = *( new af::array( *a ));
    //m.lock();
    set_Ms_field_if_m_minvalnorm_is_zero( this->m, this->Ms_field);
    //check_discretization();
    check_m_norm();
}

///< For wrapping: State method taking additional boolean array for specific mean evaluation where this array is true (==1)
State::State (Mesh mesh_in, Material param_in, long int aptr, long int evaluate_mean_ptr): mesh(mesh_in), material(param_in), verbose(false)
{
    void **a = (void **)aptr;
    m = *( new af::array( *a ));
    //m.lock();

    void **b = (void **)evaluate_mean_ptr;
    evaluate_mean_ = *( new af::array( *b ));
    //evaluate_mean_.lock();

    set_Ms_field_if_m_minvalnorm_is_zero( this->m, this->Ms_field);
    //check_discretization();
    check_m_norm();
    evaluate_mean_is_1_ = afvalue_u32(af::sum(af::sum(af::sum(evaluate_mean_, 0), 1), 2));
    evaluate_mean_ = af::tile(evaluate_mean_, 1, 1, 1, 3);// expanding to 3 vector dimensions, now calculating evaluate_mean_is_1_ would be 3 times too high
}


void State::Normalize(){
    this->m = renormalize(this->m);
}

void State::set_m(long int aptr){
    void **a = (void **)aptr;
    m = *( new af::array( *a ));
    check_m_norm();
}

long int State::get_m_addr(){
    af::array *a = new af::array(m);
    return (long int) a->get();
}

void State::set_Ms_field(long int aptr){
    void **a = (void **)aptr;
    Ms_field = *( new af::array( *a )); // TODO rename Ms_field -> micro_Ms_field
}

long int State::get_Ms_field(){
    af::array *a = new af::array(Ms_field);
    return (long int) a->get();
}

void State::write_vti(std::string outputname){
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


void State::vtr_writer(std::string outputname){
    ::magnumafcpp::vtr_writer(this->m, this->nonequimesh, outputname, false);
}
void State::vtr_reader(std::string inputname){
    ::magnumafcpp::vtr_reader(this->m, this->nonequimesh, inputname, false);
}

double State::meani(const int i){
    double *norm_host=NULL;
    if (!evaluate_mean_.isempty()){
        //af::print ("temp", temp);
        //std::cout << "tem type = "<< temp.type() << std::endl;
        ///< Calculates the mean values for the specified values given in evaluate_mean_
        norm_host = (af::sum(af::sum(af::sum((m * evaluate_mean_)(af::span, af::span, af::span, i), 0), 1), 2)/evaluate_mean_is_1_).host<double>();
    }
    else if(!Ms_field.isempty() && n_cells_ != 0){
        if (n_cells_ == 0) printf("%s State::meani: n_cells_ is empty and will be divieded by 0!\n", red("Warning:").c_str());
        norm_host = (af::sum(af::sum(af::sum(m(af::span, af::span, af::span, i), 0), 1), 2)/n_cells_).host<double>();
    }
    else{
        norm_host = af::mean(af::mean(af::mean(m(af::span, af::span, af::span, i), 0), 1), 2).host<double>();
    }
    double norm = norm_host[0];
    //std::cout << "norm_host = "<< norm << std::endl;
    af::freeHost(norm_host);
    return norm;
}

///< Writing to filestrean: state.t, <mx>,  <my>,  <mz>
void State::calc_mean_m(std::ostream& myfile ){
    if(Ms_field.isempty()){
        af::array mean_dim3 = af::mean(af::mean(af::mean(this->m, 0), 1), 2);
        myfile << std::setw(12) << this->t << "\t" << afvalue(mean_dim3(af::span, af::span, af::span, 0)) << "\t" << afvalue(mean_dim3(af::span, af::span, af::span, 1))<< "\t" << afvalue(mean_dim3(af::span, af::span, af::span, 2)) << std::endl;
    }
    else{
        af::array sum_dim3 = af::sum(af::sum(af::sum(this->m, 0), 1), 2);
        myfile << std::setw(12) << this->t << "\t" << afvalue(sum_dim3(af::span, af::span, af::span, 0))/n_cells_ << "\t" << afvalue(sum_dim3(af::span, af::span, af::span, 1))/n_cells_<< "\t" << afvalue(sum_dim3(af::span, af::span, af::span, 2))/n_cells_ << std::endl;
    }

}

///< Writing to filestrean: state.t, <mx>,  <my>,  <mz>, hzee
void State::calc_mean_m( std::ostream& myfile, double hzee){
    if(Ms_field.isempty()){
        af::array mean_dim3 = af::mean(af::mean(af::mean(this->m, 0), 1), 2);
        myfile << std::setw(12) << this->t << "\t" << afvalue(mean_dim3(af::span, af::span, af::span, 0)) << "\t" << afvalue(mean_dim3(af::span, af::span, af::span, 1))<< "\t" << afvalue(mean_dim3(af::span, af::span, af::span, 2)) <<  "\t" << hzee << std::endl;
    }
    else{
        af::array sum_dim3 = af::sum(af::sum(af::sum(this->m, 0), 1), 2);
        myfile << std::setw(12) << this->t << "\t" << afvalue(sum_dim3(af::span, af::span, af::span, 0))/n_cells_ << "\t" << afvalue(sum_dim3(af::span, af::span, af::span, 1))/n_cells_<< "\t" << afvalue(sum_dim3(af::span, af::span, af::span, 2))/n_cells_ <<  "\t" << hzee << std::endl;
    }
}

///< Writing to filestrean: state.t, <mx>,  <my>,  <mz>, hzee_x, hzee_y, hzee_z
void State::calc_mean_m( std::ostream& myfile, const af::array& hzee){
    af::array sum_dim3 = sum(sum(sum(this->m, 0), 1), 2);
    if(Ms_field.isempty()){
        af::array mean_dim3 = af::mean(af::mean(af::mean(this->m, 0), 1), 2);
        myfile << std::setw(12) << this->t << "\t" << afvalue(mean_dim3(af::span, af::span, af::span, 0)) << "\t" << afvalue(mean_dim3(af::span, af::span, af::span, 1))<< "\t" << afvalue(mean_dim3(af::span, af::span, af::span, 2)) << "\t" << afvalue(hzee(0, 0, 0, 0)) << "\t" << afvalue(hzee(0, 0, 0, 1)) << "\t" << afvalue(hzee(0, 0, 0, 2)) << std::endl;
    }
    else{
        af::array sum_dim3 = af::sum(af::sum(af::sum(this->m, 0), 1), 2);
        myfile << std::setw(12) << this->t << "\t" << afvalue(sum_dim3(af::span, af::span, af::span, 0))/n_cells_ << "\t" << afvalue(sum_dim3(af::span, af::span, af::span, 1))/n_cells_<< "\t" << afvalue(sum_dim3(af::span, af::span, af::span, 2))/n_cells_ << "\t" << afvalue(hzee(0, 0, 0, 0)) << "\t" << afvalue(hzee(0, 0, 0, 1)) << "\t" << afvalue(hzee(0, 0, 0, 2)) << std::endl;
    }
}

///< Writing to filestrean: state.steps, <mx>,  <my>,  <mz>, hzee
void State::calc_mean_m_steps( std::ostream& myfile, double hzee){
    if(Ms_field.isempty()){
        af::array mean_dim3 = af::mean(af::mean(af::mean(this->m, 0), 1), 2);
        myfile << std::setw(12) << this->steps << "\t" << afvalue(mean_dim3(af::span, af::span, af::span, 0)) << "\t" << afvalue(mean_dim3(af::span, af::span, af::span, 1))<< "\t" << afvalue(mean_dim3(af::span, af::span, af::span, 2)) <<  "\t" << hzee << std::endl;
    }
    else{
        af::array sum_dim3 = af::sum(af::sum(af::sum(this->m, 0), 1), 2);
        myfile << std::setw(12) << this->steps << "\t" << afvalue(sum_dim3(af::span, af::span, af::span, 0))/n_cells_ << "\t" << afvalue(sum_dim3(af::span, af::span, af::span, 1))/n_cells_<< "\t" << afvalue(sum_dim3(af::span, af::span, af::span, 2))/n_cells_ <<  "\t" << hzee << std::endl;
    }
}

///< Writing to filestrean: state.steps, <mx>,  <my>,  <mz>, hzee_x, hzee_y, hzee_z
void State::calc_mean_m_steps( std::ostream& myfile, const af::array& hzee){
    af::array sum_dim3 = sum(sum(sum(this->m, 0), 1), 2);
    if(Ms_field.isempty()){
        af::array mean_dim3 = af::mean(af::mean(af::mean(this->m, 0), 1), 2);
        myfile << std::setw(12) << this->steps << "\t" << afvalue(mean_dim3(af::span, af::span, af::span, 0)) << "\t" << afvalue(mean_dim3(af::span, af::span, af::span, 1))<< "\t" << afvalue(mean_dim3(af::span, af::span, af::span, 2)) << "\t" << afvalue(hzee(0, 0, 0, 0)) << "\t" << afvalue(hzee(0, 0, 0, 1)) << "\t" << afvalue(hzee(0, 0, 0, 2)) << std::endl;
    }
    else{
        af::array sum_dim3 = af::sum(af::sum(af::sum(this->m, 0), 1), 2);
        myfile << std::setw(12) << this->steps << "\t" << afvalue(sum_dim3(af::span, af::span, af::span, 0))/n_cells_ << "\t" << afvalue(sum_dim3(af::span, af::span, af::span, 1))/n_cells_<< "\t" << afvalue(sum_dim3(af::span, af::span, af::span, 2))/n_cells_ << "\t" << afvalue(hzee(0, 0, 0, 0)) << "\t" << afvalue(hzee(0, 0, 0, 1)) << "\t" << afvalue(hzee(0, 0, 0, 2)) << std::endl;
    }
}


// Calculate nonequi distant mesh integral:  integral(M * Hdemag) dx, where M = Ms * m
double State::integral_nonequimesh(const af::array& h_times_m) const{
    af::array z_spacing_afarray = af::array(1, 1, this->nonequimesh.nz, 1, this->nonequimesh.z_spacing.data());
    af::array ms_h_times_m;

    // Global or local Ms switch
    if (this->Ms_field.isempty() == true){
        ms_h_times_m = this->Ms * h_times_m;
;
    }
    else {
        ms_h_times_m = this->Ms_field * h_times_m;
    }

    af::array xy_integral = af::sum( af::sum( af::sum( ms_h_times_m, 0), 1), 3) * this->nonequimesh.dx * this->nonequimesh.dy;
    af::array xyz_integral = af::sum(xy_integral * z_spacing_afarray, 2);
    return afvalue( xyz_integral );
}
}// namespace magnumafcpp
