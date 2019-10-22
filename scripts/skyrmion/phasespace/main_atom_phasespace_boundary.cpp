#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace magnumafcpp;


using namespace af; typedef std::shared_ptr<LLGTerm> llgt_ptr;
int main(int argc, char** argv)
{
    std::string filepath(argc >= 1? argv[1]: "data");
    if( argc >= 1 ){ filepath.append("/");}
    if( argc >= 2 ){ setDevice( std::stoi( argv[2]));}
    std::cout<<"Writing into path "<<filepath.c_str()<<std::endl;
    std::string path_mrelax(argc>3? argv[3]: "");
    info();

    // Parameter initialization
    const int nx = 112, ny=112 , nz=1;
    const double dx=2.715e-10;


    double n_interp = 60;
    double string_dt=5e-14;
    const int string_max_steps = 10000;
    double rel_diff = 1e-12;
    double abs_diff = 1e-27;


    //Generating Objects
    Mesh mesh(nx, ny, nz, dx, dx, dx);
    Material material = Material();
    state.Ms    = 1.1e6;
    material.A     = 1.6e-11;
    material.alpha = 1;
    state.material.afsync  = false;
    material.D = (argc >= 3 ? std::stod(argv[3]) : 0.01152);
    material.Ku1=(argc >= 4 ? std::stod(argv[4]) : 6400000);

    material.set_atomistic_from_micromagnetic(mesh.dx);

    std::cout<<"D="<<material.D<<std::endl;
    std::cout<<"Ku1="<<material.Ku1<<std::endl;
    std::cout<<"D_atom="<<material.D_atom<<std::endl;
    std::cout<<"Ku1_atom="<<material.K_atom<<std::endl;

    State state(mesh, material, mesh.skyrmconf());
    vti_writer_micro(state.m, mesh , (filepath + "minit").c_str());

    LLGIntegrator Llg;
    Llg.llgterms.push_back( LlgTerm (new AtomisticDipoleDipoleField(mesh)));
    Llg.llgterms.push_back( LlgTerm (new AtomisticExchangeField(mesh)));
    Llg.llgterms.push_back( LlgTerm (new AtomisticDmiField(mesh, material)));
    Llg.llgterms.push_back( LlgTerm (new AtomisticUniaxialAnisotropyField(mesh, material)));

    Llg.relax(state);
    vti_writer_micro(state.m, mesh , filepath + "relax");
    state.t=0;


    array last   = constant( 0, mesh.dims, f64);
    last(span, span, span, 2)=1;

    std::vector<State> inputimages;
    for(int i=0; i < n_interp; i++){
        array mm = array(state.m);
        mm=shift(mm, 2*i);
        if(2*i < nx){
            mm(seq(0, 2*i), span, span, span)=0;
            mm(seq(0, 2*i), span, span, 2)=1.;
        }
        else {
            mm(span, span, span, span)=0;
            mm(span, span, span, 2)=1.;
        }
        inputimages.push_back(State(mesh, material, mm));
    }

    String string(state, inputimages, n_interp, string_dt , Llg.llgterms);
    string.run(filepath, rel_diff, abs_diff, string_max_steps);
    return 0;
}
