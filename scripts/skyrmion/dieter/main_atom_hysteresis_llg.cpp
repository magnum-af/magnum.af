#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace magnumafcpp;

using namespace af;
typedef std::unique_ptr<FieldTerm> llgt_ptr;

void calc_mean_m(const State& state, std::ostream& myfile, double hzee) {
    array sum_dim3 = sum(sum(sum(state.m, 0), 1), 2);
    myfile << std::setw(12) << state.t << "\t" << afvalue(sum_dim3(span, span, span, 0)) << "\t"
           << afvalue(sum_dim3(span, span, span, 1)) << "\t" << afvalue(sum_dim3(span, span, span, 2)) << "\t" << hzee
           << std::endl;
}

double rate = 0.10e6; //[T/s]
double hzee_max = 2;  //[T]

// af::array zee_func(State state){
//    double field_Tesla = 0;
//    double rate = hzee_max/quater_steps; //[T/s]
//    if(state.t < hzee_max/rate) field_Tesla = rate *state.t;
//    else if (state.t < 3*hzee_max/rate) field_Tesla = -rate *state.t +
//    2*hzee_max; else if(state.t < 4*hzee_max/rate) field_Tesla = rate*state.t
//    - 4*hzee_max; else {field_Tesla = 0; std::cout << "WARNING ZEE time out of
//    range" << std::endl;} array zee = constant(0.0, state.mesh.nx,
//    state.mesh.ny, state.mesh.nz, 3, f64); zee(span, span, span,
//    0)=constant(field_Tesla/state.constants::mu0 , state.mesh.nx,
//    state.mesh.ny, state.mesh.nz, 1, f64); return  zee;
//}

af::array zee_func(State state) {
    double field_Tesla = 0;
    if (state.t < hzee_max / rate)
        field_Tesla = rate * state.t;
    else if (state.t < 3 * hzee_max / rate)
        field_Tesla = -rate * state.t + 2 * hzee_max;
    else if (state.t < 4 * hzee_max / rate)
        field_Tesla = rate * state.t - 4 * hzee_max;
    else {
        field_Tesla = 0;
        std::cout << "WARNING ZEE time out of range" << std::endl;
    }
    array zee = constant(0.0, state.mesh.nx, state.mesh.ny, state.mesh.nz, 3, f64);
    zee(span, span, span, 0) =
        constant(field_Tesla / state.constants::mu0, state.mesh.nx, state.mesh.ny, state.mesh.nz, 1, f64);
    return zee;
}

int main(int argc, char** argv) {

    std::cout << "argc = " << argc << std::endl;
    for (int i = 0; i < argc; i++)
        cout << "Parameter " << i << " was " << argv[i] << "\n";

    std::string filepath(argc > 1 ? argv[1] : "../Data/skyrmion_stoch");
    if (argc > 0)
        filepath.append("/");
    std::cout << "Writing into path " << filepath.c_str() << std::endl;

    setDevice(argc > 2 ? std::stoi(argv[2]) : 0);
    info();

    // Parameter initialization
    double length = 90e-9; //[nm]
    const double dx = 0.5e-9;
    const int nx = (int)(length / dx);
    std::cout << "nx = " << nx << std::endl;

    // Generating Objects
    Mesh mesh(nx, nx, 1, dx, dx, dx);
    Material material = Material();
    state.Ms = 580000;
    material.A = 15e-12;
    material.alpha = 1;
    material.D = 3e-3;
    material.Ku1 = 0.6e6;

    material.J_atom = 2. * material.A * dx;
    material.D_atom = material.D * pow(dx, 2);
    material.K_atom = material.Ku1 * pow(dx, 3);
    material.p = state.Ms * pow(dx, 3); // Compensate nz=1 instead of nz=4

    // Initial magnetic field
    array m = constant(0.0, mesh.nx, mesh.ny, mesh.nz, 3, f64);
    m(span, span, span, 2) = -1;
    for (int ix = 0; ix < mesh.nx; ix++) {
        for (int iy = 0; iy < mesh.ny; iy++) {
            const double rx = double(ix) - mesh.nx / 2.;
            const double ry = double(iy) - mesh.ny / 2.;
            const double r = sqrt(pow(rx, 2) + pow(ry, 2));
            if (r > nx / 4.)
                m(ix, iy, span, 2) = 1.;
        }
    }

    State state(mesh, material, m);
    vti_writer_atom(state.m, mesh, (filepath + "minit").c_str());

    std::vector<llgt_ptr> llgterm;
    llgterm.push_back(llgt_ptr(new AtomisticDipoleDipoleField(mesh)));
    llgterm.push_back(llgt_ptr(new AtomisticExchangeField(mesh)));
    llgterm.push_back(llgt_ptr(new AtomisticDmiField(mesh, material)));
    llgterm.push_back(llgt_ptr(new AtomisticUniaxialAnisotropyField(mesh, material)));

    LLG Llg(state, llgterm);

    timer t = af::timer::start();
    double E_prev = 1e20;
    while (fabs((E_prev - Llg.E(state)) / E_prev) > 1e-10) {
        E_prev = Llg.E(state);
        for (int i = 0; i < 100; i++) {
            state.m = Llg.step(state);
        }
        if (state.steps % 1000 == 0)
            std::cout << "step " << state.steps << " rdiff= " << fabs((E_prev - Llg.E(state)) / E_prev) << std::endl;
    }
    std::cout << "time =" << state.t << " [s], E = " << Llg.E(state) << "[J]" << std::endl;
    std::cout << "timerelax [af-s]: " << af::timer::stop(t) << ", steps = " << state.steps << std::endl;
    vti_writer_atom(state.m, mesh, (filepath + "relax").c_str());

    // State state(mesh, material, m);
    // vti_writer_atom(state.m, mesh , (filepath + "minit").c_str());

    //// Relax
    // af::timer timer_llgterms = af::timer::start();
    // Minimizer minimizer("BB", 1e-10, 1e-5, 1e4, 100);
    ////minimizer.llgterms.push_back( uptr_Fieldterm (new
    /// AtomisticDipoleDipoleField(mesh)));
    // minimizer.llgterms.push_back( uptr_Fieldterm (new
    // AtomisticExchangeField(mesh))); minimizer.llgterms.push_back( uptr_Fieldterm
    // (new AtomisticDmiField(mesh, material))); minimizer.llgterms.push_back(
    // uptr_Fieldterm (new AtomisticUniaxialAnisotropyField(mesh, material)));
    // std::cout<<"Llgterms assembled in "<< af::timer::stop(timer_llgterms)
    // <<std::endl;

    ////obtaining relaxed magnetization
    // timer t = af::timer::start();
    // minimizer.minimize(state);
    // std::cout<<"timerelax [af-s]: "<< af::timer::stop(t) <<std::endl;
    // vti_writer_micro(state.m, mesh , (filepath + "relax").c_str());

    // Hysteresis
    std::ofstream stream;
    stream.precision(12);
    stream.open((filepath + "m.dat").c_str());
    stream << "# t	<mx>    <my>    <mz>    hzee" << std::endl;
    calc_mean_m(state, stream, 0);

    // timer t_hys = af::timer::start();
    // double rate = hzee_max/quater_steps; //[T/s]
    // minimizer.llgterms.push_back( uptr_Fieldterm (new ExternalField(&zee_func)));
    // while (state.t < 4* hzee_max/rate){
    //    minimizer.minimize(state);
    //    calc_mean_m(state, stream, afvalue(minimizer.llgterms[4]->h(state)(0,
    //    0, 0, 0))); state.t+=1.; state.steps++; if( state.steps % 1 == 0){
    //        vti_writer_micro(state.m, mesh , (filepath +
    //        "m_hysteresis_"+std::to_string(state.steps)).c_str());
    //    }
    //}
    // std::cout<<"time full hysteresis [af-s]: "<< af::timer::stop(t_hys)
    // <<std::endl;
    timer t_hys = af::timer::start();
    Llg.Fieldterms.push_back(llgt_ptr(new ExternalField(&zee_func))); // Rate in T/s
    while (state.t < 4 * hzee_max / rate) {
        state.m = Llg.step(state);
        if (state.steps % 2000 == 0) {
            calc_mean_m(state, stream, afvalue(Llg.Fieldterms[4]->h(state)(0, 0, 0, 0)));
            vti_writer_micro(state.m, mesh, (filepath + "m_hysteresis_" + std::to_string(state.steps)).c_str());
        }
        // std::cout << state.t << "\t" <<
        // afvalue(sum(sum(sum(sum(Llg.Fieldterms[2]->h(state), 3), 2), 1),
        // 0))<< std::endl;
    }
    stream.close();
    std::cout << "time full hysteresis [af-s]: " << af::timer::stop(t_hys) << std::endl;
    return 0;
}
