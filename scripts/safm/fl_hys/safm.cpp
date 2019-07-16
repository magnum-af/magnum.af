#include "magnum_af.hpp"

using namespace magnumaf;

void print_m(const af::array& m, std::string name, std::ostream& stream = std::cout){
    af::array m_mean = af::mean(af::mean(af::mean(m, 0), 1), 2);
    stream << name << ": " << afvalue(m_mean(0, 0, 0, 0)) << ", " << afvalue(m_mean(0, 0, 0, 1)) << ", " << afvalue(m_mean(0, 0, 0, 2)) <<  std::endl;
}

int main(int argc, char** argv)
{
    // Checking input variables and setting GPU Device
    StageTimer timer;
    for (int i=0; i<argc; i++){std::cout << "Parameter " << i << " was " << argv[i] << std::endl;}
    std::string filepath = std::string(argc>1? argv[1]: "output_magnum.af/");
    af::setDevice(argc>2? std::stoi(argv[2]):0);
    af::info();

    //std::cout.precision(18);

    // Parameter initialization
    const double x=40.e-9, y=40e-9, z=4.e-9;
    const int nx = 32, ny=32 , nz=4;
    const double Ms = 8e5;//TODO 0.5/Constants.mu0
    const double A = 15e-12;
    const double Ku1_safm=1e6; // Anisotropy constant
    const double Ku1_freelayer=79e+03; // Anisotropy constant

    //af::timer timer = af::timer::start();
    //Generating Objects
    std::vector<double> z_spacing = {z/nz, 9.84825e-10, z/nz, z/nz};
    NonequispacedMesh mesh(nx, ny, x/nx, y/ny, z_spacing);
    mesh.print();

    // Initial magnetic field
    af::array m = af::constant(0, mesh.dims, f64);
    m(af::span, af::span, 0, 2) = 1;//SAFM Layer 0 in z
    m(af::span, af::span, 1, 2) = -1;//SAFM Layer 1 in -z
    m(af::span, af::span, 2, af::span) = 0;//Vacuum Layer 2
    m(af::span, af::span, 3, 0) = 1;// Free Layer in x

    af::array Ku1_field = af::constant(0, mesh.dims, f64);
    Ku1_field(af::span, af::span, 0, af::span) = Ku1_safm;//SAFM Layer 0 in z
    Ku1_field(af::span, af::span, 1, af::span) = Ku1_safm;//SAFM Layer 1 in -z
    Ku1_field(af::span, af::span, 2, af::span) = 0;//Vacuum Layer 2
    Ku1_field(af::span, af::span, 3, af::span) = Ku1_freelayer;// Free Layer in x
    State state(mesh, Ms, m, false);
    state.mesh = Mesh(nx, ny, nz, x/nx, y/ny, z/nz);
    vtr_writer(state.m, Mesh(nx, ny, nz, x/nx, y/ny, 0) , z_spacing, filepath + "m_init");
    timer.print_stage("mesh  ");

    auto demag = LlgTerm (new NonEquiDemagField(mesh, true, false, 0));
    auto exch = LlgTerm (new NonequiExchangeField(A, mesh, true));
    auto aniso = LlgTerm (new UniaxialAnisotropyField(Ku1_field, std::array<double, 3>{0, 0, 1}));//TODO energy is wrong is its not nonequi!
    LlgTerms llgterms = {demag, exch, aniso};
    timer.print_stage("setup ");

    af::array h = demag->h(state);
    vtr_writer(h, Mesh(nx, ny, nz, x/nx, y/ny, 0) , z_spacing, filepath + "h");

    std::ofstream stream;
    stream.precision(12);
    stream.open (filepath + "m.dat");
    stream << z_spacing[1] << ", " << afvalue(h(nx/2, ny/2, 3, 0)) << ", " << afvalue(h(nx/2, ny/2, 3, 1)) << ", " << afvalue(h(nx/2, ny/2, 3, 2)) << std::endl;
    stream.close();
    timer.print_stage("calc h");

    LLGIntegrator llg(1, llgterms);
    const double time_relax = 1e-11;
    while (state.t < time_relax){
        llg.step(state);
        //state.calc_mean_m(std::cout);
        print_m(state.m(af::span, af::span, 0, af::span), "safm lay0");
        print_m(state.m(af::span, af::span, 1, af::span), "safm lay1");
        print_m(state.m(af::span, af::span, 2, af::span), "vac  lay2");
        print_m(state.m(af::span, af::span, 3, af::span), "freelayer");
        std::cout << std::endl;
    }
    timer.print_stage("relax " + std::to_string(time_relax * 1e9) + " [ns]");
    vtr_writer(state.m, Mesh(nx, ny, nz, x/nx, y/ny, 0) , z_spacing, filepath + "m_relax");
    //llg.relax(state, 0.01);
    

    timer.print_accumulated();
    return 0;
}
