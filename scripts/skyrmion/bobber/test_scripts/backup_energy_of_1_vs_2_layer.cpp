#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace magnumafcpp;

using namespace af;

int main(int argc, char** argv)
{

    std::cout<<"argc = "<<argc<<std::endl;
    for (int i=0; i<argc; i++) std::cout << "Parameter " << i << " was " << argv[i] << "\n";
    std::string filepath(argc>1? argv[1]: "./run/");
    if(argc>1)filepath.append("/");
    std::cout<<"Writing into path "<<filepath.c_str()<<std::endl;

    setDevice(argc>2? std::stoi(argv[2]):0);
    info();
    StageTimer timer;

    double E_skrym_1layer;
    double E_ferro_1layer;
    af::array m_skyrm1layer;
    af::array m_ferro1layer;
    std::vector<double> e_1layer_state_1;
    std::vector<double> e_1layer_state_2;
    std::vector<double> e_2layer_state_1;
    std::vector<double> e_2layer_state_2;
    std::array<std::string, 5> llgnames = {"demag","exch ","aniso","dmi  ","ext  "};

    {
        // Parameter initialization
        const int nx = 100, ny=100 , nz=1;
        const double x=400e-9;
        const double y=400e-9;
        const double z= nz * 1e-9;

        const double dx= x/nx;
        const double dy= y/ny;
        const double dz= z/nz;

        const double Hz = 130e-3/constants::mu0;
        //const double RKKY_val = 0.8e-3 * 1e-9* 0.5;

        // SK layer params
        const double Ms =1371e3;// A/m
        const double A = 15e-12;// J/m
        const double Ku = 1.411e6;// J/m^3
        const double D =2.5e-3;// J/m^2

        //Generating Objects
        Mesh mesh(nx, ny, nz, dx, dy, dz);

        // Initial magnetic field
        array m = constant(0.0, mesh.n0, mesh.n1, mesh.n2, 3, f64);
        m(af::span, af::span, af::span, 2) = -1;
        for(int ix=0;ix<mesh.n0;ix++){
            for(int iy=0;iy<mesh.n1;iy++){
                const double rx=double(ix)-mesh.n0/2.;
                const double ry=double(iy)-mesh.n1/2.;
                const double r = sqrt(pow(rx, 2)+pow(ry, 2));
                if(r>nx/4.) m(ix, iy, af::span, 2)=1.;
            }
        }

        // defining interactions
        auto demag = LlgTerm (new DemagField(mesh, true, true, 0));
        auto exch = LlgTerm (new SparseExchangeField(A, mesh));
        //TODO//becomes zero//auto exch = LlgTerm (new ExchangeField(A));
        auto aniso = LlgTerm (new UniaxialAnisotropyField(Ku, (std::array<double ,3>) {0, 0, 1}));
        auto dmi = LlgTerm (new DmiField(D, {0, 0, -1}));
        af::array zee = af::constant(0.0, mesh.n0, mesh.n1, mesh.n2, 3, f64);
        zee(af::span, af::span, af::span, 2) = Hz;
        auto external = LlgTerm (new ExternalField(zee));

        //af::print("dmi", dmi->h(state_1));
        //af::print("exch", exch->h(state_1));

        LLGIntegrator llg(1, {demag, exch, aniso, dmi, external});
        timer.print_stage("init ");

        State state_1(mesh, Ms, m);

        if (! exists(filepath + "m_relaxed.vti" ) ){
            std::cout << "Relaxing minit" << std::endl; state_1.write_vti(filepath + "minit");
            //LLGIntegrator Llg(1, {demag, exch, aniso, dmi, external});
            while (state_1.t < 0.5e-9){
                if (state_1.steps % 100 == 0) state_1.write_vti(filepath + "m_step" + std::to_string(state_1.steps));
                llg.step(state_1);
                std::cout << std::scientific << state_1.steps << "\t" << state_1.t << "\t" <<  state_1.meani(2) << "\t" << llg.E(state_1) << std::endl;
            }
            //Llg.relax(state_1);
            timer.print_stage("relax");
            state_1.write_vti(filepath + "m_relaxed");
        }
        else {
            std::cout << "Found m_relaxed, reading in state_1." << std::endl;
            state_1._vti_reader(filepath + "m_relaxed.vti");
            state_1.write_vti(filepath + "m_relaxed_from_read_in");

        }

        af::array m2 = af::constant(0.0, mesh.n0, mesh.n1, mesh.n2, 3, f64);
        m2(af::span, af::span, af::span, 2) = 1;

        State state_2(mesh, Ms, m2);

        vti_writer_micro(state_2.m, state_2.mesh, filepath + "m_state_2_init");
        m_skyrm1layer = state_1.m;
        m_ferro1layer = state_2.m;

        std::cout.precision(18);
        E_skrym_1layer = llg.E(state_1);
        E_ferro_1layer = llg.E(state_2);

        //for(auto const& value: llg.llgterms) {
        for(unsigned i = 0; i < llg.llgterms.size(); i++) {
            e_1layer_state_1.push_back(llg.llgterms[i]->E(state_1));
            std::cout << "1layer E(state_1) " << llgnames[i] << " =" << llg.llgterms[i]->E(state_1) << std::endl;
        }
        std::cout << std::endl;
        //for(auto const& value: llg.llgterms) {
        for(unsigned i = 0; i < llg.llgterms.size(); i++) {
            e_1layer_state_2.push_back(llg.llgterms[i]->E(state_2));
            std::cout << "1layer E(state_2) " << llgnames[i] << " =" << llg.llgterms[i]->E(state_2) << std::endl;
        }
        std::cout << std::endl;
    }

    double E_skrym_2layer;
    double E_ferro_2layer;

    {
        // Parameter initialization
        const int nx = 100, ny=100 , nz=2;
        const double x=400e-9;
        const double y=400e-9;
        const double z= nz * 1e-9;

        const double dx= x/nx;
        const double dy= y/ny;
        const double dz= z/nz;

        const double Hz = 130e-3/constants::mu0;
        //const double RKKY_val = 0.8e-3 * 1e-9* 0.5;

        // SK layer params
        const double SK_Ms =1371e3;// A/m
        const double SK_A = 15e-12;// J/m
        const double SK_Ku = 1.411e6;// J/m^3
        const double SK_D =2.5e-3;// J/m^2

        array geom = af::constant(0.0, nx, ny, nz, 3, f64);
        geom(af::span, af::span, 0, af::span) = 1;
        af::array todouble = af::constant(1., nx, ny, nz, 3, f64);
        af::array Ms = SK_Ms * (geom == 1) * todouble;
        af::array A  = SK_A  * (geom == 1) * todouble;
        af::array Ku = SK_Ku * (geom == 1) * todouble;
        af::array D  = SK_D  * (geom == 1) * todouble; // + IL_D * (geom == 2);

        //Generating Objects
        Mesh mesh(nx, ny, nz, dx, dy, dz);

        // defining interactions
        auto demag = LlgTerm (new DemagField(mesh, true, true, 0));
        auto exch = LlgTerm (new SparseExchangeField(A, mesh));
        //TODO//becomes zero//auto exch = LlgTerm (new ExchangeField(A));
        auto aniso = LlgTerm (new UniaxialAnisotropyField(Ku, (std::array<double ,3>) {0, 0, 1}));
        auto dmi = LlgTerm (new DmiField(D, {0, 0, -1}));
        af::array zee = af::constant(0.0, mesh.n0, mesh.n1, mesh.n2, 3, f64);
        zee(af::span, af::span, af::span, 2) = Hz;
        auto external = LlgTerm (new ExternalField(zee));

        //af::print("dmi", dmi->h(state_1));
        //af::print("exch", exch->h(state_1));

        LLGIntegrator llg(1, {demag, exch, aniso, dmi, external});
        timer.print_stage("init ");

        af::array m_state_2layer = af::constant(0, mesh.dims, f64);
        m_state_2layer(af::span, af::span, 0, af::span) = m_skyrm1layer;
        State state_1(mesh, Ms, m_state_2layer);

        af::array m2 = af::constant(0.0, mesh.n0, mesh.n1, mesh.n2, 3, f64);
        m2(af::span, af::span, 0, 2) = 1;

        State state_2(mesh, Ms, m2);

        vti_writer_micro(state_2.m, state_2.mesh, filepath + "2layer_m_state_2_init");

        std::cout.precision(18);
        E_skrym_2layer = llg.E(state_1);
        E_ferro_2layer = llg.E(state_2);

        for(unsigned i = 0; i < llg.llgterms.size(); i++) {
            e_2layer_state_1.push_back(llg.llgterms[i]->E(state_1));
            std::cout << "2layer E(state_1) " << llgnames[i] << " =" << llg.llgterms[i]->E(state_1) << std::endl;
        }
        std::cout << std::endl;
        for(unsigned i = 0; i < llg.llgterms.size(); i++) {
            e_2layer_state_2.push_back(llg.llgterms[i]->E(state_2));
            std::cout << "2layer E(state_2) " << llgnames[i] << " =" << llg.llgterms[i]->E(state_2) << std::endl;
        }
        //for(auto const& value: llg.llgterms) {
        //    e_2layer_state_1.push_back(value->E(state_1));
        //    std::cout << "2layer E(state_1)=" << value->E(state_1) << std::endl;
        //}
        //std::cout << std::endl;

        //for(auto const& value: llg.llgterms) {
        //    e_2layer_state_2.push_back(value->E(state_2));
        //    std::cout << "2layer E(state_2)=" << value->E(state_2) << std::endl;
        //}
        std::cout << std::endl;
    }

    std::cout << std::scientific;
    std::cout << "E skyrm 1layer=" << E_skrym_1layer << std::endl;
    std::cout << "E ferro 1layer=" << E_ferro_1layer << std::endl;
    std::cout << "E diff  1layer=" << E_skrym_1layer - E_ferro_1layer << std::endl;

    std::cout << "E skyrm 2layer=" << E_skrym_2layer << std::endl;
    std::cout << "E ferro 2layer=" << E_ferro_2layer << std::endl;
    std::cout << "E diff  2layer=" << E_skrym_2layer - E_ferro_2layer << std::endl;
    std::cout << std::endl;

    std::cout << "E diff skyrm 12=" << E_skrym_1layer - E_skrym_2layer<< std::endl;
    std::cout << std::endl;

    for(unsigned i = 0; i<e_1layer_state_1.size(); i++){
        std::cout << "diff " << llgnames[i]<< " skyrm["<<i<<"]=" << e_1layer_state_1[i] - e_2layer_state_1[i] << std::endl;
        std::cout << "diff " << llgnames[i]<< " ferro["<<i<<"]=" << e_1layer_state_2[i] - e_2layer_state_2[i] << std::endl;
        //std::cout << "e_1layer_state_1[i]" << e_1layer_state_1[i] << std::endl;
        //std::cout << "e_1layer_state_2[i]" << e_1layer_state_2[i] << std::endl;
        //std::cout << "e_2layer_state_1[i]" << e_2layer_state_1[i] << std::endl;
        //std::cout << "e_2layer_state_2[i]" << e_2layer_state_2[i] << std::endl;
    }
    //TODO check RKKYExchangeField in this setup?

    return 0;
}
