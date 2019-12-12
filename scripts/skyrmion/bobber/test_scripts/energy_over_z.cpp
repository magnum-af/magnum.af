#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace magnumafcpp;

using namespace af;

int main(int argc, char **argv)
{

    std::cout << "argc = " << argc << std::endl;
    for (int i = 0; i < argc; i++)
        std::cout << "Parameter " << i << " was " << argv[i] << "\n";
    std::string filepath(argc > 1 ? argv[1] : "./run/");
    if (argc > 1)
        filepath.append("/");
    std::cout << "Writing into path " << filepath.c_str() << std::endl;

    setDevice(argc > 2 ? std::stoi(argv[2]) : 0);
    info();
    StageTimer timer;
    std::array<std::string, 5> llgnames = {"demag", "exch ", "aniso", "dmi  ", "ext  "};

    // Calculating m_relaxed.vti
    af::array m_relaxed;
    {
        // Parameter initialization
        const int nx = 100, ny = 100, nz = 1;
        const double x = 400e-9;
        const double y = 400e-9;
        const double z = nz * 1e-9;

        const double dx = x / nx;
        const double dy = y / ny;
        const double dz = z / nz;

        const double Hz = 130e-3 / constants::mu0;
        //const double RKKY_val = 0.8e-3 * 1e-9* 0.5;

        // SK layer params
        const double Ms = 1371e3;  // A/m
        const double A = 15e-12;   // J/m
        const double Ku = 1.411e6; // J/m^3
        const double D = 2.5e-3;   // J/m^2

        //Generating Objects
        Mesh mesh(nx, ny, nz, dx, dy, dz);

        // Initial magnetic field
        array m = constant(0.0, mesh.n0, mesh.n1, mesh.n2, 3, f64);
        m(af::span, af::span, af::span, 2) = -1;
        for (int ix = 0; ix < mesh.n0; ix++)
        {
            for (int iy = 0; iy < mesh.n1; iy++)
            {
                const double rx = double(ix) - mesh.n0 / 2.;
                const double ry = double(iy) - mesh.n1 / 2.;
                const double r = sqrt(pow(rx, 2) + pow(ry, 2));
                if (r > nx / 4.)
                    m(ix, iy, af::span, 2) = 1.;
            }
        }

        // defining interactions
        auto demag = LlgTerm(new DemagField(mesh, true, true, 0));
        auto exch = LlgTerm(new ExchangeField(A));
        auto aniso = LlgTerm(new UniaxialAnisotropyField(Ku, (std::array<double, 3>){0, 0, 1}));
        auto dmi = LlgTerm(new DmiField(D, {0, 0, -1}));
        af::array zee = af::constant(0.0, mesh.n0, mesh.n1, mesh.n2, 3, f64);
        zee(af::span, af::span, af::span, 2) = Hz;
        auto external = LlgTerm(new ExternalField(zee));

        LLGIntegrator llg(1, {demag, exch, aniso, dmi, external});
        timer.print_stage("init ");

        State state_1(mesh, Ms, m);

        if (!exists(filepath + "m_relaxed.vti"))
        {
            std::cout << "Relaxing minit" << std::endl;
            state_1.write_vti(filepath + "minit");
            //LLGIntegrator Llg(1, {demag, exch, aniso, dmi, external});
            while (state_1.t < 3e-9)
            {
                if (state_1.steps % 100 == 0)
                    state_1.write_vti(filepath + "m_step" + std::to_string(state_1.steps));
                llg.step(state_1);
                std::cout << std::scientific << state_1.steps << "\t" << state_1.t << "\t" << state_1.meani(2) << "\t" << llg.E(state_1) << std::endl;
            }
            //Llg.relax(state_1);
            timer.print_stage("relax");
            state_1.write_vti(filepath + "m_relaxed");
        }
        else
        {
            std::cout << "Found m_relaxed, reading in state_1." << std::endl;
            state_1._vti_reader(filepath + "m_relaxed.vti");
            state_1.write_vti(filepath + "m_relaxed_from_read_in");
        }
        m_relaxed = state_1.m;
    }

    // Looping over nz
    std::ofstream stream(filepath + "m_and_E.dat");
    stream << "#steps, time, mz, E, demag , exch , aniso , dmi , ext" << std::endl;
    for (int nz = 1; nz < 20; nz++)
    {
        // Parameter initialization
        const int nx = 100, ny = 100;
        const double x = 400e-9;
        const double y = 400e-9;
        const double z = nz * 1e-9;

        const double dx = x / nx;
        const double dy = y / ny;
        const double dz = z / nz;

        const double Hz = 130e-3 / constants::mu0;
        //const double RKKY_val = 0.8e-3 * 1e-9* 0.5;

        // SK layer params
        const double Ms = 1371e3;      // A/m
        const double val_A = 15e-12;   // J/m
        const double val_Ku = 1.411e6; // J/m^3
        const double val_D = 2.5e-3;   // J/m^2
        array A = constant(0.0, nx, ny, nz, 3, f64);
        A(af::span, af::span, 0, af::span) = val_A;

        array Ku = constant(0.0, nx, ny, nz, 3, f64);
        Ku(af::span, af::span, 0, af::span) = val_Ku;

        array D = constant(0.0, nx, ny, nz, 3, f64);
        D(af::span, af::span, 0, af::span) = val_D;

        //Generating Objects
        Mesh mesh(nx, ny, nz, dx, dy, dz);

        // Initial magnetic field
        af::array m = constant(0.0, mesh.n0, mesh.n1, mesh.n2, 3, f64);
        m(af::span, af::span, 0, af::span) = m_relaxed;
        State state_1(mesh, Ms, m);
        state_1.write_vti(filepath + "m_relaxed_nz" + std::to_string(nz));

        // defining interactions
        auto demag = LlgTerm(new DemagField(mesh, true, true, 0));
        auto exch = LlgTerm(new ExchangeField(A));
        //auto exch = LlgTerm (new SparseExchangeField(A, mesh));
        auto aniso = LlgTerm(new UniaxialAnisotropyField(Ku, (std::array<double, 3>){0, 0, 1}));
        auto dmi = LlgTerm(new DmiField(D, {0, 0, -1}));
        af::array zee = af::constant(0.0, mesh.n0, mesh.n1, mesh.n2, 3, f64);
        zee(af::span, af::span, af::span, 2) = Hz;
        auto external = LlgTerm(new ExternalField(zee));

        LLGIntegrator llg(1, {demag, exch, aniso, dmi, external});
        timer.print_stage("init ");

        std::cout << std::scientific << nz << ", t=" << state_1.t << ", mz=" << state_1.meani(2) << ", E=" << llg.E(state_1);
        for (unsigned i = 0; i < llg.llgterms.size(); i++)
        {
            std::cout << ", " << llgnames[i] << "=" << llg.llgterms[i]->E(state_1);
        }
        std::cout << std::endl;

        stream << std::scientific << nz << "\t" << state_1.t << "\t" << state_1.meani(2) << "\t" << llg.E(state_1);
        for (unsigned i = 0; i < llg.llgterms.size(); i++)
        {
            stream << "\t" << llg.llgterms[i]->E(state_1);
        }
        stream << std::endl;
    }

    return 0;
}
