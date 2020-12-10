// layer 0 is magnetic, layer nz 1 is empty
#include "arrayfire.h"
#include "magnum_af.hpp"

#if __has_include(<filesystem>)
#include <filesystem>
#define have_filesystem 1
#elif __has_include(<experimental/filesystem>)
#include <experimental/filesystem>
#define have_filesystem 1
#define experimental_filesystem 1
#else
#define have_filesystem 0
#endif

#if __cplusplus < 201703L // If the version of C++ is less than 17
// It was still in the experimental:: namespace
namespace fs = std::experimental::filesystem;
#else
namespace fs = std::filesystem;
#endif
// https://en.cppreference.com/w/cpp/preprocessor/include
// https://stackoverflow.com/questions/20358455/cross-platform-way-to-make-a-directory
//#include <experimental/filesystem>

using namespace magnumafcpp;

using namespace af;

int main(int argc, char** argv) {

    std::cout << "argc = " << argc << std::endl;
    for (int i = 0; i < argc; i++)
        std::cout << "Parameter " << i << " was " << argv[i] << "\n";
    std::string filepath(argc > 1 ? argv[1] : "./run/");
    if (argc > 1)
        filepath.append("/");
    std::cout << "Writing into path " << filepath << std::endl;

    setDevice(argc > 2 ? std::stoi(argv[2]) : 0);
    info();
    StageTimer timer;

    std::ofstream stream(filepath + "m_and_E.dat");
    stream << "#skrym: steps, time, mz, E, exch , aniso , dmi , ext";
    stream << "#ferro: steps, time, mz, E, exch , aniso , dmi , ext" << std::endl;
    for (unsigned nz = 1; nz < 20; nz++) {
        // Parameter initialization
        const unsigned nx = 100, ny = 100;
        // const int nz(argc>3? std::stoi(argv[3]):2);
        std::cout << "Running for nz=" << nz << std::endl;
        const double x = 400e-9;
        const double y = 400e-9;
        const double z = nz * 1e-9;

        const double dx = x / nx;
        const double dy = y / ny;
        const double dz = z / nz;

        const double Hz = 130e-3 / constants::mu0;
        // const double RKKY_val = 0.8e-3 * 1e-9* 0.5;

        // SK layer params
        const double Ms = 1371e3;  // A/m
        const double A = 15e-12;   // J/m
        const double Ku = 1.411e6; // J/m^3
        const double D = 2.5e-3;   // J/m^2

        // Generating Objects
        Mesh mesh(nx, ny, nz, dx, dy, dz);

        // Initial magnetic field
        array m = constant(0.0, mesh.nx, mesh.ny, mesh.nz, 3, f64);
        for (unsigned ix = 0; ix < mesh.nx; ix++) {
            for (unsigned iy = 0; iy < mesh.ny; iy++) {
                const double rx = double(ix) - mesh.nx / 2.;
                const double ry = double(iy) - mesh.ny / 2.;
                const double r = sqrt(pow(rx, 2) + pow(ry, 2));
                if (r > nx / 4.) {
                    m(ix, iy, af::span, 2) = 1.;
                } else {
                    m(ix, iy, af::span, 2) = -1.;
                }
            }
        }

        // defining interactions
        auto exch = uptr_FieldTerm(new ExchangeField(A));
        auto aniso = uptr_FieldTerm(new UniaxialAnisotropyField(Ku, (std::array<double, 3>){0, 0, 1}));
        auto dmi = uptr_FieldTerm(new DmiField(D, {0, 0, -1}));
        af::array zee = af::constant(0.0, mesh.nx, mesh.ny, mesh.nz, 3, f64);
        zee(af::span, af::span, af::span, 2) = Hz;
        auto external = uptr_FieldTerm(new ExternalField(zee));

        // af::print("dmi", dmi->h(state_1));
        // af::print("exch", exch->h(state_1));

        LLGIntegrator llg(1, {exch, aniso, dmi, external});
        timer.print_stage("init ");

        State state_1(mesh, Ms, m);

        // llg.relax(state_1);
        while (state_1.t < 3e-9) {
            llg.step(state_1);
        }
        timer.print_stage("relax");
        state_1.write_vti(filepath + "m_relaxed_nz_" + std::to_string(nz));
        std::cout << "skrmy relaxed after [s]" << state_1.t << std::endl;

        af::array m2 = af::constant(0.0, mesh.nx, mesh.ny, mesh.nz, 3, f64);
        m2(af::span, af::span, af::span, 2) = 1;

        State state_2(mesh, Ms, m2);

        state_2.write_vti(filepath + "m_state_2_init_nz_" + std::to_string(nz));
        // llg.relax(state_2);
        while (state_2.t < 3e-9) {
            llg.step(state_2);
        }
        timer.print_stage("relax state2");

        state_2.write_vti(filepath + "m_state_2_relaxed_nz_" + std::to_string(nz));
        timer.print_stage("state2 relaxed");
        std::cout << "ferro relaxed after [s]" << state_2.t << std::endl;

        std::array<std::string, 5> llgnames = {"exch ", "aniso", "dmi  ", "ext  "};

        std::cout << std::scientific << "Skrym, " << nz << ", t=" << state_1.t << ", mz=" << state_1.meani(2)
                  << ", E=" << llg.E(state_1);
        for (unsigned i = 0; i < llg.llgterms.size(); i++) {
            std::cout << ", " << llgnames[i] << "=" << llg.llgterms[i]->Energy_in_J(state_1);
        }
        std::cout << std::endl;

        std::cout << std::scientific << "Ferro, " << nz << ", t=" << state_2.t << ", mz=" << state_2.meani(2)
                  << ", E=" << llg.E(state_2);
        for (unsigned i = 0; i < llg.llgterms.size(); i++) {
            std::cout << ", " << llgnames[i] << "=" << llg.llgterms[i]->Energy_in_J(state_2);
        }
        std::cout << std::endl;

        stream << std::scientific << nz << "\t" << state_1.t << "\t" << state_1.meani(2) << "\t" << llg.E(state_1);
        for (unsigned i = 0; i < llg.llgterms.size(); i++) {
            stream << "\t" << llg.llgterms[i]->Energy_in_J(state_1);
        }

        stream << std::scientific << "\t" << nz << "\t" << state_2.t << "\t" << state_2.meani(2) << "\t"
               << llg.E(state_2);
        for (unsigned i = 0; i < llg.llgterms.size(); i++) {
            stream << "\t" << llg.llgterms[i]->Energy_in_J(state_2);
        }
        stream << std::endl;
        timer.print_stage(std::to_string(nz) + " stateenergy");

        // StringMethod method
        double n_interp = 60;
        double string_dt = 1e-13;
        std::vector<State> inputimages;
        inputimages.push_back(state_1);
        inputimages.push_back(state_2);

        StringMethod string(state_1, inputimages, n_interp, string_dt, llg);
        fs::create_directory(filepath + std::to_string(nz) + "/");
        string.run(filepath + std::to_string(nz) + "/");
    }

    stream.close();
    timer.print_accumulated();

    return 0;
}
