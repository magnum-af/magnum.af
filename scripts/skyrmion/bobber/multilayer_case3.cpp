#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace magnumafcpp;

using namespace af;

int main(int argc, char** argv) {

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

    // Parameter initialization
    const double x = 400e-9;
    const double y = 400e-9;
    const double z = 32e-9;

    const int nx = 100, ny = 100, nz = 32;
    const double dx = x / nx;
    const double dy = y / ny;
    const double dz = z / nz;

    const double Hz = 130e-3 / constants::mu0;
    const double RKKY_val = 0.8e-3 * 1e-9; // TODO
    // Note: maybe 0.5 factor (see mumax3):
    // const double RKKY_val = 0.8e-3 * 1e-9* 0.5;//

    // SK layer params
    const double SK_Ms = 1371e3;  // A/m
    const double SK_A = 15e-12;   // J/m
    const double SK_Ku = 1.411e6; // J/m^3
    const double SK_D = 2.5e-3;   // J/m^2

    // Ferrimagnetic interface layer params
    const double IL_Ms = 488.2e3; // A/m
    const double IL_A = 4e-12;    // J/m
    const double IL_Ku = 486.6e3; // J/m^3
    const double IL_D = 0.8e-3;   // J/m^2

    array geom = af::constant(0.0, nx, ny, nz, 3, f64);
    geom(af::span, af::span, 0, af::span) = 1;
    geom(af::span, af::span, 3, af::span) = 1;
    geom(af::span, af::span, 6, af::span) = 1;
    geom(af::span, af::span, 9, af::span) = 1;
    geom(af::span, af::span, 12, af::span) = 1;

    geom(af::span, af::span, 13, af::span) = 2;
    geom(af::span, af::span, 14, af::span) = 2;
    geom(af::span, af::span, 15, af::span) = 2;
    geom(af::span, af::span, 16, af::span) = 2;

    geom(af::span, af::span, 19, af::span) = 1;
    geom(af::span, af::span, 22, af::span) = 1;
    geom(af::span, af::span, 25, af::span) = 1;
    geom(af::span, af::span, 28, af::span) = 1;
    geom(af::span, af::span, 31, af::span) = 1;

    // af::print("geom", geom);
    // af::print("geom == 1", geom == 1);
    // af::print("geom == 2", geom == 2);

    af::array todouble = af::constant(1., nx, ny, nz, 3, f64);
    af::array Ms = (SK_Ms * (geom == 1) + IL_Ms * (geom == 2)) * todouble;
    af::array A = (SK_A * (geom == 1) + IL_A * (geom == 2)) * todouble;
    af::array Ku = (SK_Ku * (geom == 1) + IL_Ku * (geom == 2)) * todouble;
    af::array D = (SK_D * (geom == 1) + IL_D * (geom == 2)) * todouble; // + IL_D * (geom == 2);
    // std::cout << "type" << geom.type() << " " << Ms.type() << std::endl;

    // af::print("Ms", Ms);
    // af::print("1e12 * A", 1e12 * A);
    // af::print("Ku", Ku);
    // af::print("D", D);

    array RKKY = af::constant(0.0, nx, ny, nz, 3, f64);
    RKKY(af::span, af::span, 12, af::span) = RKKY_val;
    RKKY(af::span, af::span, 13, af::span) = RKKY_val;

    // Generating Objects
    Mesh mesh(nx, ny, nz, dx, dy, dz);

    // Initial magnetic field
    array m = constant(0.0, mesh.nx, mesh.ny, mesh.nz, 3, f64);
    m(af::span, af::span, af::span, 2) = -1;
    for (unsigned ix = 0; ix < mesh.nx; ix++) {
        for (unsigned iy = 0; iy < mesh.ny; iy++) {
            const double rx = double(ix) - mesh.nx / 2.;
            const double ry = double(iy) - mesh.ny / 2.;
            const double r = sqrt(pow(rx, 2) + pow(ry, 2));
            if (r > nx / 16.)
                m(ix, iy, af::span, 2) = 1.;
        }
    }
    m(af::span, af::span, 1, af::span) = 0;
    m(af::span, af::span, 2, af::span) = 0;

    m(af::span, af::span, 4, af::span) = 0;
    m(af::span, af::span, 5, af::span) = 0;

    m(af::span, af::span, 7, af::span) = 0;
    m(af::span, af::span, 8, af::span) = 0;

    m(af::span, af::span, 10, af::span) = 0;
    m(af::span, af::span, 11, af::span) = 0;

    // Setting IL layer ferromagnetic
    m(af::span, af::span, 13, 2) = 1.;
    m(af::span, af::span, 14, 2) = 1.;
    m(af::span, af::span, 15, 2) = 1.;
    m(af::span, af::span, 16, 2) = 1.;

    m(af::span, af::span, 17, af::span) = 0;
    m(af::span, af::span, 18, af::span) = 0;

    m(af::span, af::span, 20, af::span) = 0;
    m(af::span, af::span, 21, af::span) = 0;

    m(af::span, af::span, 23, af::span) = 0;
    m(af::span, af::span, 24, af::span) = 0;

    m(af::span, af::span, 26, af::span) = 0;
    m(af::span, af::span, 27, af::span) = 0;

    m(af::span, af::span, 29, af::span) = 0;
    m(af::span, af::span, 30, af::span) = 0;

    // defining interactions
    auto demag = uptr_FieldTerm(new DemagField(mesh, true, true, 0));
    auto exch = uptr_FieldTerm(new RKKYExchangeField(RKKY_values(RKKY), Exchange_values(A), mesh));
    auto aniso = uptr_FieldTerm(new UniaxialAnisotropyField(Ku, (std::array<double, 3>){0, 0, 1}));

    // NOTE try//auto dmi = uptr_FieldTerm (new DmiField(SK_D, {0, 0, -1}));
    auto dmi = uptr_FieldTerm(new DmiField(D, {0, 0, 1}));

    array zee = constant(0.0, mesh.nx, mesh.ny, mesh.nz, 3, f64);
    zee(af::span, af::span, af::span, 2) = Hz;
    auto external = uptr_FieldTerm(new ExternalField(zee));

    // af::print("dmi", dmi->h(state_1));
    // af::print("exch", exch->h(state_1));

    LLGIntegrator llg(1, {demag, exch, aniso, dmi, external});
    timer.print_stage("init ");

    State state_1(mesh, Ms, m);

    if (!exists(filepath + "m_relaxed.vti")) {
        std::cout << "Relaxing minit" << std::endl;
        state_1.write_vti(filepath + "minit");
        // LLGIntegrator llg(1, {demag, exch, aniso, dmi, external});
        while (state_1.t < 3e-9) {
            if (state_1.steps % 100 == 0)
                state_1.write_vti(filepath + "m_step" + std::to_string(state_1.steps));
            llg.step(state_1);
            std::cout << std::scientific << state_1.steps << "\t" << state_1.t << "\t" << state_1.meani(2) << "\t"
                      << llg.E(state_1) << std::endl;
        }
        // llg.relax(state_1);
        timer.print_stage("relax");
        state_1.write_vti(filepath + "m_relaxed");
        vti_writer_micro(state_1.m(af::span, af::span, 0, af::span), Mesh(nx, ny, 1, dx, dy, dz),
                         filepath + "m_relaxed_bottom_ferro_layer1");
        vti_writer_micro(state_1.m(af::span, af::span, 12, af::span), Mesh(nx, ny, 1, dx, dy, dz),
                         filepath + "m_relaxed_bottom_ferro_layer5");
        vti_writer_micro(state_1.m(af::span, af::span, 13, af::span), Mesh(nx, ny, 1, dx, dy, dz),
                         filepath + "m_relaxed__ferri_layer1");
        vti_writer_micro(state_1.m(af::span, af::span, 19, af::span), Mesh(nx, ny, 1, dx, dy, dz),
                         filepath + "m_relaxed_top_ferro_layer1");
        vti_writer_micro(state_1.m(af::span, af::span, 31, af::span), Mesh(nx, ny, 1, dx, dy, dz),
                         filepath + "m_relaxed_top_ferro_layer5");
    } else {
        std::cout << "Found m_relaxed, reading in state_1." << std::endl;
        state_1._vti_reader(filepath + "m_relaxed.vti");
        state_1.write_vti(filepath + "m_relaxed_from_read_in");
    }

    // preparing string method
    double n_interp = 60;
    double string_dt = 1e-13;
    // const int string_steps = 1000;

    State state_2(mesh, Ms, state_1.m);
    state_2.m(af::span, af::span, 19, af::span) = 0;
    state_2.m(af::span, af::span, 22, af::span) = 0;
    state_2.m(af::span, af::span, 25, af::span) = 0;
    state_2.m(af::span, af::span, 28, af::span) = 0;
    state_2.m(af::span, af::span, 31, af::span) = 0;

    state_2.m(af::span, af::span, 19, 2) = 1;
    state_2.m(af::span, af::span, 22, 2) = 1;
    state_2.m(af::span, af::span, 25, 2) = 1;
    state_2.m(af::span, af::span, 28, 2) = 1;
    state_2.m(af::span, af::span, 31, 2) = 1;

    vti_writer_micro(state_2.m, state_2.mesh, filepath + "m_state_2_init");

    while (state_2.t < 3e-9) {
        if (state_2.steps % 100 == 0)
            state_2.write_vti(filepath + "m_state2_relaxing" + std::to_string(state_2.steps));
        llg.step(state_2);
        std::cout << std::scientific << state_2.steps << "\t" << state_2.t << "\t" << state_2.meani(2) << "\t"
                  << llg.E(state_2) << std::endl;
    }
    vti_writer_micro(state_2.m, state_2.mesh, filepath + "m_state_2_relaxed");
    timer.print_stage("state2 relaxed");

    std::vector<State> inputimages;
    inputimages.push_back(state_1);
    inputimages.push_back(state_2);

    StringMethod string(state_1, inputimages, n_interp, string_dt, llg);
    string.run(filepath);
    timer.print_stage("string relaxed");
    return 0;
}
