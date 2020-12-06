#include "arrayfire.h"
#include "magnum_af.hpp"
#include <iostream>

using namespace magnumafcpp;

using namespace af;

// Layer numeration
// 0 BL1 (Bottom Layer, FeCo)
// 1 vacuum
// 2 vacuum
// 3 BL2
// 4 vacuum
// 5 vacuum
// 6 BL3
// 7 vacuum
// 8 vacuum
// 9 BL4
// 10 vacuum
// 11 vacuum
// 12 BL5
// 13 vacuum Pt intermediate layer
// 14 IL (Ferrimagnetic layer)
// 15 IL
// 16 IL
// 17 IL
// 18 vacuum
// 19 vacuum
// 20 TL1 (Top layer, FeCo)
// 21 vacuum
// 22 vacuum
// 23 TL2
// 24 vacuum
// 25 vacuum
// 26 TL3
// 27 vacuum
// 28 vacuum
// 29 TL4
// 30 vacuum
// 31 vacuum
// 32 TL5
// 33 vacuum
// 34 vacuum
// 35 vacuum
// 36 vacuum
// 37 vacuum
// 38 vacuum
// 39 vacuum
// 40 vacuum
// 41 vacuum
// 42 vacuum
// 43 vacuum
// 44 vacuum
// 45 vacuum
// 46 vacuum
// 47 vacuum
// 48 vacuum
// 49 vacuum
// 50 vacuum
// 51 vacuum
// 52 vacuum
// 53 vacuum
void set_air_to_zero(af::array& m) {
    m(af::span, af::span, 1, af::span) = 0;
    m(af::span, af::span, 2, af::span) = 0;

    m(af::span, af::span, 4, af::span) = 0;
    m(af::span, af::span, 5, af::span) = 0;

    m(af::span, af::span, 7, af::span) = 0;
    m(af::span, af::span, 8, af::span) = 0;

    m(af::span, af::span, 10, af::span) = 0;
    m(af::span, af::span, 11, af::span) = 0;

    m(af::span, af::span, 18, af::span) = 0;
    m(af::span, af::span, 19, af::span) = 0;

    m(af::span, af::span, 21, af::span) = 0;
    m(af::span, af::span, 22, af::span) = 0;

    m(af::span, af::span, 24, af::span) = 0;
    m(af::span, af::span, 25, af::span) = 0;

    m(af::span, af::span, 27, af::span) = 0;
    m(af::span, af::span, 28, af::span) = 0;

    m(af::span, af::span, 30, af::span) = 0;
    m(af::span, af::span, 31, af::span) = 0;

    m(af::span, af::span, af::seq(33, af::end), af::span) = 0;
}

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

    // icase 1: DMI in IL = 0, minit skyrm in all layers
    // icase 2: DMI in IL = 0.8e-3, minit skyrm in all layers
    // icase 2: DMI in IL = 0.8e-3, minit skyrm in all layers but IL
    const int icase(argc > 3 ? std::stoi(argv[3]) : 1);
    std::cout << "case=" << icase << std::endl;
    const double r_ratio_to_sample(argc > 4 ? std::stod(argv[4]) : 16.); // r=4 means radius of init skryrm is 1/4 x
    std::cout << "r_ratio_to_sample=" << r_ratio_to_sample << std::endl;
    const double t_relax_state1_sec(argc > 5 ? std::stod(argv[5]) : 3e-9); // relaxtime for state1 in seconds
    std::cout << "t_relax_state1_sec=" << t_relax_state1_sec << std::endl;
    const bool int_over_relax = true; // true== interate, false relax
    // Parameter initialization
    const double x = 400e-9;
    const double y = 400e-9;
    const double z = 32e-9;

    const unsigned ixy(argc > 6 ? std::stoi(argv[6]) : 100);
    const unsigned nx = ixy, ny = ixy, nz = 32;
    const unsigned layers_above_top = 22; // adds 22 nm vacuum above top layer to evaluate demag. We evaluate
                                          // at 21 nm, but need 22 layers to prevent prime number of layers
                                          // conflicting with opencl fft.
    std::cout << "nx = " << ixy << std::endl;
    const double dx = x / nx;
    const double dy = y / ny;
    const double dz = z / nz;
    std::cout << "nx = " << ixy << " " << dx << " " << dy << " " << dz << std::endl;

    const double Hz = 130e-3 / constants::mu0;
    const double RKKY_val = 0.8e-3 * 1e-9; // TODO
    // Note: maybe 0.5 factor (see mumax3):
    // const double RKKY_val = 0.8e-3 * 1e-9* 0.5;//

    // SK layer params
    const double SK_Ms = 1371e3;  // A/m
    const double SK_A = 15e-12;   // J/m
    const double SK_Ku = 1.411e6; // J/m^3
    const double SK_D = -2.5e-3;  // J/m^2

    // Ferrimagnetic interface layer params
    const double IL_Ms = 488.2e3;                   // A/m
    const double IL_A = 4e-12;                      // J/m
    const double IL_Ku = 486.6e3;                   // J/m^3
    const double IL_D = (icase == 1 ? 0. : 0.8e-3); // J/m^2

    // Pt layer is vacuum like, strongly exchange coupled to layer 12
    // Layer 13 (i.e. 1nm Pt interlayer) is strongly exchange coupled to layer
    // 12 with very low Js. js 1e-4 T, A = 10 A_rest, K=0; gros/ RKKY
    //    const double Pt_Ms = 1e-3;  // A/m
    //    const double Pt_A = 10 * 15e-12;   // J/m
    //    const double Pt_Ku = 1.411e6; // J/m^3
    //    const double Pt_D = 1e-3 * -2.5e-3;   // J/m^2
    //    //const double Pt_D = 0;   // J/m^2
    const double Pt_Ms = 5e3;       // A/m //NOTE: Small Ms highly decreases integration speed!!! Ms
                                    // 1e3 or below causes smaller than 1e-15 s steps.
    const double Pt_A = 2 * 15e-12; // J/m
    const double Pt_Ku = 0;         // J/m^3
    // const double Pt_Ku = 1.411e6; // J/m^3
    // const double Pt_D = 1e-3 * -2.5e-3;   // J/m^2
    const double Pt_D = 0; // J/m^2

    // const double Pt_Ms = 1371e0;  // A/m
    // const double Pt_A = 15e-12;   // J/m
    // const double Pt_Ku = 1.411e6; // J/m^3
    // const double Pt_D = -2.5e-3;   // J/m^2

    array geom = af::constant(0.0, nx, ny, nz + layers_above_top, 3, f64);
    geom(af::span, af::span, 0, af::span) = 1;
    geom(af::span, af::span, 3, af::span) = 1;
    geom(af::span, af::span, 6, af::span) = 1;
    geom(af::span, af::span, 9, af::span) = 1;
    geom(af::span, af::span, 12, af::span) = 1;

    geom(af::span, af::span, 13, af::span) = 3;

    geom(af::span, af::span, 14, af::span) = 2;
    geom(af::span, af::span, 15, af::span) = 2;
    geom(af::span, af::span, 16, af::span) = 2;
    geom(af::span, af::span, 17, af::span) = 2;

    geom(af::span, af::span, 20, af::span) = 1;
    geom(af::span, af::span, 23, af::span) = 1;
    geom(af::span, af::span, 26, af::span) = 1;
    geom(af::span, af::span, 29, af::span) = 1;
    geom(af::span, af::span, 32, af::span) = 1;

    // af::print("geom", geom);
    // af::print("geom == 1", geom == 1);
    // af::print("geom == 2", geom == 2);

    af::array todouble = af::constant(1., nx, ny, nz + layers_above_top, 3, f64);
    af::array Ms = (SK_Ms * (geom == 1) + IL_Ms * (geom == 2) + Pt_Ms * (geom == 3)) * todouble;
    af::array A = (SK_A * (geom == 1) + IL_A * (geom == 2) + Pt_A * (geom == 3)) * todouble;
    af::array Ku = (SK_Ku * (geom == 1) + IL_Ku * (geom == 2) + Pt_Ku * (geom == 3)) * todouble;
    af::array D = (SK_D * (geom == 1) + IL_D * (geom == 2) + Pt_D * (geom == 3)) * todouble; // + IL_D * (geom == 2);
    // std::cout << "type" << geom.type() << " " << Ms.type() << std::endl;

    // af::print("Ms", Ms);
    // af::print("1e12 * A", 1e12 * A);
    // af::print("Ku", Ku);
    // af::print("D", D);

    array RKKY = af::constant(0.0, nx, ny, nz + layers_above_top, 3, f64);
    // RKKY(af::span, af::span, 12, af::span) = RKKY_val;
    RKKY(af::span, af::span, 13, af::span) = RKKY_val;
    RKKY(af::span, af::span, 14, af::span) = RKKY_val;

    // Generating Objects
    Mesh mesh(nx, ny, nz + layers_above_top, dx, dy, dz);

    vti_writer_micro(Ms, mesh, filepath + "Ms_field");
    vti_writer_micro(A, mesh, filepath + "A_field");
    vti_writer_micro(Ku, mesh, filepath + "Ku_field");
    vti_writer_micro(D, mesh, filepath + "D_field");
    vti_writer_micro(RKKY, mesh, filepath + "RKKY_field");

    // Initial magnetic field
    array m = constant(0.0, mesh.nx, mesh.ny, mesh.nz, 3, f64);
    m(af::span, af::span, af::span, 2) = -1;
    for (unsigned ix = 0; ix < mesh.nx; ix++) {
        for (unsigned iy = 0; iy < mesh.ny; iy++) {
            const double rx = double(ix) - mesh.nx / 2.;
            const double ry = double(iy) - mesh.ny / 2.;
            const double r = sqrt(pow(rx, 2) + pow(ry, 2));
            if (r > nx / r_ratio_to_sample)
                m(ix, iy, af::span, 2) = 1.;
        }
    }

    if (icase == 3) {
        // Setting IL layer ferromagnetic
        m(af::span, af::span, 14, 2) = 1.;
        m(af::span, af::span, 15, 2) = 1.;
        m(af::span, af::span, 16, 2) = 1.;
        m(af::span, af::span, 17, 2) = 1.;
    }
    set_air_to_zero(m);

    // defining interactions
    auto demag = uptr_Fieldterm(new DemagField(mesh, true, true, 0));
    auto exch = uptr_Fieldterm(new RKKYExchangeField(RKKY_values(RKKY), Exchange_values(A), mesh));
    auto aniso = uptr_Fieldterm(new UniaxialAnisotropyField(Ku, (std::array<double, 3>){0, 0, 1}));

    // NOTE try//auto dmi = uptr_Fieldterm (new DmiField(SK_D, {0, 0, -1}));
    auto dmi = uptr_Fieldterm(new DmiField(D, {0, 0, -1})); // TODO current definition, will change sign with update

    array zee = constant(0.0, mesh.nx, mesh.ny, mesh.nz, 3, f64);
    zee(af::span, af::span, af::span, 2) = Hz;
    auto external = uptr_Fieldterm(new ExternalField(zee));

    LLGIntegrator llg(1, {demag, exch, aniso, dmi, external});
    timer.print_stage("init ");

    State state_1(mesh, Ms, m);

    af::print("m", af::mean(af::mean(af::mean(state_1.m, 0), 1), 2));
    af::print("Ms", af::mean(af::mean(af::mean(state_1.Ms_field, 0), 1), 2));
    af::print("demag", af::mean(af::mean(af::mean(demag->h(state_1), 0), 1), 2));
    af::print("exch", af::mean(af::mean(af::mean(exch->h(state_1), 0), 1), 2));
    af::print("aniso", af::mean(af::mean(af::mean(aniso->h(state_1), 0), 1), 2));
    af::print("dmi", af::mean(af::mean(af::mean(dmi->h(state_1), 0), 1), 2));
    // af::print("dmi", dmi->h(state_1));
    // af::print("exch", exch->h(state_1));

    std::cout << "Relaxing minit" << std::endl;
    state_1.write_vti(filepath + "minit");
    std::cout << "Prestart:" << std::scientific << state_1.steps << "\t" << state_1.t << "\t" << state_1.meani(2)
              << "\t" << llg.E(state_1) << std::endl;
    // auto demag_init = demag->h(state_1) * 1e3 * constants::mu0; // in mT
    // vti_writer_micro(demag_init, mesh, filepath + "demag_init");
    // vti_writer_micro(demag_init(af::span, af::span, 53, af::span), mesh,
    // filepath + "demag_init_21nm_above_top");

    // write_ascii(demag_init, mesh, filepath + "demag_init.txt");
    // write_ascii(demag_init(af::span, af::span, 53, af::span), mesh, filepath
    // + "demag_init_21nm_above_top.txt");

    while (state_1.t < t_relax_state1_sec) {
        if (state_1.steps % 100 == 0) {
            // state_1.write_vti(filepath + "m_step" +
            // std::to_string(state_1.steps));
            state_1.write_vti(filepath + "m_current_step");
            af::array m_current_cleanPtLayer = state_1.m;
            m_current_cleanPtLayer(af::span, af::span, 13, af::span) = 0;
            vti_writer_micro(m_current_cleanPtLayer, mesh, filepath + "m_current_step_cleanPtlayer");
            vti_writer_micro(m_current_cleanPtLayer(af::span, af::span, af::seq(nz + 1), af::span),
                             Mesh(nx, ny, nz + 1, dx, dy, dz), filepath + "m_current_step_cleanPtlayer_no_airbox");

            auto demag_current = demag->h(state_1) * 1e3 * constants::mu0; // in mT
            vti_writer_micro(demag_current(af::span, af::span, 53, af::span), mesh,
                             filepath + "demag_current_21nm_above_top");
            // vti_writer_micro(demag_current, mesh, filepath +
            // "demag_current");

            // write_ascii(demag_current, mesh, filepath + "demag_current.txt");
            // write_ascii(demag_current(af::span, af::span, 53, af::span),
            // mesh, filepath + "demag_current_21nm_above_top.txt");
        }
        llg.step(state_1);
        std::cout << std::scientific << state_1.steps << "\t" << state_1.t << "\t" << state_1.meani(2) << "\t"
                  << llg.E(state_1) << std::endl;
    }

    timer.print_stage("relax");
    auto demag_relaxed = demag->h(state_1) * 1e3 * constants::mu0; // in mT
    vti_writer_micro(demag_relaxed, mesh, filepath + "demag_relaxed");
    vti_writer_micro(demag_relaxed(af::span, af::span, 53, af::span), mesh, filepath + "demag_relaxed_21nm_above_top");

    write_ascii(demag_relaxed, mesh, filepath + "demag_relaxed.txt");
    write_ascii(demag_relaxed(af::span, af::span, 53, af::span), mesh, filepath + "demag_relaxed_21nm_above_top.txt");

    state_1.write_vti(filepath + "m_relaxed");

    af::array m_relaxed_cleanPtLayer = state_1.m;
    m_relaxed_cleanPtLayer(af::span, af::span, 13, af::span) = 0;
    vti_writer_micro(m_relaxed_cleanPtLayer, mesh,
                     filepath + "m_relaxed_cleanPtlayer"); // For plotting
    vti_writer_micro(m_relaxed_cleanPtLayer(af::span, af::span, af::seq(nz + 1), af::span),
                     Mesh(nx, ny, nz + 1, dx, dy, dz), filepath + "m_relaxed_step_cleanPtlayer_no_airbox");

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

    return 0;
}
