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
}

void write_vti_and_cleanlayer(const State& state, std::string fullfilepath, unsigned cleanlayer = 13) {
    state.write_vti(fullfilepath);
    af::array m_current_cleanPtLayer = state.m;
    m_current_cleanPtLayer(af::span, af::span, cleanlayer, af::span) = 0;
    vti_writer_micro(m_current_cleanPtLayer, state.mesh, fullfilepath + "cleanPtlayer");
}

// write_vti_and_cleanlayer(state, filepath + );

int main(int argc, char** argv) {

    std::cout << "argc = " << argc << std::endl;
    for (int i = 0; i < argc; i++)
        std::cout << "Parameter " << i << " was " << argv[i] << "\n";

    std::string filepath(argc > 1 ? argv[1] : "./run/");
    if (argc > 1)
        filepath.append("/");
    std::cout << "Writing into path " << filepath.c_str() << std::endl;

    af::setDevice(argc > 2 ? std::stoi(argv[2]) : 0);
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

    const unsigned ixy(argc > 6 ? std::stoi(argv[6]) : 100);
    const double Factor_IL_MS_over_Pt_Ms(argc > 7 ? std::stod(argv[7]) : 100.);
    std::cout << "Factor_IL_MS_over_Pt_Ms = " << Factor_IL_MS_over_Pt_Ms << std::endl;
    const unsigned nx = ixy, ny = ixy, nz = 33;
    std::cout << "nx = " << ixy << std::endl;
    const double dx = x / nx;
    const double dy = y / ny;
    const double dz = 1e-9;
    std::cout << "nx = " << ixy << " " << dx << " " << dy << " " << dz << std::endl;

    // const double Hz = 130e-3 / constants::mu0;
    const double RKKY_val = 0.8e-3 * dz; // TODO
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
    const double Pt_Ms = IL_Ms / Factor_IL_MS_over_Pt_Ms; // A/m //NOTE: Small Ms highly decreases integration speed!!!
                                                          // Ms 1e3 or below causes smaller than 1e-15 s steps.
    // const double Pt_Ms = 5e3;       // A/m //NOTE: Small Ms highly decreases integration speed!!! Ms
    //                                // 1e3 or below causes smaller than 1e-15 s steps.
    const double Pt_A = 2 * 15e-12; // J/m
    const double Pt_Ku = 0;         // J/m^3
    const double Pt_D = 0;          // J/m^2

    // 1: Bottom and Top SK layer
    // 2: IL ferrimagnetic interlayer
    // 3: Pt intermediate layer
    array geom = af::constant(0.0, nx, ny, nz, 3, f64);
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

    af::array todouble = af::constant(1., nx, ny, nz, 3, f64);
    af::array Ms = (SK_Ms * (geom == 1) + IL_Ms * (geom == 2) + Pt_Ms * (geom == 3)) * todouble;
    af::array A = (SK_A * (geom == 1) + IL_A * (geom == 2) + Pt_A * (geom == 3)) * todouble;
    af::array Ku = (SK_Ku * (geom == 1) + IL_Ku * (geom == 2) + Pt_Ku * (geom == 3)) * todouble;
    af::array D = (SK_D * (geom == 1) + IL_D * (geom == 2) + Pt_D * (geom == 3)) * todouble; // + IL_D * (geom == 2);
    // std::cout << "type" << geom.type() << " " << Ms.type() << std::endl;

    // af::print("Ms", Ms);
    // af::print("1e12 * A", 1e12 * A);
    // af::print("Ku", Ku);
    // af::print("D", D);

    array RKKY = af::constant(0.0, nx, ny, nz, 3, f64);
    RKKY(af::span, af::span, 13, af::span) = RKKY_val;
    RKKY(af::span, af::span, 14, af::span) = RKKY_val;

    // Generating Objects
    Mesh mesh(nx, ny, nz, dx, dy, dz);

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
    set_air_to_zero(m);

    if (icase == 3) {
        // Setting IL layer ferromagnetic
        m(af::span, af::span, 14, 2) = 1.;
        m(af::span, af::span, 15, 2) = 1.;
        m(af::span, af::span, 16, 2) = 1.;
        m(af::span, af::span, 17, 2) = 1.;
    }

    // defining interactions
    auto demag = LlgTerm(new DemagField(mesh, true, true, 0));
    auto exch = LlgTerm(new RKKYExchangeField(RKKY_values(RKKY), Exchange_values(A), mesh));
    auto aniso = LlgTerm(new UniaxialAnisotropyField(Ku, std::array<double, 3>{0, 0, 1}));

    // NOTE try//auto dmi = LlgTerm (new DmiField(SK_D, {0, 0, -1}));
    auto dmi = LlgTerm(new DmiField(D, {0, 0, -1})); // TODO current definition, will change sign with update

    // af::print("dmi", dmi->h(state_1));
    // af::print("exch", exch->h(state_1));

    LLGIntegrator llg(1, {demag, exch, aniso, dmi});
    State state_1(mesh, Ms, m);

    timer.print_stage("init ");

    if (!exists(filepath + "m_relaxed.vti")) {
        std::cout << "Relaxing minit" << std::endl;
        double Hz_init = 130e-3 / constants::mu0;
        af::array zee_init = constant(0.0, mesh.nx, mesh.ny, mesh.nz, 3, f64);
        zee_init(af::span, af::span, af::span, 2) = Hz_init;
        llg.llgterms.push_back(LlgTerm(new ExternalField(zee_init)));

        state_1.write_vti(filepath + "minit");
        if (int_over_relax) {
            while (state_1.t < t_relax_state1_sec) {
                llg.step(state_1);
                if (state_1.steps % 1000 == 0) {
                    write_vti_and_cleanlayer(state_1, filepath + "m_step" + std::to_string(state_1.steps));
                    // state_1.write_vti(filepath + "m_step" + std::to_string(state_1.steps));
                    // af::array m_current_cleanPtLayer = state_1.m;
                    // m_current_cleanPtLayer(af::span, af::span, 13, af::span) = 0;
                    // vti_writer_micro(m_current_cleanPtLayer, mesh, filepath + "m_current_step_cleanPtlayer");
                }
                std::cout << std::scientific << state_1.steps << "\t" << state_1.t << "\t" << state_1.meani(2) << "\t"
                          << llg.E(state_1) << std::endl;
            }
        } else {
            llg.relax(state_1, 1e-12);
        }
        state_1.t = 0;
        llg.llgterms.pop_back(); // remove zee
        timer.print_stage("relax");
        write_vti_and_cleanlayer(state_1, filepath + "m_relaxed");

        vti_writer_micro(state_1.m(af::span, af::span, 0, af::span), Mesh(nx, ny, 1, dx, dy, dz),
                         filepath + "m_relaxed_bottom_ferro_layer1");
        vti_writer_micro(state_1.m(af::span, af::span, 12, af::span), Mesh(nx, ny, 1, dx, dy, dz),
                         filepath + "m_relaxed_bottom_ferro_layer5");
        vti_writer_micro(state_1.m(af::span, af::span, 13, af::span), Mesh(nx, ny, 1, dx, dy, dz),
                         filepath + "m_relaxed__Pt_layer");
        vti_writer_micro(state_1.m(af::span, af::span, 14, af::span), Mesh(nx, ny, 1, dx, dy, dz),
                         filepath + "m_relaxed__ferri_layer1");
        vti_writer_micro(state_1.m(af::span, af::span, 20, af::span), Mesh(nx, ny, 1, dx, dy, dz),
                         filepath + "m_relaxed_top_ferro_layer1");
        vti_writer_micro(state_1.m(af::span, af::span, 32, af::span), Mesh(nx, ny, 1, dx, dy, dz),
                         filepath + "m_relaxed_top_ferro_layer5");
    } else {
        std::cout << "Found m_relaxed, reading in state_1." << std::endl;
        state_1._vti_reader(filepath + "m_relaxed.vti");
        state_1.write_vti(filepath + "m_relaxed_from_read_in");
    }

    bool relax_per_step_else_continuous = true; // switch
    // bool relax_per_step_else_continuous = false; // switch
    if (relax_per_step_else_continuous) {
        const int iHz_max = 200; // e-3 / constants::mu0;
        // const int iHz_max = 1000; //e-3 / constants::mu0;
        // const int iHz_max = 250; //e-3 / constants::mu0;
        // const double Hz_max = 250e-3 / constants::mu0;
        // const int nsteps = 100;
        std::ofstream stream(filepath + "hys.dat");
        std::cout << "Hz\t<mz> " << std::endl;
        af::array eval_mean_region = af::constant(1, nx, ny, nz, 1, f64);
        set_air_to_zero(eval_mean_region);
        eval_mean_region(af::span, af::span, 13, af::span) = 0;

        // for (int iHz_current = 130; iHz_current < iHz_max; iHz_current++) {
        for (int iHz_current = 130; iHz_current < iHz_max; iHz_current += 2) {
            double Hz_current = iHz_current * 1e-3 / constants::mu0;
            af::array zee = constant(0.0, mesh.nx, mesh.ny, mesh.nz, 3, f64);
            zee(af::span, af::span, af::span, 2) = Hz_current;
            llg.llgterms.push_back(LlgTerm(new ExternalField(zee)));
            llg.relax(state_1, 1e-10);
            // TODO mean must account for empty layers
            auto mean = spacial_mean_in_region(state_1.m, eval_mean_region);
            std::cout << Hz_current * constants::mu0 << "\t" << mean[2] << "\t"
                      << af::mean(af::mean(af::mean(state_1.m, 2), 1), 0)(0, 0, 0, 2).scalar<double>() << "\t"
                      << Hz_current << "\t" << meani(state_1.m, 0) << "\t" << meani(state_1.m, 1) << "\t"
                      << meani(state_1.m, 2) << std::endl;
            stream << Hz_current * constants::mu0 << "\t" << mean[2] << "\t"
                   << af::mean(af::mean(af::mean(state_1.m, 2), 1), 0)(0, 0, 0, 2).scalar<double>() << "\t"
                   << Hz_current << "\t" << meani(state_1.m, 0) << "\t" << meani(state_1.m, 1) << "\t"
                   << meani(state_1.m, 2) << std::endl;
            llg.llgterms.pop_back();
            // state_1.write_vti(filepath + "m_hys" + std::to_string(iHz_current));
            write_vti_and_cleanlayer(state_1, filepath + "m_hys" + std::to_string(iHz_current));
            // llg.llgterms[llg.llgterms.size() - 1] ->
            // external->set_homogeneous_field(0, 0, Hz_current);
            // llg.llgterms[llg.llgterms.size() - 1] ->
        }
        stream.close();
        timer.print_stage("string relaxed");
    } else {
        const double H_strt = 130e-3 / constants::mu0;
        const double H_stop = 200e-3 / constants::mu0;
        const double inttime_in_sec = 200e-9;
        const double H_xz_plane_tilt_in_degree = 1;

        auto calc_Hext = [H_strt, H_stop, inttime_in_sec, H_xz_plane_tilt_in_degree](State state) {
            af::array field = af::constant(0.0, state.mesh.nx, state.mesh.ny, state.mesh.nz, 3, f64);
            const double H_current = (H_stop - H_strt) * state.t / inttime_in_sec + H_strt;
            const double H_x = H_current * std::sin(2 * M_PI * H_xz_plane_tilt_in_degree / 360.);
            const double H_z = H_current * std::cos(2 * M_PI * H_xz_plane_tilt_in_degree / 360.);
            field(af::span, af::span, af::span, 0) = H_x;
            field(af::span, af::span, af::span, 2) = H_z;
            return field;
        };
        llg.llgterms.push_back(LlgTerm(new ExternalField(calc_Hext)));
        while (state_1.t < inttime_in_sec) {
            llg.step(state_1);
            if (state_1.steps % 100 == 0) {
                af::array Hext = af::mean(af::mean(af::mean(llg.llgterms.back()->h(state_1), 0), 1), 2);
                double Hx_in_mT = Hext(0, 0, 0, 0).scalar<double>() * constants::mu0 * 1e3;
                double Hz_in_mT = Hext(0, 0, 0, 2).scalar<double>() * constants::mu0 * 1e3;
                std::cout << state_1 << " " << Hx_in_mT << " " << Hz_in_mT << std::endl;
            }
            if (state_1.steps % 1000 == 0) {
                write_vti_and_cleanlayer(state_1, filepath + "m_hys_step" + std::to_string(state_1.steps));
                write_vti_and_cleanlayer(state_1, filepath + "m_hys_ns" + std::to_string(1e9 * state_1.t));
            }
        }
        timer.print_stage("string relaxed");
    }

    return 0;
}
