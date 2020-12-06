#include "arrayfire.h"
#include "llg.hpp"
#include "micro/external_field.hpp"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
//#include "micro/demag_field.hpp"
//#include "micro/exchange_field.hpp"
#include "atom/atomistic_uniaxial_anisotropy_field.hpp"
#include "atom/atomistic_dipole_dipole_field.hpp"
#include "atom/atomistic_dmi_field.hpp"
#include "atom/atomistic_exchange_field.hpp"
#include "solvers/string_method.hpp"
#include "vtk_writer.hpp"

using namespace magnumafcpp;

using namespace af;

int main(int argc, char** argv) {
    if (argc > 1)
        setDevice(std::stoi(argv[2]));
    info();

    std::cout << "argc" << argc << std::endl;
    for (int i = 0; i < argc; i++)
        cout << "Parameter " << i << " was " << argv[i] << "\n";
    // char* charptr;
    // std::cout<<"argv"<<std::strtod(argv[1], &charptr)<<std::endl;
    // std::cout<<"argv"<<std::strtod(argv[2], &charptr)<<std::endl;

    // Parameter initialization
    const int nx = 112, ny = 112, nz = 1; // nz=5 -> lz=(5-1)*dx
    const double dx = 2.715e-10;

    // Simulation Parameters
    double hmax = 3.5e-10;
    double hmin = 1.0e-15;

    double atol = 1e-6;

    double rtol = atol;

    double n_interp = 60;
    double string_dt = 5e-13;
    const int string_steps = 10000;
    std::string filepath(argc > 0 ? argv[1] : "../Data/Testing/");
    if (argc > 0)
        filepath.append("/");
    // else std::string filepath("../Data/Testing/");
    std::cout << "Writing into path " << filepath.c_str() << std::endl;

    // Generating Objects
    Mesh mesh(nx, ny, nz, dx, dx, dx);
    Material material = Material();
    state.Ms = 1.1e6;
    material.A = 1.5e-11; // TODO why this value? Check!
    material.alpha = 1;
    state.material.afsync = false;

    // material.D=2*5.76e-3;
    // material.D=3.0e-3;
    material.D = std::stod(argv[3]);
    std::cout << "D=" << material.D << std::endl;
    // material.D_axis[2]=-1;

    // material.Ku1=510e3;
    // material.Ku1=800399;
    material.Ku1 = std::stod(argv[4]);
    std::cout << "Ku1=" << material.Ku1 << std::endl;
    // material.Ku1_axis[2]=1;

    material.J_atom = 2. * material.A * dx;
    material.D_atom = material.D * pow(dx, 2);
    std::cout << "D_atom=" << material.D_atom << std::endl;
    material.K_atom = material.Ku1 * pow(dx, 3);
    std::cout << "Ku1_atom=" << material.K_atom << std::endl;
    material.p = state.Ms * pow(dx, 3); // Compensate nz=1 instead of nz=4

    // Initial magnetic field
    array m = constant(0.0, mesh.nx, mesh.ny, mesh.nz, 3, f64);
    m(span, span, span, 2) = -1;
    for (int ix = 0; ix < mesh.nx; ix++) {
        for (int iy = 0; iy < mesh.ny; iy++) {
            const double rx = double(ix) - mesh.nx / 2.;
            const double ry = double(iy) - mesh.ny / 2.;
            const double r = sqrt(pow(rx, 2) + pow(ry, 2));
            if (r > ny / 4.)
                m(ix, iy, span, 2) = 1.;
        }
    }

    State state(mesh, material, m);
    af_to_vti(state.m, mesh, (filepath + "minit").c_str());

    std::vector<uptr_FieldTerm> llgterm;
    // llgterm.push_back( uptr_FieldTerm (new DemagField(mesh, material)));
    // llgterm.push_back( uptr_FieldTerm (new ExchangeField(mesh, material)));
    // llgterm.push_back( uptr_FieldTerm (new DmiField(mesh, material)));
    // llgterm.push_back( uptr_FieldTerm (new UniaxialAnisotropyField(mesh,
    // material)));

    llgterm.push_back(uptr_FieldTerm(new AtomisticDipoleDipoleField(mesh)));
    llgterm.push_back(uptr_FieldTerm(new AtomisticExchangeField(mesh)));
    llgterm.push_back(uptr_FieldTerm(new AtomisticDmiField(mesh, material)));
    llgterm.push_back(uptr_FieldTerm(new AtomisticUniaxialAnisotropyField(mesh, material)));

    LLG Llg(state, atol, rtol, hmax, hmin, llgterm);

    double energy_n0 = 1.e20;
    double energy_n1 = 1.e30;
    int irel = 0;
    timer t = af::timer::start();
    while (state.t < 1.e-8) {
        state.m = Llg.step(state);
        irel++;
        if (irel % 100 == 0)
            std::cout << "irel: " << irel << " at time: " << state.t << " Energy: " << Llg.E(state) << std::endl;
        energy_n1 = Llg.E(state);
        if (fabs(2 * (energy_n1 - energy_n0) / (energy_n1 + energy_n0)) < 1e-10) {
            std::cout << "Exiting loop: initial image relaxated,  relative "
                         "difference smaller than 1e-10"
                      << std::endl;
            break;
        }
        energy_n0 = Llg.E(state);
    }
    double timerelax = af::timer::stop(t);
    af_to_vti(state.m, mesh, (filepath + "relax").c_str());

    std::cout << "Relaxed after:" << state.t << " timerelax [af-s]: " << timerelax << " for "
              << Llg.counter_accepted + Llg.counter_reject << " steps, thereof " << Llg.counter_accepted
              << " Steps accepted, " << Llg.counter_reject << " Steps rejected" << std::endl;

    array last = constant(0, dims_vector(mesh), f64);
    last(span, span, span, 2) = 1;

    std::vector<State> inputimages;
    inputimages.push_back(state);
    inputimages.push_back(State(mesh, material, last));

    StringMethod string(state, inputimages, n_interp, string_dt, llgterm);
    // StringMethod* string = new StringMethod(state, inputimages, n_interp , llgterm);
    std::cout.precision(12);

    std::ofstream stream_E_barrier;
    stream_E_barrier.precision(12);
    stream_E_barrier.open((filepath + "E_barrier.dat").c_str());

    std::ofstream stream_steps;
    stream_steps.precision(12);
    stream_steps.open((filepath + "steps.dat").c_str());

    std::ofstream stream_E_curves;
    stream_E_curves.precision(12);
    stream_E_curves.open((filepath + "E_curves.dat").c_str());

    double max_lowest = 1e100;
    double max_prev_step = 1e100;
    int i_max_lowest = -1;
    std::vector<State> images_max_lowest;
    std::vector<double> E_max_lowest;
    for (int i = 0; i < string_steps; i++) {
        string.step();
        string.calc_E();
        auto max = std::max_element(std::begin(string.E), std::end(string.E));
        if (*max - string.E[0] < max_lowest) {
            max_lowest = *max - string.E[0];
            i_max_lowest = i;
            images_max_lowest = string.images;
            E_max_lowest = string.E;
        }
        // Wrong approach
        // else if(i>50){
        //  std::cout   << "Exiting loop: Energy barrier after 50step relaxation
        //  becomes bigger "<<std::endl; stream_steps<<"#Exiting loop: Energy
        //  barrier after 50step relaxation becomes bigger "<<std::endl; break;
        //}

        // std::cout<<"Test: fabs=
        // "<<fabs(2*(*max-string.E[0]-max_prev_step)/(*max-string.E[0]+max_prev_step))<<std::endl;

        if (i > 25 && fabs(2 * (*max - string.E[0] - max_prev_step) / (*max - string.E[0] + max_prev_step)) < 1e-6) {
            std::cout << "Exiting loop: Energy barrier relative difference "
                         "smaller than 1e-6"
                      << std::endl;
            stream_steps << "#Exiting loop: Energy barrier relative difference "
                            "smaller than 1e-6"
                         << std::endl;
            break;
        }
        if (i > 25 && fabs(*max - string.E[0] - max_prev_step) < 1e-27) {
            std::cout << "Exiting loop: Energy barrier difference smaller than 1e-27" << std::endl;
            stream_steps << "#Exiting loop: Energy barrier difference smaller than 1e-27" << std::endl;
            break;
        }
        std::cout << i << "\t" << *max - string.E[0] << "\t" << string.E[0] << "\t" << *max - string.E[-1] << "\t"
                  << *max << "\t"
                  << fabs(2 * (*max - string.E[0] - max_prev_step) / (*max - string.E[0] + max_prev_step)) << std::endl;
        stream_steps << i << "\t" << *max - string.E[0] << "\t" << string.E[0] << "\t" << *max - string.E[-1] << "\t"
                     << *max << "\t"
                     << fabs(2 * (*max - string.E[0] - max_prev_step) / (*max - string.E[0] + max_prev_step))
                     << std::endl;
        for (unsigned j = 0; j < string.E.size(); ++j) {
            stream_E_curves << i << " " << j << " " << string.E[j] - string.E[0] << " " << string.E[j] - string.E[-1]
                            << " " << string.E[j] << std::endl;
        }
        stream_E_curves << i << "\n\n" << std::endl;
        max_prev_step = *max - string.E[0];
        if (i % 20 == 1) {
            std::cout << "Writing current skyrm images for iteration" << i << std::endl;
            for (unsigned j = 0; j < string.images.size(); j++) {
                std::string name = filepath;
                name.append("current_skyrm_image");
                name.append(std::to_string(j));
                af_to_vti(string.images[j].m, mesh, name.c_str());
            }
        }
    }
    std::cout << "#i , lowest overall:   max-[0], max-[-1] max [J]: " << i_max_lowest << "\t" << max_lowest << "\t"
              << max_lowest + E_max_lowest[0] - E_max_lowest[-1] << "\t" << max_lowest + E_max_lowest[0] << std::endl;
    stream_steps << "#i , lowest overall:   max-[0], max-[-1] max [J]: " << i_max_lowest << "\t" << max_lowest << "\t"
                 << max_lowest + E_max_lowest[0] - E_max_lowest[-1] << "\t" << max_lowest + E_max_lowest[0]
                 << std::endl;
    stream_E_barrier << max_lowest << "\t" << nx << "\t" << dx << "\t" << material.D << "\t" << material.Ku1 << "\t"
                     << material.D_atom << "\t" << material.K_atom << "\t" << std::endl;

    std::ofstream myfileE;
    myfileE.precision(12);
    myfileE.open((filepath + "E_last_step.dat").c_str());

    std::ofstream stream_max_lowest;
    stream_max_lowest.precision(12);
    stream_max_lowest.open((filepath + "E_max_lowest.dat").c_str());

    std::cout << string.E.size() << "\t" << string.images.size() << "\t" << std::endl;
    for (unsigned i = 0; i < string.images.size(); i++) {
        std::cout << "i=" << i << "\t"
                  << "E= " << string.E[i] << std::endl;
        myfileE << i << "\t" << string.E[i] << "\t" << string.E[i] - string.E[0] << "\t" << string.E[i] - string.E[-1]
                << std::endl;
        std::string name = filepath;
        name.append("skyrm_image");
        name.append(std::to_string(i));
        af_to_vti(string.images[i].m, mesh, name.c_str());
        // af_to_vtk(string.images_interp[i], name.c_str());
        stream_max_lowest << i << "\t" << E_max_lowest[i] << "\t" << E_max_lowest[i] - E_max_lowest[0] << "\t"
                          << E_max_lowest[i] - E_max_lowest[-1] << std::endl;
        name = filepath;
        name.append("skyrm_image_max_lowest");
        name.append(std::to_string(i));
        af_to_vti(images_max_lowest[i].m, mesh, name.c_str());
    }

    for (unsigned i = 0; i < Llg.Fieldterms.size(); ++i) {
        std::cout << "get_cpu_time()" << std::endl;
        std::cout << i << "\t" << *Llg.Fieldterms[i]->get_cpu_time() << std::endl;
        stream_steps << "#"
                     << "get_cpu_time()" << std::endl;
        stream_steps << "#" << i << "\t" << *Llg.Fieldterms[i]->get_cpu_time() << std::endl;
    }

    myfileE.close();
    stream_E_barrier.close();
    stream_steps.close();
    stream_E_curves.close();
    stream_max_lowest.close();
    // delete[] string;

    return 0;
}
