#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace magnumafcpp;

using namespace af;
typedef std::shared_ptr<LLGTerm> llgt_ptr;

int main(int argc, char** argv) {

    std::cout << "argc = " << argc << std::endl;
    for (int i = 0; i < argc; i++)
        cout << "Parameter " << i << " was " << argv[i] << "\n";

    std::string filepath(argc > 1 ? argv[1] : "../Data/skyrmion_stoch");
    if (argc > 0)
        filepath.append("/");
    std::cout << "Writing into path " << filepath.c_str() << std::endl;

    setDevice(argc > 2 ? std::stoi(argv[2]) : 0);
    // if(argc>1) setDevice(std::stoi(argv[2]));
    info();

    // Parameter initialization
    const int nx = 30, ny = 30, nz = 1;
    const double dx = 1e-9;

    double n_interp = 60;
    double string_dt = 1e-13;
    const int string_steps = 10000;

    // Generating Objects
    Mesh mesh(nx, ny, nz, dx, dx, dx);
    Material material = Material();
    state.Ms = 1.1e6;
    material.A = 1.6e-11;
    material.alpha = 1;
    material.D = 2 * 5.76e-3;
    material.Ku1 = 6.4e6;

    material.J_atom = 2. * material.A * dx;
    material.D_atom = material.D * pow(dx, 2);
    material.K_atom = material.Ku1 * pow(dx, 3);
    material.p = state.Ms * pow(dx, 3); // Compensate nz=1 instead of nz=4

    // Initial magnetic field
    array m = constant(0.0, mesh.n0, mesh.n1, mesh.n2, 3, f64);
    m(span, span, span, 2) = -1;
    for (int ix = 0; ix < mesh.n0; ix++) {
        for (int iy = 0; iy < mesh.n1; iy++) {
            const double rx = double(ix) - mesh.n0 / 2.;
            const double ry = double(iy) - mesh.n1 / 2.;
            const double r = sqrt(pow(rx, 2) + pow(ry, 2));
            if (r > nx / 4.)
                m(ix, iy, span, 2) = 1.;
        }
    }

    State state(mesh, material, m);
    vti_writer_atom(state.m, mesh, (filepath + "minit").c_str());

    std::vector<llgt_ptr> llgterm;
    // llgterm.push_back( llgt_ptr (new DemagField(mesh, material)));
    // llgterm.push_back( llgt_ptr (new ExchangeField(mesh, material)));
    // llgterm.push_back( llgt_ptr (new DmiField(mesh, material)));
    // llgterm.push_back( llgt_ptr (new UniaxialAnisotropyField(mesh,
    // material)));

    // TODO//run with demag?
    // llgterm.push_back( llgt_ptr (new AtomisticDipoleDipoleField(mesh)));
    llgterm.push_back(llgt_ptr(new AtomisticExchangeField(mesh)));
    llgterm.push_back(llgt_ptr(new AtomisticDmiField(mesh, material)));
    llgterm.push_back(
        llgt_ptr(new AtomisticUniaxialAnisotropyField(mesh, material)));

    LLG Llg(state, llgterm);

    timer t = af::timer::start();
    while (state.t < 8.e-10) {
        state.m = Llg.step(state);
    }
    double timerelax = af::timer::stop(t);
    vti_writer_atom(state.m, mesh, filepath + "relax");

    std::cout << "timerelax [af-s]: " << timerelax << " for "
              << Llg.counter_accepted + Llg.counter_reject << " steps, thereof "
              << Llg.counter_accepted << " Steps accepted, "
              << Llg.counter_reject << " Steps rejected" << std::endl;

    array last = constant(0, mesh.dims, f64);
    last(span, span, span, 2) = 1;

    std::vector<State> inputimages;
    // inputimages.push_back(state);

    for (int i = 0; i < mesh.n0; i++) {
        array mm = array(state.m);
        mm = shift(mm, i);
        mm(seq(0, i), span, span, span) = 0;
        mm(seq(0, i), span, span, 2) = 1.;
        inputimages.push_back(State(mesh, material, mm));
    }

    // inputimages.push_back(State(mesh, material, last));

    String string(state, inputimages, n_interp, string_dt, llgterm);

    string.write_vti(filepath + "init_string");
    // for(unsigned j = 0; j < string.images.size(); j++){
    //    std::string name = filepath;
    //    name.append("init_string");
    //    name.append(std::to_string(j));
    //    vti_writer_atom(string.images[j].m, mesh , name.c_str());
    //}

    // String* string = new String(state, inputimages, n_interp , llgterm);
    std::cout.precision(12);

    std::ofstream stream_E_barrier;
    stream_E_barrier.precision(12);

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
        // af::printMemInfo();
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

        if (i > 25 && fabs(2 * (*max - string.E[0] - max_prev_step) /
                           (*max - string.E[0] + max_prev_step)) < 1e-6) {
            std::cout << "Exiting loop: Energy barrier relative difference "
                         "smaller than 1e-6"
                      << std::endl;
            stream_steps << "#Exiting loop: Energy barrier relative difference "
                            "smaller than 1e-6"
                         << std::endl;
            break;
        }
        if (i > 25 && fabs(*max - string.E[0] - max_prev_step) < 1e-27) {
            std::cout
                << "Exiting loop: Energy barrier difference smaller than 1e-27"
                << std::endl;
            stream_steps
                << "#Exiting loop: Energy barrier difference smaller than 1e-27"
                << std::endl;
            break;
        }
        std::cout << i << "\t" << *max - string.E[0] << "\t" << string.E[0]
                  << "\t" << *max - string.E[-1] << "\t" << *max << "\t"
                  << fabs(2 * (*max - string.E[0] - max_prev_step) /
                          (*max - string.E[0] + max_prev_step))
                  << std::endl;
        stream_steps << i << "\t" << *max - string.E[0] << "\t" << string.E[0]
                     << "\t" << *max - string.E[-1] << "\t" << *max << "\t"
                     << fabs(2 * (*max - string.E[0] - max_prev_step) /
                             (*max - string.E[0] + max_prev_step))
                     << std::endl;
        stream_E_barrier.open((filepath + "E_barrier.dat").c_str());
        stream_E_barrier << max_lowest << "\t" << nx << "\t" << dx << "\t"
                         << material.D << "\t" << material.Ku1 << "\t"
                         << material.K_atom << "\t" << material.D_atom
                         << std::endl;
        stream_E_barrier.close();
        for (unsigned j = 0; j < string.E.size(); ++j) {
            stream_E_curves << i << " " << j << " " << string.E[j] - string.E[0]
                            << " " << string.E[j] - string.E[-1] << " "
                            << string.E[j] << std::endl;
        }
        stream_E_curves << i << "\n\n" << std::endl;
        max_prev_step = *max - string.E[0];
        if (i % 20 == 1) {
            string.write_vti(filepath + "current_skyrm_image");
            // std::cout<<"Writing current skyrm images for
            // iteration"<<i<<std::endl; for(unsigned j = 0; j <
            // string.images.size(); j++){
            //    std::string name = filepath;
            //    name.append("current_skyrm_image");
            //    name.append(std::to_string(j));
            //    vti_writer_atom(string.images[j].m, mesh , name.c_str());
            //}
        }
    }
    std::cout << "#i , lowest overall:   max-[0], max-[-1] max [J]: "
              << i_max_lowest << "\t" << max_lowest << "\t"
              << max_lowest + E_max_lowest[0] - E_max_lowest[-1] << "\t"
              << max_lowest + E_max_lowest[0] << std::endl;
    stream_steps << "#i , lowest overall:   max-[0], max-[-1] max [J]: "
                 << i_max_lowest << "\t" << max_lowest << "\t"
                 << max_lowest + E_max_lowest[0] - E_max_lowest[-1] << "\t"
                 << max_lowest + E_max_lowest[0] << std::endl;
    stream_E_barrier.open((filepath + "E_barrier.dat").c_str());
    stream_E_barrier << max_lowest << "\t" << nx << "\t" << dx << "\t"
                     << material.D << "\t" << material.Ku1 << "\t"
                     << material.K_atom << "\t" << material.D_atom << std::endl;
    stream_E_barrier.close();
    // OLD:
    // stream_E_barrier<<max_lowest<<"\t"<<nx<<"\t"<<dx<<"\t"<<material.D<<"\t"<<material.Ku1<<"\t"<<material.D_atom<<"\t"<<material.K_atom<<"\t"<<std::endl;

    std::ofstream myfileE;
    myfileE.precision(12);
    myfileE.open((filepath + "E_last_step.dat").c_str());

    std::ofstream stream_max_lowest;
    stream_max_lowest.precision(12);
    stream_max_lowest.open((filepath + "E_max_lowest.dat").c_str());

    std::cout << string.E.size() << "\t" << string.images.size() << "\t"
              << std::endl;
    string.write_vti(filepath + "final_skyrm_image");
    for (unsigned i = 0; i < string.images.size(); i++) {
        std::cout << "i=" << i << "\t"
                  << "E= " << string.E[i] << std::endl;
        myfileE << i << "\t" << string.E[i] << "\t" << string.E[i] - string.E[0]
                << "\t" << string.E[i] - string.E[-1] << std::endl;
        stream_max_lowest << i << "\t" << E_max_lowest[i] << "\t"
                          << E_max_lowest[i] - E_max_lowest[0] << "\t"
                          << E_max_lowest[i] - E_max_lowest[-1] << std::endl;
        std::string name = filepath;
        name.append("skyrm_image_max_lowest");
        name.append(std::to_string(i));
        vti_writer_atom(images_max_lowest[i].m, mesh, name.c_str());
    }

    for (unsigned i = 0; i < Llg.Fieldterms.size(); ++i) {
        std::cout << "get_cpu_time()" << std::endl;
        std::cout << i << "\t" << Llg.cpu_time() << std::endl;
        stream_steps << "#"
                     << "get_cpu_time()" << std::endl;
        stream_steps << "#" << i << "\t" << Llg.cpu_time() << std::endl;
    }

    myfileE.close();
    stream_steps.close();
    stream_E_curves.close();
    stream_max_lowest.close();
    // delete[] string;

    return 0;
}
