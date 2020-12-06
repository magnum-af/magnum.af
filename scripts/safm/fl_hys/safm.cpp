#include "magnum_af.hpp"

using namespace magnumafcpp;

void print_m(const af::array& m, std::string name, std::ostream& stream = std::cout) {
    af::array mean_m = af::mean(af::mean(af::mean(m, 0), 1), 2);
    stream << name << ": " << afvalue(mean_m(0, 0, 0, 0)) << ", " << afvalue(mean_m(0, 0, 0, 1)) << ", "
           << afvalue(mean_m(0, 0, 0, 2)) << std::endl;
}

void print_m(const af::array& m, std::ostream& stream = std::cout) {
    af::array mean_m = af::mean(af::mean(af::mean(m, 0), 1), 2);
    stream << std::fixed << afvalue(mean_m(0, 0, 0, 0)) << "\t" << std::fixed << afvalue(mean_m(0, 0, 0, 1)) << "\t"
           << std::fixed << afvalue(mean_m(0, 0, 0, 2));
}

void print_layers(const af::array& m) {
    print_m(m(af::span, af::span, 0, af::span), "safm lay0");
    print_m(m(af::span, af::span, 1, af::span), "vac  lay1");
    print_m(m(af::span, af::span, 2, af::span), "safm lay2");
    print_m(m(af::span, af::span, 3, af::span), "vac  lay3");
    print_m(m(af::span, af::span, 4, af::span), "freelayer");
    std::cout << std::endl;
}

int main(int argc, char** argv) {
    // Checking input variables and setting GPU Device
    StageTimer timer;
    for (int i = 0; i < argc; i++) {
        std::cout << "Parameter " << i << " was " << argv[i] << std::endl;
    }
    std::string filepath = std::string(argc > 1 ? argv[1] : "output_magnum.af/");
    af::setDevice(argc > 2 ? std::stoi(argv[2]) : 0);
    af::info();

    std::cout.precision(8);

    // Parameter initialization
    const double x = 40.e-9, y = 40e-9, z = 4.e-9;
    const int nx = 32, ny = 32, nz = 5;
    const double Ms = 8e5; // TODO 0.5/Constants.mu0
    const double A = 15e-12;
    const double Ku1_safm = 1e6;         // Anisotropy constant
    const double Ku1_freelayer = 79e+03; // Anisotropy constant

    const int zee_dir = (argc > 3 ? std::stoi(argv[3]) : 2); // 0 == x, 2 == z
    const std::string int_over_min =
        (argc > 4 ? std::string(argv[4]) : "int"); // true: use integrate, false: use minimizer
    // const bool int_over_min = (argc > 4 ? std::stob(argv[4]) : true);// true:
    // use integrate, false: use minimizer
    const double hzee_max = (argc > 5 ? std::stod(argv[5]) : 0.12); //[Tesla]
    // const double inttime = (argc > 6 ? std::stod(argv[6]): 10e-9); //[s]
    const int quater_steps = (argc > 7 ? std::stoi(argv[7]) : 10);
    // af::timer timer = af::timer::start();
    // Generating Objects
    std::vector<double> z_spacing = {z / 4., 1e-10, 9.84825e-10, z / 4., z / 4.};
    // std::vector<double> z_spacing = {z/nz, 9.84825e-10, z/nz, z/nz};
    NonequiMesh mesh(nx, ny, x / nx, y / ny, z_spacing);
    std::cout << mesh << std::endl;

    // Initial magnetic field
    af::array m = af::constant(0, dims_vector(mesh), f64);
    m(af::span, af::span, 0, 2) = 1;        // SAFM Layer 0 in z
    m(af::span, af::span, 1, 2) = 0;        // Vacuum Layer 1
    m(af::span, af::span, 2, 2) = -1;       // SAFM Layer 1 in -z
    m(af::span, af::span, 3, af::span) = 0; // Vacuum Layer 2
    m(af::span, af::span, 4, 0) = 1;        // Free Layer in x

    af::array Ku1_field = af::constant(0, dims_vector(mesh), f64);
    Ku1_field(af::span, af::span, 0, af::span) = Ku1_safm;      // SAFM Layer 0 in z
    Ku1_field(af::span, af::span, 1, af::span) = 0;             // Vacuum Layer
    Ku1_field(af::span, af::span, 2, af::span) = Ku1_safm;      // SAFM Layer 1 in -z
    Ku1_field(af::span, af::span, 3, af::span) = 0;             // Vacuum Layer 2
    Ku1_field(af::span, af::span, 4, af::span) = Ku1_freelayer; // Free Layer in x
    State state(mesh, Ms, m, false);
    state.mesh = Mesh(nx, ny, nz, x / nx, y / ny, z / nz);
    vtr_writer(state.m, Mesh(nx, ny, nz, x / nx, y / ny, 0), z_spacing, filepath + "m_init");
    timer.print_stage("mesh  ");

    auto demag = LlgTerm(new NonequiDemagField(mesh, true, true, 0));
    auto exch = LlgTerm(new NonequiExchangeField(mesh, A, true));
    auto aniso = LlgTerm(new NonequiUniaxialAnisotropyField(
        mesh, Ku1_field, std::array<double, 3>{0, 0, 1})); // TODO energy is wrong is its not nonequi!
    vec_uptr_Fieldterm llgterms = {demag, exch, aniso};
    timer.print_stage("setup ");

    af::array h = demag->h(state);
    vtr_writer(h, Mesh(nx, ny, nz, x / nx, y / ny, 0), z_spacing, filepath + "h");

    std::ofstream stream;
    stream.precision(12);
    stream.open(filepath + "h.dat");
    stream << z_spacing[1] << ", " << afvalue(h(nx / 2, ny / 2, 3, 0)) << ", " << afvalue(h(nx / 2, ny / 2, 3, 1))
           << ", " << afvalue(h(nx / 2, ny / 2, 3, 2)) << std::endl;
    stream.close();
    timer.print_stage("calc h");

    // integrate
    if (int_over_min == "int") {
        const double time_relax = 1e-10;

        LLGIntegrator llg(1, llgterms);
        // double relax_prec = 0.001;
        // llg.relax(state, relax_prec);
        // timer.print_stage("relax " + std::to_string(relax_prec));
        // std::cout << "relaxed" << std::endl;
        const bool relax = false;
        if (relax) {
            while (state.t < time_relax) {
                llg.step(state);
                std::cout << state.t << "\t";
                print_m(state.m(af::span, af::span, -1, af::span), std::cout);
                std::cout << std::endl;
            }
            timer.print_stage("relax " + std::to_string(time_relax * 1e9) + " [ns]");
            vtr_writer(state.m, Mesh(nx, ny, nz, x / nx, y / ny, 0), z_spacing, filepath + "m_relax");
            state.t = 0;
            state.steps = 0;
        }
        {
            // double rate = 2e9 ; //[T/s]
            double int_time_per_quater = 10e-11;
            double rate = hzee_max / int_time_per_quater; //[T/s]
            std::cout << "hzee_max= " << hzee_max << ", rate=" << rate
                      << ", int_time_per_quater=" << int_time_per_quater << std::endl;
            // double rate = 0.34e6 ; //[T/s]
            auto zee_func_llg = [hzee_max, rate, zee_dir](State state) -> af::array {
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
                af::array zee = af::constant(0.0, state.mesh.nx, state.mesh.ny, state.mesh.nz, 3, f64);
                zee(af::span, af::span, af::span, zee_dir) =
                    af::constant(field_Tesla / constants::mu0, state.mesh.nx, state.mesh.ny, state.mesh.nz, 1, f64);
                return zee;
            };

            llg.llgterms.push_back(LlgTerm(new NonequiExternalField(mesh, zee_func_llg)));
            stream.open(filepath + "m_int.dat");
            stream << "# step	<mx>    <my>    <mz>    hx      hy      hz" << std::endl;
            while (state.t < 4 * hzee_max / rate) {
                llg.step(state);
                if (state.steps % 100 == 0) {
                    af::array h_zee = llg.llgterms.back()->h(state) * constants::mu0; // in Tesla only for output
                    std::cout << "t[ns]= " << std::fixed << state.t * 1e9 << ", mx= " << std::fixed;
                    print_m(state.m(af::span, af::span, -1, af::span), std::cout);
                    std::cout << ", zee=" << std::fixed << afvalue(h_zee(0, 0, 0, 0)) << "\t" << std::fixed
                              << afvalue(h_zee(0, 0, 0, 1)) << "\t" << std::fixed << afvalue(h_zee(0, 0, 0, 2));
                    std::cout << std::endl;
                    // stream << state.steps << "\t";
                    stream << state.t << "\t";
                    print_m(state.m, stream);
                    stream << "\t" << afvalue(h_zee(0, 0, 0, 0)) << ", " << afvalue(h_zee(0, 0, 0, 1)) << ", "
                           << afvalue(h_zee(0, 0, 0, 2)) << std::endl;
                }
                // if(state.steps % 1 == 0){vtr_writer(state.m, Mesh(nx, ny, nz,
                // x/nx, y/ny, 0) , z_spacing, filepath + "m_int" +
                // std::to_string(state.steps));}
            }
            timer.print_stage("hys " + std::to_string(state.t) + " [ns]");
            state.t = 0;
            state.steps = 0;
        }
    }

    // minimize
    else {
        LBFGS_Minimizer minimizer(llgterms, 1e-6, 230, 0);
        minimizer.Minimize(state);
        print_layers(state.m);
        vtr_writer(state.m, Mesh(nx, ny, nz, x / nx, y / ny, 0), z_spacing, filepath + "m_minimzed");

        auto zee_func = [hzee_max, quater_steps, zee_dir](State state) -> af::array {
            double field_Tesla = 0;
            double rate = hzee_max / quater_steps; //[T/s]
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
            af::array zee = af::constant(0.0, state.mesh.nx, state.mesh.ny, state.mesh.nz, 3, f64);
            zee(af::span, af::span, af::span, zee_dir) =
                af::constant(field_Tesla / constants::mu0, state.mesh.nx, state.mesh.ny, state.mesh.nz, 1, f64);
            return zee;
        };

        // hys loop
        state.t = 0;
        state.steps = 0;
        minimizer.llgterms_.push_back(LlgTerm(new NonequiExternalField(mesh, zee_func)));
        stream.precision(12);
        stream.open(filepath + "m.dat");
        stream << "# step	<mx>    <my>    <mz>    hx      hy      hz" << std::endl;
        double rate = hzee_max / quater_steps; //[T/s]
        while (state.t < 4 * hzee_max / rate) {
            minimizer.Minimize(state);
            // state.calc_mean_m(stream,
            // afvalue(minimizer.llgterms_[minimizer.llgterms_.size()-1]->h(state)(0,
            // 0, 0, 2)));
            af::array h_zee = minimizer.llgterms_[minimizer.llgterms_.size() - 1]->h(state) * constants::mu0;
            std::cout << state.steps << "\t";
            print_m(state.m(af::span, af::span, -1, af::span), std::cout);
            std::cout << ", zee=" << afvalue(h_zee(0, 0, 0, 0)) << ", " << afvalue(h_zee(0, 0, 0, 1)) << ", "
                      << afvalue(h_zee(0, 0, 0, 2)) << std::endl;
            stream << state.steps << "\t";
            print_m(state.m, stream);
            stream << "\t" << afvalue(h_zee(0, 0, 0, 0)) << ", " << afvalue(h_zee(0, 0, 0, 1)) << ", "
                   << afvalue(h_zee(0, 0, 0, 2)) << std::endl;
            if (state.steps % 1 == 0) {
                vtr_writer(state.m, Mesh(nx, ny, nz, x / nx, y / ny, 0), z_spacing,
                           filepath + "m_minimzed" + std::to_string(state.steps));
            }
            state.steps++;
            state.t += 1.;
        }
    }

    timer.print_accumulated();
    return 0;
}
