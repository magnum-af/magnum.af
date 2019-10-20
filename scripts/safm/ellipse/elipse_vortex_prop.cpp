#include "magnum_af.hpp"

using namespace magnumaf;

void print_m(const af::array& m, std::string name, std::ostream& stream = std::cout){
    af::array m_mean = af::mean(af::mean(af::mean(m, 0), 1), 2);
    stream << name << ": " << afvalue(m_mean(0, 0, 0, 0)) << ", " << afvalue(m_mean(0, 0, 0, 1)) << ", " << afvalue(m_mean(0, 0, 0, 2)) <<  std::endl;
}
void print_m(const af::array& m, std::ostream& stream = std::cout, float factor = 1){
    af::array m_mean = af::mean(af::mean(af::mean(m, 0), 1), 2); stream << std::fixed << factor * afvalue(m_mean(0, 0, 0, 0)) << "\t" << std::fixed << factor * afvalue(m_mean(0, 0, 0, 1)) << "\t" << std::fixed << factor * afvalue(m_mean(0, 0, 0, 2));
}

void print_layers(const af::array& m){
    print_m(m(af::span, af::span, 0, af::span), "safm lay0");
    print_m(m(af::span, af::span, 1, af::span), "vac  lay1");
    print_m(m(af::span, af::span, 2, af::span), "safm lay2");
    print_m(m(af::span, af::span, 3, af::span), "vac  lay3");
    print_m(m(af::span, af::span, 4, af::span), "freelayer");
    std::cout << std::endl;
}

int main(int argc, char** argv)
{
    // Checking input variables and setting GPU Device
    StageTimer timer;
    for (int i=0; i<argc; i++){std::cout << "Parameter " << i << " was " << argv[i] << std::endl;}
    std::string filepath = std::string(argc>1? argv[1]: "output_magnum.af/");
    af::setDevice(argc>2? std::stoi(argv[2]):0);
    af::info();

    std::cout.precision(8);

    // Parameter initialization
    // TODO vals?
    const float x=5000.e-9, y=1000e-9, z=4.e-9;
    const int nx = 5*32, ny=32 , nz=5;
    //const float Ms = 8e5;//TODO 0.5/Constants.mu0
    //const float A = 15e-12;
    float Ms    = 1.393e6;//[J/T/m^3] == [Joule/Tesla/meter^3] = 1.75 T/mu_0
    float A     = 1.5e-11;//[J/m]
    //TODO
    const float Ku1_safm=1e6; // Anisotropy constant
    const float Ku1_freelayer=79e+03; // Anisotropy constant

    const int zee_dir = (argc > 3 ? std::stoi(argv[3]) : 2);// 0 == x, 2 == z
    const std::string int_over_min = (argc > 4 ? std::string(argv[4]) : "int");// true: use integrate, false: use minimizer
    //const bool int_over_min = (argc > 4 ? std::stob(argv[4]) : true);// true: use integrate, false: use minimizer
    const float hzee_max = (argc > 5 ? std::stod(argv[5]): 0.12); //[Tesla]
    const float integr_time = (argc > 6 ? std::stod(argv[6]): 4 * 10e-11); //[s]
    const int quater_steps =(argc > 7 ? std::stoi(argv[7]) : 10);
    //af::timer timer = af::timer::start();
    //Generating Objects
    //std::vector<float> z_spacing = {z/4., 1e-10, 9.84825e-10, z/4., z/4.};
    float dz = 10e-9;
    std::vector<float> z_spacing = {dz, 1e-10, dz, dz, dz};
    //std::vector<float> z_spacing = {z/nz, 9.84825e-10, z/nz, z/nz};
    NonequispacedMesh mesh(nx, ny, x/nx, y/ny, z_spacing);
    mesh.print();

    // Initial magnetic field
    // Returns an initial elliptical magnetization
    // n_cells gives number of cells with non-zero Ms
    // xyz gives direction of initial magnetization direction,
    // positive_direction true points +, false in - direction
    af::array m = af::constant(0.0, nx, ny, nz, 3, f32);
    af::array Ku1_field = af::constant(0, mesh.dims, f32);
    for(int ix=0;ix<nx;ix++){
        for(int iy=0;iy<ny;iy++){
            const float a= (float)(nx/2);
            const float b= (float)(ny/2);
            const float rx=float(ix)-nx/2.;
            const float ry=float(iy)-ny/2.;
            const float r = pow(rx, 2)/pow(a, 2)+pow(ry, 2)/pow(b, 2);
            if(r<1){
                for(int iz=0;iz<nz;iz++){
                }
                //if(positive_direction) m(ix, iy, af::span, xyz)=1;
                //else m(ix, iy, af::span, xyz)=-1;

                m(ix, iy, 0, 2) = 1;//SAFM Layer 0 in z
                m(ix, iy, 2, 2) = -1;//SAFM Layer 1 in -z
                //TD//m(ix, iy, 4, 2) = 1;// Free Layer in z
                m(ix, iy, 4, 0) = 1;// Free Layer in x
                Ku1_field(ix, iy, 0, af::span) = Ku1_safm;//SAFM Layer 0 in z
                Ku1_field(ix, iy, 2, af::span) = Ku1_safm;//SAFM Layer 1 in -z
                Ku1_field(ix, iy, 4, af::span) = Ku1_freelayer;// Free Layer in x
            }
        }
    }
    std::cout << "Info: Mesh::ellipse(): n_cells should be approx a*b*M_PI*this->n2= " << nx/2*ny/2*M_PI*nz << std::endl;
    //af::array m = af::constant(0, mesh.dims, f32);

    State state(mesh, Ms, m, false);
    state.mesh = Mesh(nx, ny, nz, x/nx, y/ny, z/nz);//TODO necessary for external field?
    //vtr_writer(state.m, Mesh(nx, ny, nz, x/nx, y/ny, 0) , z_spacing, filepath + "m_init");
    state.vtr_writer(filepath + "m_init");
    vtr_writer(Ku1_field, state.nonequimesh, filepath + "m_init");
    timer.print_stage("mesh  ");

    auto demag = LlgTerm (new NonEquiDemagField(mesh, true, true, 0));
    auto exch = LlgTerm (new NonequiExchangeField(A, mesh, true));
    //TODO//auto aniso = LlgTerm (new NonequiUniaxialAnisotropyField(Ku1_field, std::array<float, 3>{0, 0, 1}));
    //LlgTerms llgterms = {demag, exch, aniso};
    LlgTerms llgterms = {demag, exch};
    timer.print_stage("setup ");

    af::array h = demag->h(state);
    //vtr_writer(h, Mesh(nx, ny, nz, x/nx, y/ny, 0) , z_spacing, filepath + "h");
    vtr_writer(h, state.nonequimesh, filepath + "h_init");

    std::ofstream stream;
    stream.precision(12);
    stream.open (filepath + "h.dat");
    stream << z_spacing[1] << ", " << afvalue(h(nx/2, ny/2, 3, 0)) << ", " << afvalue(h(nx/2, ny/2, 3, 1)) << ", " << afvalue(h(nx/2, ny/2, 3, 2)) << std::endl;
    stream.close();
    timer.print_stage("calc h");

    //integrate
    if(int_over_min == "int")
    {
        const float time_relax = 1e-10;

        LLGIntegrator llg(1, llgterms);
        //float relax_prec = 0.001;
        //llg.relax(state, relax_prec);
        //timer.print_stage("relax " + std::to_string(relax_prec));
        //std::cout << "relaxed" << std::endl;
        const bool relax = false;
        if(relax)
        {
            while (state.t < time_relax){
                llg.step(state);
                std::cout << state.t << "\t";
                print_m(state.m(af::span, af::span, -1, af::span), std::cout);
                std::cout << std::endl;
            }
            timer.print_stage("relax " + std::to_string(time_relax * 1e9) + " [ns]");
            //vtr_writer(state.m, Mesh(nx, ny, nz, x/nx, y/ny, 0) , z_spacing, filepath + "m_relax");
            state.vtr_writer(filepath + "m_relax");
            state.t = 0;
            state.steps = 0;
        }
        {
            //float rate = 2e9 ; //[T/s]
            float integr_time_per_quater = integr_time/4.;
            float rate = hzee_max/integr_time_per_quater; //[T/s]
            std::cout << "hzee_max= " << hzee_max << ", rate=" << rate << ", integr_time_per_quater=" << integr_time_per_quater << std::endl;
            //float rate = 0.34e6 ; //[T/s]
            auto zee_func_llg= [ hzee_max, rate, zee_dir ] ( State state ) -> af::array {
                float field_Tesla = 0;
                if(state.t < hzee_max/rate) field_Tesla = rate *state.t;
                else if (state.t < 3*hzee_max/rate) field_Tesla = -rate *state.t + 2*hzee_max;
                else if(state.t < 5*hzee_max/rate) field_Tesla = rate*state.t - 4*hzee_max;
                //else if(state.t < 4*hzee_max/rate) field_Tesla = rate*state.t - 4*hzee_max;
                else {field_Tesla =  rate*state.t - 4*hzee_max; std::cout << "NOTE: zee time out of range" << std::endl;}
                af::array zee = af::constant(0.0, state.mesh.n0, state.mesh.n1, state.mesh.n2, 3, f32);
                zee(af::span, af::span, af::span, zee_dir) = af::constant(field_Tesla/constants::mu0 , state.mesh.n0, state.mesh.n1, state.mesh.n2, 1, f32);
                return  zee;
            };

            llg.llgterms.push_back( LlgTerm (new ExternalField(zee_func_llg)));
            //llg.llgterms.push_back( LlgTerm (new NonequiExternalField(zee_func_llg)));
            stream.open(filepath + "m.dat");
            stream << "# step	<mx>    <my>    <mz>    hx      hy      hz" << std::endl;
            while (state.t < 5* hzee_max/rate){
            //while (state.t < 4* hzee_max/rate){
                llg.step(state);
                if (state.steps % 1 == 0){
                    af::array h_zee = llg.llgterms.back()->h(state) * constants::mu0;//in Tesla only for output
                    std::cout << "t[ns]= " << std::fixed << state.t * 1e9 << ", mx= "<< std::fixed ;
                    print_m(state.m(af::span, af::span, -1, af::span), std::cout, 1/(0.25 * M_PI));
                    std::cout << ", safm_l2=";
                    print_m(state.m(af::span, af::span, 2, af::span), std::cout, 1/(0.25 * M_PI));
                    std::cout << ", zee=" << std::fixed << afvalue(h_zee(0, 0, 0, 0)) << "\t" << std::fixed  << afvalue(h_zee(0, 0, 0, 1)) << "\t" << std::fixed << afvalue(h_zee(0, 0, 0, 2));
                    std::cout << std::endl;
                    //stream << state.steps << "\t";
                    stream << state.t << "\t";
                    print_m(state.m(af::span, af::span, -1, af::span), stream, 1/(0.25 * M_PI));
                    stream << "\t" << afvalue(h_zee(0, 0, 0, 0)) << ", "  << afvalue(h_zee(0, 0, 0, 1)) << ", " << afvalue(h_zee(0, 0, 0, 2)) << std::endl;
                }
                if (state.steps % 100 == 0){
                    state.vtr_writer(filepath + "m_int" + std::to_string(state.steps));
                }
                //if(state.steps % 1 == 0){vtr_writer(state.m, Mesh(nx, ny, nz, x/nx, y/ny, 0) , z_spacing, filepath + "m_int" + std::to_string(state.steps));}
            }
            timer.print_stage("hys " + std::to_string(state.t) + " [ns]");
            state.t = 0;
            state.steps = 0;
        }
    }

    // minimize
    else
    {
        LBFGS_Minimizer minimizer(llgterms, 1e-6, 230, 0);
        minimizer.Minimize(state);
        print_layers(state.m);
        //vtr_writer(state.m, Mesh(nx, ny, nz, x/nx, y/ny, 0) , z_spacing, filepath + "m_minimzed");
        state.vtr_writer(filepath + "m_minimized");

        auto zee_func= [ hzee_max, quater_steps, zee_dir ] ( State state ) -> af::array {
            float field_Tesla = 0;
            float rate = hzee_max/quater_steps; //[T/s]
            if(state.t < hzee_max/rate) field_Tesla = rate *state.t;
            else if (state.t < 3*hzee_max/rate) field_Tesla = -rate *state.t + 2*hzee_max;
            else if(state.t < 4*hzee_max/rate) field_Tesla = rate*state.t - 4*hzee_max;
            else {field_Tesla = 0; std::cout << "WARNING ZEE time out of range" << std::endl;}
            af::array zee = af::constant(0.0, state.mesh.n0, state.mesh.n1, state.mesh.n2, 3, f32);
            zee(af::span, af::span, af::span, zee_dir) = af::constant(field_Tesla/constants::mu0 , state.mesh.n0, state.mesh.n1, state.mesh.n2, 1, f32);
            return  zee;
        };

        //hys loop
        state.t = 0;
        state.steps = 0;
        minimizer.llgterms_.push_back( LlgTerm (new ExternalField(zee_func)));
        stream.precision(12);
        stream.open(filepath + "m.dat");
        stream << "# step	<mx>    <my>    <mz>    hx      hy      hz" << std::endl;
        float rate = hzee_max/quater_steps; //[T/s]
        while (state.t < 4* hzee_max/rate){
            minimizer.Minimize(state);
            //state.calc_mean_m(stream, afvalue(minimizer.llgterms_[minimizer.llgterms_.size()-1]->h(state)(0, 0, 0, 2)));
            af::array h_zee = minimizer.llgterms_[minimizer.llgterms_.size()-1]->h(state) * constants::mu0;
            std::cout << state.steps << "\t";
            print_m(state.m(af::span, af::span, -1, af::span), std::cout);
            std::cout << ", zee=" << afvalue(h_zee(0, 0, 0, 0)) << ", "  << afvalue(h_zee(0, 0, 0, 1)) << ", " << afvalue(h_zee(0, 0, 0, 2)) << std::endl;
            stream << state.steps << "\t";
            print_m(state.m, stream);
            stream << "\t" << afvalue(h_zee(0, 0, 0, 0)) << ", "  << afvalue(h_zee(0, 0, 0, 1)) << ", " << afvalue(h_zee(0, 0, 0, 2)) << std::endl;
            if(state.steps % 1 == 0){
                state.vtr_writer(filepath + "m_minimzed" + std::to_string(state.steps));
                //vtr_writer(state.m, Mesh(nx, ny, nz, x/nx, y/ny, 0) , z_spacing, filepath + "m_minimzed" + std::to_string(state.steps));
            }
            state.steps++;
            state.t+=1.;
        }
    }


    timer.print_accumulated();
    return 0;
}
