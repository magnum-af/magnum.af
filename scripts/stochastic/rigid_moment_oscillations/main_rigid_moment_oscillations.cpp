//Following paper Journal of Magnetism and Magnetic Materials 233 (2001) 296–304 from Scholz, Schrefl, Fidler
//ATTENTION, ERROR IN PAPER EQ 22 second expression, should be: (epx(chi)-1)/(sqrt(pi*chi)*erfi(sqrt(chi)))
//Integrals in EQ 22 only consider positive z, so we take fabs(mz) ! (this compensates switching which would lead to an average around zero)
#include <iostream>
#include <complex>
#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace magnumafcpp;

using namespace af;
using namespace std::complex_literals;
using Faddeeva::erfi;
typedef std::shared_ptr<LLGTerm> llgt_ptr;

//Mathematica:
//(e^x-1)/(sqrt(pi*x)*erfi(sqrt(x))) =(int(exp(x * z^2)*z) dz from 0 to 1 )/(int(exp(x * z^2)) dz from 0 to 1)
double mean_mz_analytical(double chi)
{
    return (exp(chi) - 1.) / (sqrt(M_PI) * sqrt(chi) * erfi(sqrt(chi))); //TODO erfi
}

void set_m_to_z(State &state)
{
    state.m(span, span, span, 0) = 0.;
    state.m(span, span, span, 1) = 0.;
    state.m(span, span, span, 2) = 1.;
}
void calcm(State state, std::ostream &myfile)
{
    myfile << std::setw(12) << state.t << "\t" << meani(state.m, 0) << "\t" << meani(state.m, 1) << "\t" << meani(state.m, 2) << "\t" << std::endl;
}

//This should handle other cases as well
bool signbit(double a)
{
    if (a >= 0)
        return true;
    if (a < 0)
        return false;
}

int main(int argc, char **argv)
{
    std::cout << "argc" << argc << std::endl;
    for (int i = 0; i < argc; i++)
    {
        std::cout << "Parameter " << i << " was " << argv[i] << "\n";
    }
    std::string filepath(argc > 1 ? argv[1] : "../Data/rigid");
    if (argc > 0)
        filepath.append("/");
    std::cout << "Writing into path " << filepath.c_str() << std::endl;
    std::ofstream stream;
    stream.precision(12);
    stream.open((filepath + "m.dat").c_str());

    setDevice(argc > 2 ? std::stoi(argv[2]) : 0);
    info();

    // Parameter initialization
    const double x = 1.e-9, y = 1.e-9, z = 1.e-9;
    const int nx = 1, ny = 1, nz = 1;

    //Generating Objects
    Mesh mesh(nx, ny, nz, x / nx, y / ny, z / nz);
    Material material = Material();
    state.Ms = 1281197;
    material.Ku1 = 6.9e6;
    material.alpha = 0.1;
    material.T = 0;

    // Initial magnetic field
    array m = constant(0.0, mesh.n0, mesh.n1, mesh.n2, 3, f64);
    const double tile = M_PI / 4.;
    m(0, 0, 0, 0) = sin(tile);
    m(0, 0, 0, 1) = 0.;
    m(0, 0, 0, 2) = cos(tile);
    State state(mesh, material, m); //ATTENTION, to be set in each loop
    std::vector<llgt_ptr> llgterm;
    Stochastic_LLG Stoch(state, llgterm, 0., "Heun"); //ATTENTION, to be set in each loop

    double dt = 5e-16;

    std::cout << "T_analytic (=1/f) = " << 2 * M_PI / material.gamma / (2 * material.Ku1 / (constants::mu0 * state.Ms)) << " [s]" << std::endl;

    int count_zero_x = 0;
    double prev_x = 0;
    double prev_x_t = 0;
    int count_zero_y = 0;
    double prev_y = 0;
    double prev_y_t = 0;

    state = State(mesh, material, m);
    llgterm.push_back(llgt_ptr(new UniaxialAnisotropyField(mesh, material)));
    Stoch = Stochastic_LLG(state, llgterm, dt, "Heun");
    while (state.t < 30e-12)
    {
        prev_x = meani(state.m, 0);
        prev_y = meani(state.m, 1);
        Stoch.step(state, dt);
        calcm(state, stream);
        if (signbit(meani(state.m, 0)) != signbit(prev_x))
        {
            count_zero_x++;
            if (count_zero_x % 2 == 0)
            {
                if (count_zero_x > 2)
                    std::cout << count_zero_x / 2 - 1 << "-" << count_zero_x / 2 << " x crossing with  " << state.t - prev_x_t << " [s]" << std::endl;
                prev_x_t = state.t;
            }
            //std::cout << "x Zerocrossing at t= " << state.t << std::endl;
        }
        if (signbit(meani(state.m, 1)) != signbit(prev_y))
        {
            count_zero_y++;
            if (count_zero_y % 2 == 0)
            {
                std::cout << count_zero_y / 2 - 1 << "-" << count_zero_y / 2 << " y crossing with  " << state.t - prev_y_t << " [s]" << std::endl;
                prev_y_t = state.t;
            }
            //std::cout << "y Zerocrossing at t= " << state.t << std::endl;
        }
        //mean_mz+=afvalue(state.m(0, 0, 0, 2));
        //abs_mean_mz+=fabs(afvalue(state.m(0, 0, 0, 2)));
    }
    std::cout << "Zerocrossings in x = " << count_zero_x << std::endl;
    std::cout << "Zerocrossings in y = " << count_zero_y << std::endl;
    stream.close();
    return 0;
}
