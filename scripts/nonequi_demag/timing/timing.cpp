#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace magnumafcpp;

int main(int argc, char **argv)
{
    af::timer timer = af::timer::start();
    // Checking input variables and setting GPU Device
    for (int i = 0; i < argc; i++)
    {
        std::cout << "Parameter " << i << " was " << argv[i] << std::endl;
    }
    std::string filepath(argc > 1 ? argv[1] : "output_magnum.af/");
    af::setDevice(argc > 2 ? std::stoi(argv[2]) : 0);
    af::info();
    std::cout << "af setup [af-s]: " << af::timer::stop(timer) << std::endl;

    const bool sync = true;

    timer = af::timer::start();
    // Parameter initialization
    const double x = 5.e-7, y = 1.25e-7, z = 3.e-9;
    const int nx = 100, ny = 25;
    const int loops(argc > 3 ? std::stoi(argv[3]) : 100);
    const int nnz = (argc > 4 ? std::stoi(argv[4]) : 10);

    //Generating Objects

    std::ofstream stream;
    stream.precision(12);
    stream.open(filepath + "timing.dat");
    stream << "# info: timings for heff evaluations non-equidistant and equidistant evaluations" << std::endl;
    stream << "#" << loops << " evaluations of heff" << std::endl;
    stream << "# nz     mean_eq          stdev_eq        mean_ne         stdev_ne" << std::endl;

    //equidistant mesh
    for (int nz = 1; nz < nnz; nz++)
    {
        std::cout << "nz = " << nz << ", ";
        stream << nz << "\t";
        af::array m = Mesh(nx, ny, nz, x / nx, y / ny, z / nz).init_sp4();
        {
            Mesh mesh(nx, ny, nz, x / nx, y / ny, z / nz);
            State state(mesh, 8e5, m);
            DemagField demag(mesh, false, true, 0);

            timer = af::timer::start();
            demag.h(state); //init run
            if (sync)
                af::sync();
            std::cout << "eqJIT [s] = " << af::timer::stop(timer) << ", ";

            std::vector<double> times;
            timer = af::timer::start();
            for (int i = 0; i < loops; i++)
            {
                af::timer timer_loop = af::timer::start();
                af::array h = demag.h(state);
                if (sync)
                    af::eval(h);
                if (sync)
                    af::sync();
                times.push_back(af::timer::stop(timer_loop));
            }
            if (sync)
                af::sync();
            std::cout << loops << " loops in [s]: " << af::timer::stop(timer);
            std::cout << ", max [s] " << *std::max_element(times.begin(), times.end());
            std::cout << ", min [s] " << *std::min_element(times.begin(), times.end());

            auto m_stdev = mean_stdev_w_minus(times);
            std::cout << ", mean = " << m_stdev.first;
            std::cout << ", stdev = " << m_stdev.second;
            stream << m_stdev.first << "\t" << m_stdev.second << "\t";
        }

        //nonequidistant mesh
        {
            std::vector<double> z_spacing;
            for (int i = 0; i < nz; ++i)
            {
                z_spacing.push_back(z / nz);
            }
            NonequispacedMesh mesh_ne(nx, ny, x / nx, y / ny, z_spacing);
            State state_ne(mesh_ne, 8e5, m, false, true);
            NonEquiDemagField demag_ne = NonEquiDemagField(mesh_ne, false, false, 0);
            std::vector<double> times;
            const bool sync = true;

            timer = af::timer::start();
            demag_ne.h(state_ne); //init run
            if (sync)
                af::sync();
            std::cout << ", neJIT [s] = " << af::timer::stop(timer) << ", ";

            timer = af::timer::start();
            for (int i = 0; i < loops; i++)
            {
                af::timer timer_loop = af::timer::start();
                af::array h = demag_ne.h(state_ne);
                if (sync)
                    af::eval(h);
                if (sync)
                    af::sync();
                times.push_back(af::timer::stop(timer_loop));
            }
            if (sync)
                af::sync();
            std::cout << loops << " loops in [af-s]: " << af::timer::stop(timer);
            std::cout << ", max [s] " << *std::max_element(times.begin(), times.end());
            std::cout << ", min [s] " << *std::min_element(times.begin(), times.end());

            auto m_stdev = mean_stdev_w_minus(times);
            std::cout << ", mean = " << m_stdev.first;
            std::cout << ", stdev = " << m_stdev.second << std::endl;
            stream << m_stdev.first << "\t" << m_stdev.second << std::endl;
        }
    } //for
    stream.close();

    return 0;
}
