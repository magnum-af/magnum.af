#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace magnumaf;


int main(int argc, char** argv)
{
    af::timer timer = af::timer::start();
    // Checking input variables and setting GPU Device
    for (int i=0; i<argc; i++){std::cout << "Parameter " << i << " was " << argv[i] << std::endl;}
    std::string filepath(argc>1? argv[1]: "output_magnum.af/");
    af::setDevice(argc>2? std::stoi(argv[2]):0);
    af::info();
    std::cout << "af setup [af-s]: " << af::timer::stop(timer) << std::endl;

    timer = af::timer::start();
    // Parameter initialization
    const double x=5.e-7, y=1.25e-7, z=3.e-9;
    const int nx = 100, ny=25 , nz=1;

    //Generating Objects
    Mesh mesh(nx, ny, nz, x/nx, y/ny, z/nz);
    State state(mesh, 8e5, mesh.init_sp4());
    DemagField demag(mesh, true, true, 0);
    std::cout << "demag setup [af-s]: " << af::timer::stop(timer) << std::endl;
    //af::print("demag", demag.h(state));

    //const int loops = 100;
    //std::array<double, loops> times = {0};
    //double times[loops] = {0};
    const int loops(argc>3? std::stoi(argv[3]):100);
    std::vector<double> times;
    const bool sync = true;
    for(int i = 0; i<loops; i++){
        af::timer timer_loop = af::timer::start();
        af::array h = demag.h(state);
        if(sync) af::sync();
        times.push_back(af::timer::stop(timer_loop));
        //times[i] = af::timer::stop(timer_loop);
        std::cout << i <<", h timing [af-s]: " << times[i] << std::endl;
    }
    if(sync) af::sync();
    std::cout << "h total [af-s]: " << af::timer::stop(timer) << std::endl;
    std::cout << "max [s] " << *std::max_element(times.begin(), times.end()) << std::endl;
    std::cout << "min [s] " << *std::min_element(times.begin(), times.end()) << std::endl;


    auto m_stdev = mean_stdev_w_minus(times);
    std::cout << "mean = " << m_stdev.first << std::endl;
    std::cout << "stdev = " << m_stdev.second << std::endl;

    return 0;
}
