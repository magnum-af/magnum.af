#include <list>
#include <sstream>
#include "arrayfire.h"
#include "magnum_af.hpp"

int main(int argc, char **argv)
{
    std::cout << "argc = " << argc << std::endl;
    for (int i = 0; i < argc; i++)
    {
        cout << "Parameter " << i << " was " << argv[i] << "\n";
    }

    std::string inputfile(argc > 1 ? argv[1] : "../Data/skyrmion_stoch/Testing");
    std::string outputfile(argc > 2 ? argv[2] : "");
    std::string outputfile2(argc > 3 ? argv[3] : "");

    array m;
    Mesh testmesh(0, 0, 0, 0, 0, 0);
    vti_reader(m, testmesh, inputfile);

    pgfplot_mz(m, testmesh, outputfile, 1);
    pgfplot_nz_quiver(m, testmesh, outputfile2, 1, 6);
    //pgfplot_mi(m, testmesh, outputfile, 0, 0, 2);
    //pgfplot_nz_quiver(m, testmesh, outputfile2, 0, 8);

    return 0;
}
