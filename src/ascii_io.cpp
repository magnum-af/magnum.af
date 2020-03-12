#include "ascii_io.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>

namespace magnumafcpp
{

//TODO .scalar<> is slow, use raw pointer
void write_ascii(const af::array& a, const Mesh& mesh, std::string filename, bool verbose, int precision)
{
    std::ofstream stream(filename);
    if( stream.is_open())
    {
        stream.precision(precision);
        int val_setw = precision + 4; // formatted at least for 0.*** values
        if (verbose) {printf("writing to ascii file %s with precision of %d", filename.c_str(), precision);}
        stream << "# " << mesh.n0 << " " << mesh.n1 << " " << mesh.n2 << " " << mesh.dx << " " << mesh.dy << " " << mesh.dz << std::endl;
        stream << "# ascii encoded vector field" << std::endl;
        stream << "# plot e.g. with gnuplot: splot '" << filename << "' using 1:2:3:4:5:6 with vectors" << std::endl;
        stream << "#                       : splot '" << filename << "' using 1:2:3:4:5:6:(sqrt($4**2+$5**2+$6**2)) with vectors linecolor palette z" << std::endl;
        stream << "#\n# ";
        stream << std::setw(val_setw) << std::left << "dx" << "\t";
        stream << std::setw(val_setw) << std::left << "dy" << "\t";
        stream << std::setw(val_setw) << std::left << "dz" << "\t";
        stream << std::setw(val_setw) << std::left << "mx" << "\t";
        stream << std::setw(val_setw) << std::left << "my" << "\t";
        stream << std::setw(val_setw) << std::left << "mz" << "\n";


        for(int ix = 0; ix < a.dims(0); ix++)
        {
            for(int iy = 0; iy < a.dims(1); iy++)
            {
                for(int iz = 0; iz < a.dims(2); iz++)
                {
                    stream << std::setw(val_setw) << std::left << ix * mesh.dx << "\t";
                    stream << std::setw(val_setw) << std::left << iy * mesh.dy << "\t";
                    stream << std::setw(val_setw) << std::left << iz * mesh.dz << "\t";
                    stream << std::setw(val_setw) << std::left << a(ix, iy, iz, 0).scalar<double>() << "\t";
                    stream << std::setw(val_setw) << std::left << a(ix, iy, iz, 1).scalar<double>() << "\t";
                    stream << std::setw(val_setw) << std::left << a(ix, iy, iz, 2).scalar<double>() << "\n";
                }
            }
        }
    }

}
void read_ascii(af::array& a, Mesh& mesh, std::string filename, bool verbose)
{
    std::ifstream infile(filename);
    if( infile.is_open())
    {
        if (verbose){
            printf("reading from ascii file %s\n", filename.c_str());
        }

        std::string line;

        // Getting mesh size and discretization
        uint32_t nx, ny, nz;
        double dx, dy, dz;
        std::getline(infile, line);
        line.erase(0, 2);
        std::istringstream iss1(line);
        if (!(iss1 >> nx >> ny >> nz >> dx >> dy >> dz)){printf("Error while reading line!\n");}
        Mesh read_mesh(nx, ny, nz, dx, dy, dz);
        mesh = read_mesh;
        af::array read_a = af::constant(0, nx, ny, nz, 3, f64);

        // skipping all further lines starting with '#'
        while (std::getline(infile, line))
        {
                if ( line.at(0) == '#' ) {
                    continue;
                }

            // reading in with loop
            for(uint32_t ix = 0; ix < nx; ix++)
            {
                for(uint32_t iy = 0; iy < ny; iy++)
                {
                    for(uint32_t iz = 0; iz < nz; iz++)
                    {
                        std::istringstream iss(line);
                        double ddx, ddy, ddz, mx, my, mz;
                        if (!(iss >> ddx >> ddy >> ddz >> mx >> my >> mz)){printf("Error while reading line!\n"); break;}
                        read_a(ix, iy, iz, 0) = mx;
                        read_a(ix, iy, iz, 1) = my;
                        read_a(ix, iy, iz, 2) = mz;
                        std::getline(infile, line);
                    }
                }
            }
            a = read_a;
        }
    }
    else
    {
        printf("Could not read file!");
    }
}


} // namespace magnumafcpp
