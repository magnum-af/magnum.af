#include "ascii_io.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>

namespace magnumafcpp
{


uint32_t colum_major_stride(const uint32_t i, const uint32_t j, const uint32_t k, const uint32_t l, const uint32_t ni, const uint32_t nj, const uint32_t nk) {
     return i + ni * (j + nj * (k + nk * l));
 }


void write_ascii(const af::array& a, const Mesh& mesh, std::string filename, bool verbose, int precision)
{
    const uint32_t nx = a.dims(0);
    const uint32_t ny = a.dims(1);
    const uint32_t nz = a.dims(2);
    const uint32_t n_scalar = a.dims(3);

    std::ofstream stream(filename);
    if( stream.is_open())
    {
        af::timer timer = af::timer::start();
        stream.precision(precision);
        int val_setw = precision + 4; // formatted at least for 0.*** values
        if (verbose) {printf("writing to ascii file %s with precision of %d\n", filename.c_str(), precision);}
        stream << "# " << a.dims(0) << " " <<  a.dims(1) << " " << a.dims(2) << " " << a.dims(3) << " " << mesh.dx << " " << mesh.dy << " " << mesh.dz << std::endl;
        stream << "# ascii encoded vector field" << std::endl;
        stream << "# plot e.g. with gnuplot: set view equal xyz; splot '" << filename << "' using 1:2:3:4:5:6 with vectors" << std::endl;
        stream << "#                       : set view equal xyz; splot '" << filename << "' using ($1*1e9):($2*1e9):($3*1e9):4:5:6 with vectors" << std::endl;
        stream << "#                       : set view equal xyz; splot '" << filename << "' using 1:2:3:4:5:6:(sqrt($4**2+$5**2+$6**2)) with vectors linecolor palette z" << std::endl;
        stream << "#                       : set view equal xyz; splot '" << filename << "' using ($1*1e9):($2*1e9):($3*1e9):4:5:6:4 with vectors linecolor palette" << std::endl;
        stream << "#\n# ";
        stream << std::setw(val_setw) << std::left << "x[m]" << "\t";
        stream << std::setw(val_setw) << std::left << "y[m]" << "\t";
        stream << std::setw(val_setw) << std::left << "z[m]" << "\t";
        stream << std::setw(val_setw) << std::left << "mx[a.u.]" << "\t";
        stream << std::setw(val_setw) << std::left << "my[a.u.]" << "\t";
        stream << std::setw(val_setw) << std::left << "mz[a.u.]" << "\n";

        // copying raw data to host is factor ~10 faster than accessing af::array with indices
        double *host_a = a.host<double>();

        for(uint32_t ix = 0; ix < nx; ix++)
        {
            for(uint32_t iy = 0; iy < ny; iy++)
            {
                for(uint32_t iz = 0; iz < nz; iz++)
                {
                    stream << std::setw(val_setw) << std::left << ix * mesh.dx << "\t";
                    stream << std::setw(val_setw) << std::left << iy * mesh.dy << "\t";
                    stream << std::setw(val_setw) << std::left << iz * mesh.dz << "\t";

                    for(uint32_t i_scalar = 0; i_scalar < n_scalar; i_scalar++)
                    {
                        stream << std::setw(val_setw) << std::left << host_a[colum_major_stride(ix, iy, iz, i_scalar, nx, ny, nz)];
                        if (i_scalar != (n_scalar - 1))
                        {
                            stream << "\t";
                        }
                    }
                    stream << "\n";
                    //stream << std::setw(val_setw) << std::left << host_a[colum_major_stride(ix, iy, iz, 0, nx, ny, nz)] << "\t";
                    //stream << std::setw(val_setw) << std::left << host_a[colum_major_stride(ix, iy, iz, 1, nx, ny, nz)] << "\t";
                    //stream << std::setw(val_setw) << std::left << host_a[colum_major_stride(ix, iy, iz, 2, nx, ny, nz)] << "\n";

                    //stream << std::setw(val_setw) << std::left << a(ix, iy, iz, 0).scalar<double>() << "\t";
                    //stream << std::setw(val_setw) << std::left << a(ix, iy, iz, 1).scalar<double>() << "\t";
                    //stream << std::setw(val_setw) << std::left << a(ix, iy, iz, 2).scalar<double>() << "\n";
                }
            }
        }
        af::freeHost(host_a);
        if (verbose){ printf("Wrote ascii file in %f [s]\n", timer.stop()); }
    }

}
void read_ascii(af::array& a, Mesh& mesh, std::string filename, bool verbose)
{
    std::ifstream infile(filename);
    if( infile.is_open())
    {
        af::timer timer = af::timer::start();
        if (verbose){
            printf("reading from ascii file %s\n", filename.c_str());
        }

        std::string line;

        // Getting mesh size and discretization
        uint32_t nx, ny, nz, n_scalar;
        double dx, dy, dz;
        std::getline(infile, line);
        line.erase(0, 2);
        std::istringstream iss1(line);
        if (!(iss1 >> nx >> ny >> nz >> n_scalar >> dx >> dy >> dz)){printf("Error while reading line!\n");}
        Mesh read_mesh(nx, ny, nz, dx, dy, dz);
        mesh = read_mesh;
        double* raw_read_a = new double[nx * ny * nz * n_scalar];
        //af::array read_a = af::constant(0, nx, ny, nz, 3, f64);

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
                        double ddx, ddy, ddz;
                        if (!(iss >> ddx >> ddy >> ddz)){printf("Error while reading line!\n"); break;}
                        for(uint32_t i_scalar = 0; i_scalar < n_scalar; i_scalar++)
                        {
                            double value;
                            if (!(iss >> value)){printf("Error while reading line!\n"); break;}
                            raw_read_a[colum_major_stride(ix, iy, iz, i_scalar, nx, ny, nz)] = value;
                            //if (!(iss >> ddx >> ddy >> ddz >> mx >> my >> mz)){printf("Error while reading line!\n"); break;}
                            //raw_read_a[colum_major_stride(ix, iy, iz, 0, nx, ny, nz)] = mx;
                            //raw_read_a[colum_major_stride(ix, iy, iz, 1, nx, ny, nz)] = my;
                            //raw_read_a[colum_major_stride(ix, iy, iz, 2, nx, ny, nz)] = mz;
                        }

                        //read_a(ix, iy, iz, 0) = mx;
                        //read_a(ix, iy, iz, 1) = my;
                        //read_a(ix, iy, iz, 2) = mz;
                        std::getline(infile, line);
                    }
                }
            }
            af::array read_a = af::array(nx, ny, nz, n_scalar, raw_read_a);
            delete[] raw_read_a;
            raw_read_a = NULL;
            a = read_a;
            //a = read_a;
        }
        if (verbose){ printf("Wrote ascii file in %f [s]\n", timer.stop()); }
    }
    else
    {
        printf("Could not read file!\n");
    }
}


} // namespace magnumafcpp
