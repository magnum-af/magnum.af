#include "constants.hpp"
#include "field_terms/micro/demag_field_pbc.hpp"
#include "util/arg_parser.hpp"
#include "util/prime_factors.hpp"
#include "util/vtk_io.hpp"
#include <chrono>

using namespace magnumaf;

int main(int argc, char** argv) {
    const auto [outdir, posargs] = ArgParser(argc, argv).outdir_posargs;
    af::info();

    auto run_for_mesh = [outdir](std::size_t N, af::dtype type) {
        DemagFieldPBC demag_pbc;
        Mesh mesh{N, N, N, 10e-9, 10e-9, 10e-9};

        const double Ms = 8e5;
        af::array m = af::constant(0, mesh::dims_v(mesh), type);
        m(af::span, af::span, af::span, 0) = 1;
        m(af::seq(3, 6), af::seq(3, 6), af::seq(3, 6), 0) = -1;
        State state(mesh, Ms, m);

        auto start0 = std::chrono::high_resolution_clock::now();
        const af::array H_in_Apm = demag_pbc.H_in_Apm(state);
        auto end0 = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double> elapsed0 = end0 - start0;
        std::cout << "Warmup in " << elapsed0.count() << std::endl;

        if (false) {
            try {
                auto [m_ref, mesh_ref] = vtr_reader(outdir / "m_ref.vtr");
                auto [h_ref, mesh_href] = vtr_reader(outdir / "h_ref.vtr");
                af::array diff_h = H_in_Apm(af::span, af::span, af::span, 0) - h_ref;
                auto diff_h_sum = af::sum(af::sum(af::sum(diff_h, 0), 1), 2);
                std::cout << "diff_h_sum    " << diff_h_sum.scalar<double>() << std::endl;
                diff_h = diff_h * diff_h;
                diff_h_sum = af::sum(af::sum(af::sum(diff_h, 0), 1), 2);
                std::cout << "diff_h_sum**2 " << diff_h_sum.scalar<double>() << std::endl;
            } catch (const std::exception& e) {
                std::cout << "Exception while comparing with reference solution:" << std::endl;
                std::cout << e.what() << std::endl;
            }
        }

        if (false) { // Switch VTK on/off
            vti_writer_micro(m, mesh, outdir / "m.vti");
            vti_writer_micro(H_in_Apm, mesh, outdir / "H_in_Apm.vti");
        }

        af::sync();
        std::ofstream stream(outdir / "timings.dat", std::ios::app);

        auto start = std::chrono::high_resolution_clock::now();
        const auto ii = 100;
        if (false) { // switch accumulate on/off
            auto accum_H_in_Apm = demag_pbc.H_in_Apm(state);
            af::eval(accum_H_in_Apm);
            for (int i = 0; i < ii; ++i) {
                // std::cout << "step " << i << std::endl;
                accum_H_in_Apm += demag_pbc.H_in_Apm(state); // gets stuck with 220x220x220, propably due to mem-cap
                af::eval(accum_H_in_Apm);
            }
        } else {
            for (int i = 0; i < ii; ++i) {
                auto tmp = demag_pbc.H_in_Apm(state);
                af::eval(tmp);
            }
        }
        af::sync();

        auto end = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double> elapsed = end - start;
        auto average_time = elapsed.count() / static_cast<double>(ii);
        const auto NNN = mesh.nx * mesh.ny * mesh.nz;
        std::cout << "Cells " << NNN << ", run " << ii << " H_demag_pbc evals in " << elapsed.count()
                  << " with average of " << elapsed.count() / static_cast<double>(ii) << std::endl;
        auto size_in_GB = static_cast<double>(NNN * sizeof(double)) / 1.e9;

        std::cout << "N= " << N << ", sizeof array in GB: " << size_in_GB << std::endl;
        stream << N << "\t" << std::setw(8) << NNN << "\t" << std::setw(10) << std::left << size_in_GB << "\t"
               << average_time << std::endl;
    };

    std::ofstream ofs(outdir / "timings.dat", std::ios::app);
    ofs << "#N\t"
        << "NNN\t"
        << "size[GB]\t"
        << "average_time[s]" << std::endl;
    for (std::size_t i = 10; i <= 260; i += 10) {
        if (af::getActiveBackend() == AF_BACKEND_OPENCL and util::max_of_prime_factors(i) > 13) {
            std::cout << "Skipping step as N=" << i << " dim not valid for OpenCL FFT" << std::endl;
        } else {
            run_for_mesh(i, f64);
        }
    }
}
