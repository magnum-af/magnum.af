#pragma once
#include "arrayfire.h"
#include "magnumafConfig_git.hpp"
#include <algorithm>
#include <boost/program_options.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

namespace magnumaf{


template <class T> std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
    copy(std::begin(v), std::end(v), std::ostream_iterator<T>(std::cout, " "));
    return os;
}

namespace argparser {

inline std::pair<std::filesystem::path, std::vector<std::string>> setup_and_get_posargs(int argc, char** argv) {
    namespace po = boost::program_options;
    namespace fs = std::filesystem;
    std::vector<std::string> posargs{};
    fs::path outdir{};

    const auto filename = fs::path(argv[0]).filename().string();
    po::options_description desc{"magnum.af ArgParser\nUsage: " + filename + " [options]"};
    desc.add_options()("help,h", "Produce help message.");
    desc.add_options()("outdir,o", po::value<fs::path>(&outdir)->default_value("output_" + filename),
                       "Output directory. Will be created and is accessible via ArgParser.outdir(). Defaults to "
                       "'output_<binaryname>'.");
    desc.add_options()("no-overwrite,n", "Abort if outdir already exists. Prevents file overwriting.");
    desc.add_options()("backend,b", po::value<std::string>(),
                       "'cuda', 'opencl' or 'cpu'. Select arrayfire backend via 'af::setBackend(b)'.");
    desc.add_options()("device,d", po::value<unsigned>(), "Set 'af::setDevice(d)', e.g. used for selecting a GPU.");
    desc.add_options()("verbose,v", "Make this parser verbose, printing parsing steps.");
    desc.add_options()("posargs", po::value<std::vector<std::string>>(&posargs),
                       "Positional arguments, access via 'ArgParser.posargs()'.");

    po::positional_options_description p;
    p.add("posargs", -1);

    try {
        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
        po::notify(vm);

        if (vm.count("help")) {
            std::cout << desc << "\n";
            std::exit(0);
        }

        if (vm.count("backend")) {
            const auto backend = vm["backend"].as<std::string>();
            if (vm.count("verbose")) {
                std::cout << "--backend: Setting af::setBackend() to " << backend << ".\n";
            }

            if (backend == "cuda" or backend == "opencl" or backend == "cpu") {
                try {
                    if (backend == "cuda") {
                        af::setBackend(AF_BACKEND_CUDA);
                    }
                    if (backend == "opencl") {
                        af::setBackend(AF_BACKEND_OPENCL);
                    }
                    if (backend == "cpu") {
                        af::setBackend(AF_BACKEND_CPU);
                    }

                } catch (const std::exception& e) {
                    std::cerr << "--backend: Error while calling af::setBackend() with '" << backend << "'.\n"
                              << "Maybe your hardware does not support the selected backend or drivers are not "
                                 "installed.\n";
                    throw;
                }
            } else {
                throw po::validation_error(po::validation_error::invalid_option_value, "backend");
            }
        }

        if (vm.count("device")) {
            if (vm.count("verbose")) {
                std::cout << "--device: Setting af::setDevice(" << vm["device"].as<unsigned>() << ").\n";
            }
            try {
                af::setDevice(vm["device"].as<unsigned>());
            } catch (const std::exception& e) {
                std::cerr << "--device: Error while calling af::setDevice(" << vm["device"].as<unsigned>()
                          << "). Device number not valid.\n";
                throw;
            }
        }

        // Location dependent, don't move after --outdir option.
        if (vm.count("no-overwrite")) {
            if (vm.count("verbose")) {
                std::cout << "--no-overwrite: Checking if --outdir " << fs::absolute(outdir) << "  exists.\n";
            }
            if (fs::exists(fs::path(outdir))) {
                throw std::runtime_error(
                    "--no-overwrite: Error, active flag '--no-overwrite' prevents writing into already "
                    "existing outdir" +
                    fs::absolute(outdir).string());
            }
        }

        if (vm.count("outdir")) {
            bool mkdirs = fs::create_directories(outdir);
            if (vm.count("verbose") and mkdirs) {
                std::cout << "--outdir: Created outdir " << fs::absolute(outdir) << '\n';
            } else if (vm.count("verbose")) {
                std::cout << "--outdir: Outputdir " << fs::absolute(outdir)
                          << " already exists, files could be overwritten.\n";
            }
        }

        if (vm.count("posargs")) {
            if (vm.count("verbose")) {
                std::cout << "--posargs: Positional arguments are: " << vm["posargs"].as<std::vector<std::string>>()
                          << '\n';
            }
        } else {
            if (vm.count("verbose")) {
                std::cout << "--posargs: No positional arguments, ArgParser.posargs() is empty.\n";
            }
        }

        // Writing logging info
        const auto logfile = outdir / "log.txt";
        if (vm.count("verbose")) {
            std::cout << "Logging posargs info for script reproducibility to" << fs::absolute(logfile) << '\n';
        }
        std::ofstream ofst(logfile);
        for (const auto& elem : posargs) {
            const auto i = &elem - &posargs[0];
            ofst << "posargs[" << i << "] = " << elem << '\n';
        }
        ofst << "GIT_HEAD_SHA1=" << GIT_HEAD_SHA1 << '\n';

        if (vm.count("verbose")) {
            af::info();
            std::cout << "Argument parsing finished.\n";
        }

    } catch (const std::exception& e) {
        std::cerr << "Exception occurred while parsing command line arguments:\n"
                  << e.what() << '\n'
                  << "Showing help (" << filename << " -h):\n"
                  << desc << '\n'
                  << "Aborting...\n";
        std::exit(1);
    }

    return {outdir, posargs};
}
} // namespace argparser

/// Parser for command line options.
class ArgParser {
  public:
    ArgParser(int argc, char** argv) : outdir_posargs(argparser::setup_and_get_posargs(argc, argv)) {}
    std::vector<std::string> posargs() const {
        return outdir_posargs.second;
    } ///< Positional arguments for usage after parsing e.g. simulation parameters
    std::filesystem::path outdir() const { return outdir_posargs.first; }            ///< Output directory.
    const std::pair<std::filesystem::path, std::vector<std::string>> outdir_posargs; ///< Returns outdir and posargs.

  private:
};
}
