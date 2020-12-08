#pragma once
#include "arrayfire.h"
#include "magnumafConfig_git.hpp"
#include <algorithm>
#include <boost/program_options.hpp>
#include <filesystem>
#include <iostream>
#include <string>

namespace magnumafcpp{


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
    desc.add_options()("help,h", "Produce help message");
    desc.add_options()("outdir,o", po::value<fs::path>(&outdir)->default_value("output_" + filename),
                       "Output directory. Will be created and is accessible via ArgParser.outdir(). Defaults to "
                       "'output_<binaryname>'.");
    desc.add_options()("no-overwrite,n", "Abort if outdir already exists. Prevents file overwriting.");
    desc.add_options()("gpu,g", po::value<unsigned>(), "set GPU number");
    desc.add_options()("verbose,v", "print parsing steps");
    desc.add_options()("posargs", po::value<std::vector<std::string>>(&posargs),
                       "positional arguments, access via ArgParser.posargs()");

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

        if (vm.count("gpu")) {
            if (vm.count("verbose")) {
                std::cout << "Setting GPU to " << vm["gpu"].as<unsigned>() << ".\n";
            }
            try {
                af::setDevice(vm["gpu"].as<unsigned>());
            } catch (const std::exception& e) {
                std::cerr << "af::setDevice(" << vm["gpu"].as<unsigned>() << ") caused: " << e.what() << std::endl;
                std::cerr << "GPU number not valid, skipping af::setDevice()." << std::endl;
            }
        }

        if (vm.count("outdir")) {
            if (vm.count("no-overwrite") and fs::exists(fs::path(outdir))) {
                std::cout << "active flag -n (--no-overwrite) prevents writing into already existing outdir " << outdir
                          << std::endl;
                std::cout << "Aborting..." << std::endl;
                std::exit(1);
            }

            bool mkdirs = fs::create_directories(outdir);

            if (vm.count("verbose") and mkdirs) {
                std::cout << "Created outdir " << fs::absolute(outdir) << std::endl;
            } else if (vm.count("verbose")) {
                std::cout << "Outputdir " << fs::absolute(outdir) << " already exists, files could be overwritten."
                          << std::endl;
            }
        }

        if (vm.count("posargs")) {
            if (vm.count("verbose")) {
                std::cout << "positional arguments are: " << vm["posargs"].as<std::vector<std::string>>() << std::endl;
            }
        } else {
            if (vm.count("verbose")) {
                std::cout << "No positional arguments, ArgParser.posargs() is empty." << std::endl;
            }
        }

        // Writing logging info
        const auto logfile = outdir / "log.txt";
        if (vm.count("verbose")) {
            std::cout << "Logging posargs info for script reproducability to" << fs::absolute(logfile) << std::endl;
        }
        std::ofstream ofst(logfile);
        for (const auto& elem : posargs) {
            const auto i = &elem - &posargs[0];
            ofst << "posargs[" << i << "] = " << elem << '\n';
        }
        ofst << "GIT_HEAD_SHA1=" << GIT_HEAD_SHA1 << std::endl;

        if (vm.count("verbose")) {
            af::info();
            std::cout << "Argument parsing finished." << std::endl;
        }

    } catch (const std::exception& e) {
        std::cerr << "Exception while parsing command line arguments: " << e.what() << '\n';
        std::cerr << desc << '\n';
        std::cerr << "Aborting..." << '\n';
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
