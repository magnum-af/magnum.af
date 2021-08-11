#include <cmath>

#include "state.hpp"
#include "util/arg_parser.hpp"
#include "util/geometry.hpp"
#include <iomanip>

// from:
// https://stackoverflow.com/questions/669438/how-to-get-memory-usage-at-runtime-using-c
#include <fstream>
#include <ios>
#include <iostream>
#include <string>
#include <unistd.h>

//////////////////////////////////////////////////////////////////////////////
//
// process_mem_usage(double &, double &) - takes two doubles by reference,
// attempts to read the system-dependent data for a process' virtual memory
// size and resident set size, and return the results in KB.
//
// On failure, returns 0.0, 0.0

void process_mem_usage(double& vm_usage, double& resident_set) {
    using std::ifstream;
    using std::ios_base;
    using std::string;

    vm_usage = 0.0;
    resident_set = 0.0;

    // 'file' stat seems to give the most reliable results
    //
    ifstream stat_stream("/proc/self/stat", ios_base::in);

    // dummy vars for leading entries in stat that we don't care about
    //
    string pid, comm, state, ppid, pgrp, session, tty_nr;
    string tpgid, flags, minflt, cminflt, majflt, cmajflt;
    string utime, stime, cutime, cstime, priority, nice;
    string O, itrealvalue, starttime;

    // the two fields we want
    //
    unsigned long vsize = 0;
    long rss = 0;

    stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr >> tpgid >> flags >> minflt >> cminflt >>
        majflt >> cmajflt >> utime >> stime >> cutime >> cstime >> priority >> nice >> O >> itrealvalue >> starttime >>
        vsize >> rss; // don't care about the rest

    stat_stream.close();

    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    vm_usage = vsize / 1024.0;
    resident_set = rss * page_size_kb;
}

using namespace magnumaf;

int main(int argc, char** argv) {
    const auto [outdir, posargs] = ArgParser(argc, argv).outdir_posargs;
    const int n_writes(posargs.size() > 0 ? std::stoi(posargs[0]) : 200);
    std::cout << "n_writes = " << n_writes << std::endl;

    Mesh mesh(100, 25, 1, 1e-9, 2e-9, 3e-9);

    // Initial magnetic field
    State state(mesh, 8e5, util::init_sp4(mesh));
    for (int i = 0; i < n_writes; i++) {
        double vm = std::numeric_limits<double>::quiet_NaN(), rss = std::numeric_limits<double>::quiet_NaN();
        process_mem_usage(vm, rss);
        std::cout << std::left << "i=" << std::setw(6) << i << "[%]=" << std::setw(5)
                  << 100. * (static_cast<double>(i)) / (static_cast<double>(n_writes)) << "VIRT[GB]=" << std::setw(10) << vm / 1e6
                  << "RES[GB]=" << std::setw(9) << rss / 1e6 << "VIRT=" << std::setw(8) << vm << "RES=" << std::setw(6)
                  << rss << std::endl;
        state.write_vti(outdir / "minit");
    }
    return 0;
}
