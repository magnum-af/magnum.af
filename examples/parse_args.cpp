#include "magnum_af.hpp"
#include "util/arg_parser.hpp"

using namespace magnumafcpp;

int main(int argc, char** argv) {
    const auto [outdir, posargs] = ArgParser(argc, argv).outdir_posargs;

    std::cout << "outdir = " << std::filesystem::absolute(outdir) << std::endl;

    for (const auto& elem : posargs) {
        const auto i = &elem - &posargs[0];
        std::cout << "posargs[" << i << "]: " << elem << std::endl;
    }

    //// defining values from positional arguments while providing default values
    // direct
    const double val0{posargs.size() > 0 ? std::stod(posargs[0]) : 1.0};
    const std::size_t val1{posargs.size() > 1 ? std::stoul(posargs[1]) : 42};
    const int val2{posargs.size() > 2 ? std::stoi(posargs[2]) : -5};

    // via lamda
    auto stod = [&posargs](std::size_t i, auto default_value) {
        return posargs.size() > i ? std::stod(posargs[i]) : default_value;
    };
    const double val3{stod(4, -4.)};
    const double val4{stod(5, -5.)};
    const double val5{stod(6, -6.)};

    std::cout << "val0: " << val0 << std::endl;
    std::cout << "val1: " << val1 << std::endl;
    std::cout << "val2: " << val2 << std::endl;
    std::cout << "val3: " << val3 << std::endl;
    std::cout << "val4: " << val4 << std::endl;
    std::cout << "val5: " << val5 << std::endl;
    return 0;
}
