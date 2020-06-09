#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace magnumafcpp;

template <typename T> class AdaptiveRK {
  public:
    void step(double& t, T&);
    virtual T f(const double t, const T&) = 0;
    virtual ~AdaptiveRK(){};
};

template <typename T> class TempSquareFormula : public AdaptiveRK<T> {
  public:
    T f(const double t, const T& m) { return t * sqrt(m); }
    TempSquareFormula(){};
};

class SquareFormula : public AdaptiveRK<af::array> {
  public:
    af::array f(const double t, const af::array& m) { return t * af::sqrt(m); }
};

class TestState {
  public:
    af::array m{};
    double t{0};
    TestState(af::array m) : m(m) {}
};

class EquationState : public AdaptiveRK<TestState> {
  public:
    TestState heff(TestState state) { return state; }
    TestState f(const double t, const TestState& m) {
        auto temp = m;
        temp.m = t * af::sqrt(m.m) + heff(m).m;
        return temp;
    }
};

int main(int argc, char** argv) {
    for (int i = 0; i < argc; i++) {
        std::cout << "Parameter " << i << " was " << argv[i] << std::endl;
    }
    std::string filepath(argc > 1 ? argv[1] : "output_magnum.af/");
    af::setDevice(argc > 2 ? std::stoi(argv[2]) : 0);
    af::info();
    // auto squaref = SquareFormula();
    double t = 0;
    SquareFormula squaref{};
    af::array m = af::constant(1, 1, f64);
    af::print("", squaref.f(t, m));
    TempSquareFormula<af::array> afsquare{};
    af::print("", afsquare.f(t, m));
    TestState teststate(m);
    EquationState equationstate{};
    af::print("", equationstate.f(t, m).m);
    return 0;
}
