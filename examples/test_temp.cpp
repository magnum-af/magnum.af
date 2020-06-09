#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace magnumafcpp;

template <typename Input, typename Return> class AdaptiveRK {
  public:
    void step(double& t, Input&);
    virtual Return f(const double t, const Input&) = 0;
    virtual ~AdaptiveRK(){};
};

template <typename T> class TempSquareFormula : public AdaptiveRK<T, T> {
  public:
    T f(const double t, const T& m) { return t * sqrt(m); }
    TempSquareFormula(){};
};

class SquareFormula : public AdaptiveRK<af::array, af::array> {
  public:
    af::array f(const double t, const af::array& m) { return t * af::sqrt(m); }
};

class TestState {
  public:
    af::array m{};
    double t{0};
    TestState(af::array m) : m(m) {}
};

class EquationState : public AdaptiveRK<TestState, TestState> {
  public:
    // TestState heff(TestState state) { return state; }
    af::array heff(TestState state) { return state.m; }
    TestState f(const double t, const TestState& m) {
        auto temp = m;
        temp.m = t * af::sqrt(m.m) + heff(m);
        return temp;
    }
};

class EquationStateArray : public AdaptiveRK<TestState, af::array> {
  public:
    af::array heff(TestState state) { return state.m; }
    af::array f(const double t, const TestState& m) {
        return t * af::sqrt(m.m) + heff(m);
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
    EquationStateArray equationstatearray{};
    af::print("", equationstatearray.f(t, m));
    return 0;
}
