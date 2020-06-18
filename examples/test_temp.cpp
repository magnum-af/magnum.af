#include "arrayfire.h"
#include "magnum_af.hpp"

using namespace magnumafcpp;

class StateInterface {
  public:
    af::array m;
    af::array Ms_field;
    double Ms;
    double t;
    // virtual int specific_type() = 0; // makes base class abstract
    virtual ~StateInterface() {} // Makes Interface Class abstract
};

class MyState : public StateInterface {
  public:
    int nz;
    // int specific_type() { return 2; }
};

template <typename Input, typename Return> class AdaptiveRK {
  public:
    // note: input is treated as const by a workaround
    Return FakeRKF(const double t, Input& input, const double dt) {
        // Input& should be const, we have to manipulate input.m, however
        // following is not working as abstract class cant be instantiated
        // Se we have to work on Input& and reset m aferwards to immitate const
        // Input& behaviour:
        // Input temp = input; auto temp = input; workaround:
        // save input.m and at the end reapply it
        const Return minit = input.m;
        // stage1
        Return k1 = dt * f(t, input);

        // stage2
        double t_current_step = t + 1. / 4. * dt;
        input.m = input.m + 1. / 4. * k1;
        Return k2 = dt * f(t_current_step, input);
        // stage3
        // TODO

        auto result = k1 + k2;
        input.m = minit; // Note: workaround to achive Input& as const
        // Preserving input.m is important if step is not accepted by
        // controller.
        return result;
    }
    void step(double& t, Input& input) {
        input.m = FakeRKF(t, input, 1);
        t += 1;
    }
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

class EquationStateInterface : public AdaptiveRK<StateInterface, af::array> {
  public:
    // TODO nz not accessible
    // af::array heff(const StateInterface& state) { return state.nz * state.m;}
    af::array heff(const StateInterface& state) { return state.m; }
    af::array f(const double t, const StateInterface& m) {
        return t * af::sqrt(m.m) + heff(m);
    }
};

class EquationHybridState : public AdaptiveRK<StateInterface, af::array> {
  public:
    // TODO nz not accessible
    af::array heff(const MyState& state) { return state.nz * state.m; }
    // TODO no const possible:
    af::array f(const double t, const StateInterface& m) {
        return t * af::sqrt(m.m) + heff(dynamic_cast<const MyState&>(m));
    }
};

class EquationMyState : public AdaptiveRK<MyState, af::array> {
  public:
    af::array heff(const MyState& state) { return state.nz * state.m; }
    af::array f(const double t, const MyState& m) {
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
    double t = 1;
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

    MyState mystate{};
    mystate.m = af::constant(1, 1, f64);

    // Passing MyState as StateInterface&
    EquationStateInterface equationstateinterface{};
    af::print("", equationstateinterface.f(t, mystate));
    equationstateinterface.step(t, mystate);
    af::print("", equationstateinterface.f(t, mystate));

    EquationHybridState hybridstate{};
    af::print("", hybridstate.f(t, mystate));
    hybridstate.step(t, mystate);
    af::print("", hybridstate.f(t, mystate));

    EquationMyState equationmystate{};
    af::print("", equationmystate.f(t, mystate));
    equationmystate.step(t, mystate);
    af::print("", equationmystate.f(t, mystate));
    return 0;
}
