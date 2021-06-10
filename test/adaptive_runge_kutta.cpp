#include "integrators/adaptive_runge_kutta.hpp"
#include "util/util.hpp"
#include <gtest/gtest.h>

using namespace magnumafcpp;

class RK : public AdaptiveRungeKutta {
  public:
    RK(std::string scheme, Controller controller) : AdaptiveRungeKutta(scheme, controller, false){};

  private:
    af::array f(const State& state) const { return state.t * sqrt(state.m); }
};

double analytic_result(double time) { return 1. / 16. * pow(pow(time, 2) + 4, 2); }

TEST(AdaptiveRungeKutta, BS23IntegrationTest) {
    RK rk("BS23", Controller(1e-15, 1e15, 1e-10, 1e-10));
    af::array m = af::constant(0.0, 1, 1, 1, 3, f64);
    m(0, 0, 0, 0) = 1;
    State state(Mesh(0, 0, 0, 0, 0, 0), 1, m);
    for (int i = 0; i < 200; i++) {
        rk.step(state);
        EXPECT_NEAR(state.m.scalar<double>(), analytic_result(state.t), 1e-12);
    }
}

TEST(AdaptiveRungeKutta, BS45IntegrationTest) {
    RK callback("BS45", Controller(1e-15, 1e15, 1e-10, 1e-10));
    af::array m = af::constant(0.0, 1, 1, 1, 3, f64);
    m(0, 0, 0, 0) = 1;
    State state(Mesh(0, 0, 0, 0, 0, 0), 1, m);
    for (int i = 0; i < 100; i++) {
        callback.step(state);
        EXPECT_NEAR(state.m.scalar<double>(), analytic_result(state.t), 1e-8);
    }
}

TEST(AdaptiveRungeKutta, DP45IntegrationTest) {
    RK callback("DP45", Controller(1e-15, 1e15, 1e-10, 1e-10));
    af::array m = af::constant(0.0, 1, 1, 1, 3, f64);
    m(0, 0, 0, 0) = 1;
    State state(Mesh(0, 0, 0, 0, 0, 0), 1, m);
    for (int i = 0; i < 100; i++) {
        callback.step(state);
        EXPECT_NEAR(state.m.scalar<double>(), analytic_result(state.t), 1e-8);
    }
}

TEST(AdaptiveRungeKutta, RKF45IntegrationTest) {
    RK callback("RKF45", Controller(1e-15, 1e15, 1e-10, 1e-10));
    af::array m = af::constant(0.0, 1, 1, 1, 3, f64);
    m(0, 0, 0, 0) = 1;
    State state(Mesh(0, 0, 0, 0, 0, 0), 1, m);
    for (int i = 0; i < 100; i++) {
        callback.step(state);
        EXPECT_NEAR(state.m.scalar<double>(), analytic_result(state.t), 1e-8);
        // std::cout << "i=" << i << ", t= " << state.t << ", m=" << state.m.scalar<double>() << std::endl;
    }
}

TEST(AdaptiveRungeKutta, DP78IntegrationTest) {
    RK callback("DP78", Controller(1e-15, 1e15, 1e-14, 1e-14));
    af::array m = af::constant(0.0, 1, 1, 1, 3, f64);
    m(0, 0, 0, 0) = 1;
    State state(Mesh(0, 0, 0, 0, 0, 0), 1, m);
    for (int i = 0; i < 100; i++) {
        callback.step(state);
        EXPECT_NEAR(state.m.scalar<double>(), analytic_result(state.t), 1e-8);
        // std::cout << "i=" << i << ", state.m= " << state.m.scalar<double>() << "; ";
    }
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
