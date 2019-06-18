#include <gtest/gtest.h>
#include "../../../src/integrators/adaptive_runge_kutta.cpp"
#include "../../../src/integrators/controller.cpp"
#include "../../../src/state.cpp"
#include "../../../src/mesh.cpp"
#include "../../../src/func.cpp"
#include "../../../src/misc.cpp"
#include "../../../src/vtk_IO.cpp"

class RK : public AdaptiveRungeKutta{
    public:
        RK(std::string scheme, Controller controller) : AdaptiveRungeKutta(scheme, controller, false){};
    private:
        af::array f(const State& state){return state.t*sqrt(state.m);}
};


double analytic_result(double time){
    return 1./16. * pow(pow(time, 2) + 4, 2);
}

TEST(AdaptiveRungeKutta, BS23IntegrationTest) {
    RK rk("BS23", Controller(1e-15, 1e15, 1e-10, 1e-10));
    af::array m = af::constant(1.0, 1, f64);
    State state(Mesh(0, 0, 0, 0, 0, 0), Material(), m);
    for (int i=0; i<200; i++){
         rk.step(state);
         ASSERT_NEAR(afvalue(state.m), analytic_result(state.t), 1e-12 );
    }
}


TEST(AdaptiveRungeKutta, BS45IntegrationTest) {
    RK callback("BS45", Controller(1e-15, 1e15, 1e-10, 1e-10));
    af::array m = af::constant(1.0, 1, f64);
    State state(Mesh(0, 0, 0, 0, 0, 0), Material(), m);
    for (int i=0; i<100; i++){
         callback.step(state);
         ASSERT_NEAR(afvalue(state.m), analytic_result(state.t), 1e-8 );
    }
}


TEST(AdaptiveRungeKutta, DP45IntegrationTest) {
    RK callback("DP45", Controller(1e-15, 1e15, 1e-10, 1e-10));
    af::array m = af::constant(1.0, 1, f64);
    State state(Mesh(0, 0, 0, 0, 0, 0), Material(), m);
    for (int i=0; i<100; i++){
         callback.step(state);
         ASSERT_NEAR(afvalue(state.m), analytic_result(state.t), 1e-8 );
    }
}

TEST(AdaptiveRungeKutta, RKF45IntegrationTest) {
    RK callback("RKF45", Controller(1e-15, 1e15, 1e-10, 1e-10));
    af::array m = af::constant(1.0, 1, f64);
    State state(Mesh(0, 0, 0, 0, 0, 0), Material(), m);
    for (int i=0; i<100; i++){
         callback.step(state);
         ASSERT_NEAR(afvalue(state.m), analytic_result(state.t), 1e-8 );
    }
}


TEST(AdaptiveRungeKutta, DP78IntegrationTest) {
    RK callback("DP78", Controller(1e-15, 1e15, 1e-14, 1e-14));
    af::array m = af::constant(1.0, 1, f64);
    State state(Mesh(0, 0, 0, 0, 0, 0), Material(), m);
    for (int i=0; i<100; i++){
         callback.step(state);
         ASSERT_NEAR(afvalue(state.m), analytic_result(state.t), 1e-8 );
         //std::cout << "i=" << i << ", state.m= " << afvalue(state.m) << "; ";
    }
}


int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
