#include <gtest/gtest.h>
#include "../../../src/integrators/adaptive_runge_kutta.cpp"
#include "../../../src/integrators/controller.cpp"
#include "../../../src/func.cpp"
 
af::array f(const double t, const af::array& m){ 
    return t*sqrt(m);
}

TEST(RKF45IntegrationTest, n) {
    AdaptiveRungeKutta rkf(&f, "RKF45", Controller(1e-15, 1e4, 1e-10, 1e-10));
    double t=0;
    array m = constant(1.0,1,f64);
    for (int i=0; i<100; i++){
         rkf.step(m,t);
         ASSERT_NEAR(afvalue(m), 1./16. * pow(pow(t,2)+4,2), 1e-8 );
    }
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
