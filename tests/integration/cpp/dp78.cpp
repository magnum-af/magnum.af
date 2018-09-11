#include <gtest/gtest.h>
#include "../../../src/integrators/adaptive_runge_kutta.cpp"
#include "../../../src/integrators/controller.cpp"
#include "../../../src/state.cpp"
#include "../../../src/mesh.cpp"
#include "../../../src/func.cpp"
#include "../../../src/vtk_IO.cpp"
 
class Callback : public AdaptiveRungeKutta{
    public:
        Callback(std::string scheme = "DP78", Controller controller = Controller());
    private:
        af::array f(const State& state){
            return state.t*sqrt(state.m);
        }
};

Callback::Callback(std::string scheme, Controller controller) : AdaptiveRungeKutta(scheme, controller, false) {
};

TEST(DP78IntegrationTest, n) {
    Callback callback = Callback("DP78",Controller(1e-15, 1e15, 1e-14, 1e-14));
    array m = constant(1.0,1,f64);
    State state(Mesh(0,0,0,0,0,0), Param(), m);
    for (int i=0; i<100; i++){
         callback.step(state);
         std::cout << "i=" << i << ", state.m= " << afvalue(state.m) << "; ";
         ASSERT_NEAR(afvalue(state.m), 1./16. * pow(pow(state.t,2)+4,2), 1e-8 );
    }
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
