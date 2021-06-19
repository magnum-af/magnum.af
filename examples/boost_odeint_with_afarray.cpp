#include "arrayfire.h"
#include <boost/numeric/odeint.hpp>
#include <iostream>

namespace boost::numeric::odeint {
// specialization for af::array
// Adapted form boost/numeric/odeint/external/vexcl/vexcl_norm_inf.hpp
// This works:
template <> struct vector_space_norm_inf<af::array> {
    typedef double result_type; // Needed for internals
    result_type operator()(const af::array& x) const {
        return af::max(af::max(af::max(af::max(af::abs(x), 0), 1), 2), 3).as(f64).scalar<double>();
    }
};

} // namespace boost::numeric::odeint

double analytic_result(double time) { return 1. / 16. * std::pow(std::pow(time, 2) + 4, 2); }

int main() {
    std::cout << "Start" << std::endl;
    namespace ode = boost::numeric::odeint;
    auto boost_callback = [](const af::array& y0, af::array& dxdt, const double t) { dxdt = t * af::sqrt(y0); };

    struct cout_m_and_time {
        void operator()(const af::array& y, double t) {
            const auto y_ = y.scalar<double>();
            const auto analyt = analytic_result(t);
            std::cout << t << "\t" << y_ << "\t" << analyt << "\t" << std::abs(y_ - analyt) << std::endl;
        };
    };

    using stepper_type = ode::runge_kutta_dopri5<af::array, double, af::array, double, ode::vector_space_algebra>;
    // typedef ode::runge_kutta_cash_karp54_classic<af::array, double, af::array, double, ode::vector_space_algebra>
    // stepper_type; typedef ode::runge_kutta_fehlberg78<af::array, double, af::array, double,
    // ode::vector_space_algebra> stepper_type;

    const double eps_abs = 1e-6;
    const double eps_rel = 1e-6;
    const double start_time = 0;
    const double end_time = 1e2;
    const double dt = 1e-3;
    af::array y = af::constant(1.0, af::dim4(1, 1, 1, 1), f64);
    size_t steps = integrate_adaptive(make_controlled(eps_abs, eps_rel, stepper_type()), boost_callback, y, start_time,
                                      end_time, dt, cout_m_and_time{});
    std::cout << "Finished after steps=" << steps << std::endl;
    return 0;
}
