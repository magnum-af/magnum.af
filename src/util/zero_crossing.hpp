#pragma once
#include <array>
#include <functional>

namespace magnumaf{

/// Class calculating the zero crossing of a monotonic function f(x) up to a given precision.

/// Returns x and f(x) at the closest point to the zero crossing.
/// In one run, f(x) is evaluated for x in [x_min, x_max] on ix_n + 1 equidistant points.
/// If a sign change in f(x) is encountert, the x-iteration is stopped and the sarch interval is updated to the two nearest points arount the sign change. This interval is used for the next run.
/// \param f callback function
/// \param precision sets precision which f(x) must fulfill
/// \param max_runs sets numer of runs
/// \param ix_min sets initial minimal x
/// \param ix_max sets initial maximal x
/// \param ix_n sets discretization of x interval
/// \param max_runs sets the numer of runs
/// \param verbose manages output levels:\n
///             false (==0): no output\n
///             true  (==1): one status output after execution\n
///             2:           also print values for every loop\n
///             3:           also print values for every evaluation of f(x)

class ZeroCrossing{
    public:
        ZeroCrossing(std::function<float (float)> f, float precision = 1e-6, int max_runs = 5, float ix_min = 0, float ix_max = 1, int ix_n = 100, int verbose = true);
        std::pair<float, float> calc_x_and_f(); ///< Calculates x and f(x)
    private:
        std::function<float (float)> f;
        const float precision;
        const int max_runs;
        float ix_min; // current x range min
        float ix_max; // current x range max
        const int ix_n; // x range discretization
        const int verbose;

        std::array<float, 4> run_loop();
};

}// namespace magnumaf
