#ifndef ZERO_CROSSING_H
#define ZERO_CROSSING_H

#include <array>
#include <iostream>
#include <functional>
#include <math.h>

/// Class calculating the zero crossing of a monotonic function f(x) up to a given precision.

/// Returns x and f(x) at the closest point to the zero crossing.
class ZeroCrossing{
    public:
        ZeroCrossing(std::function<double (double)> f, double precision = 1e-6, int max_runs = 5, double ix_min = 0, double ix_max = 1, int ix_n = 100, bool verbose = true);
        std::pair<double, double> calc_x_and_f(); ///< Calculates x and f(x)
    private:
        std::function<double (double)> f;
        const double precision;
        const int max_runs;
        double ix_min; // current x range min
        double ix_max; // current x range max
        const int ix_n; // x range discretization
        const bool verbose;

        std::array<double, 4> run_loop();
};

#endif
