#pragma once
#include "arrayfire.h"
#include <iostream>
#include <ostream>

namespace magnumafcpp {

class StageTimer {
  public:
    double get_time_and_restart_timer() {
        if (sync) {
            af::sync();
        }
        double time = af::timer::stop(timer);
        accumulated_time += time;
        timer = af::timer::start();
        // timer.start();
        stage++;
        return time;
    }
    void print_stage(std::string stagename = "", std::ostream& stream = std::cout) {
        stream << "timing ";
        stagename != "" ? stream << "'" << stagename << "'" : stream << "stage " << stage;
        stream << ": " << get_time_and_restart_timer() << " [s]" << std::endl;
    }
    void print_accumulated(std::ostream& stream = std::cout) {
        stream << "timing accumulated: " << accumulated_time << " [s]" << std::endl;
    }

  private:
    af::timer timer{af::timer::start()};
    double accumulated_time{0};
    int stage{0};
    bool sync{true};
};

} // namespace magnumafcpp
