#pragma once
#include "arrayfire.h"
#include <ostream>
#include <iostream>

namespace magnumaf{

class StageTimer{
  public:
    double get_time_and_restart(){
        if (sync) af::sync();
        double time = af::timer::stop(timer);
        accumulated_time += time;
        timer = af::timer::start();
        //timer.start();
        stage++;
        return time;
    }
    void print_stage(std::ostream& stream, std::string stagename = ""){
        stream << "timing ";
        stagename != "" ? stream << "'" << stagename << "'" : stream << "stage " << stage;
        stream << ": " << get_time_and_restart() << " [s]"<< std::endl;
    }
    void print_accumulated(std::ostream& stream){
        stream << "timing accumulated: " << accumulated_time << " [s]"<< std::endl;
    }

  private:
    af::timer timer{af::timer::start()};
    int stage{0};
    bool sync {true};
    double accumulated_time{0};
};

}// namespace magnumaf
