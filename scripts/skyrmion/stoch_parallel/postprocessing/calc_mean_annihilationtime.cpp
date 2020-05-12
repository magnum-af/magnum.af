// Usage:
// g++ -std=c++14 -o calc_mean_anihilation calc_mean_anihilationtime.cpp
//./calc_mean_anihilation $PWD/anihilationtime.dat $PWD
#include <cmath>
#include <fstream>
#include <iostream>
#include <list>
#include <sstream>
#include <string>
#include <vector>

using namespace magnumafcpp;

using namespace std;

int main(int argc, char** argv) {

    string indatapath(argv[1]);
    ifstream stream((indatapath).c_str());
    double detect_time;
    vector<double> vec_detect_time;
    double state_t;
    int i;
    int reverse;
    int ID;
    if (stream.is_open()) {
        while (stream >> detect_time >> state_t >> i >> reverse >> ID) {
            // cout << detect_time << "\t" << state_t << "\t" << i <<"\t" <<
            // reverse << "\t" << ID << std::endl;
            vec_detect_time.push_back(detect_time);
        }
        double mean_detect_time = 0;
        for (auto const& value : vec_detect_time) {
            mean_detect_time += value;
            // std::cout << value << std::endl;
        }
        mean_detect_time /= (double)vec_detect_time.size();

        double unbiased_sample_variance =
            0; // s^2= 1/(n-1) sum(y_i - y_mean)^2 from i = 1 to n
        for (double n : vec_detect_time) {
            unbiased_sample_variance += pow(n - mean_detect_time, 2);
        }
        unbiased_sample_variance /= ((double)vec_detect_time.size() - 1);
        double unbiased_sample_sigma = sqrt(unbiased_sample_variance);

        std::cout << "#mean_detect_time <<<< unbiased_sample_sigma <<<< "
                     "unbiased_sample_variance <<<< .size()"
                  << std::endl;
        std::cout << mean_detect_time << "\t" << unbiased_sample_sigma << "\t"
                  << unbiased_sample_variance << "\t" << vec_detect_time.size()
                  << "\n"
                  << std::endl;

        ofstream stream_out(
            (string(argv[2]) + "/mean_annihilationtime.dat").c_str(),
            std::ios_base::trunc);
        stream_out << mean_detect_time << "\t" << unbiased_sample_sigma << "\t"
                   << unbiased_sample_variance << "\t" << vec_detect_time.size()
                   << std::endl;
        stream_out << "#mean_detect_time <<<< unbiased_sample_sigma <<<< "
                      "unbiased_sample_variance <<<< vec_detect_time.size() "
                   << std::endl;
        stream_out.close();
    } else
        std::cout << "Could not open file" << std::endl;
    stream.close();

    return 0;
}
