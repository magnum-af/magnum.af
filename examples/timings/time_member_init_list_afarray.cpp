#include "arrayfire.h"
#include <iostream>
#include <utility>


// struct naming: if callsite passes by lval: first term (e.g. CpCp is copy and copy)
// if callsite passes by rval (via std::move): second term (e.g. MvCp is move and copy)

struct CpCp_or_MvCp {
    af::array data;
    explicit CpCp_or_MvCp(af::array data) : data(data) {}
    // Backup (clang-tidy performance):
    // Reference: CpCp_or_MvCp(af::array data) : data(data) {}
};

struct CpMv_or_MvMv {
    af::array data;
    explicit CpMv_or_MvMv(af::array data) : data(std::move(data)) {}
};

struct Cp_or_Mv {
    af::array data;
    explicit Cp_or_Mv(af::array const& data) : data(data) {}
    explicit Cp_or_Mv(af::array&& data) : data(std::move(data)) {}
};

int main() {
    // int main(int argc, char** argv) {

    // const auto N = 1;
    // const auto dims = af::dim4(1, 1, 1, 1);
    const auto N = 6e5;
    const auto dims = af::dim4(1000, 1000, 10, 1);
    const auto type = af::dtype::f64;

    // Warmup
    {
        auto timer = af::timer::start();
        for (auto n = 0; n < N; ++n) {
            auto data = af::constant(42, dims, type);
        }
        double time_warm = af::timer::stop(timer);
        std::cout << "Warmup Time Ctor only= " << time_warm << " [s]\n";
    }

    // Ref timing
    double time_ctor;
    {
        auto timer = af::timer::start();
        for (auto n = 0; n < N; ++n) {
            auto data = af::constant(42, dims, type);
        }
        time_ctor = af::timer::stop(timer);
        std::cout << "Reference Time Ctor only = " << time_ctor << " [s]\n";
    }

    // lambda timings
    // using [](auto&& data){...} forwards arugments perfectly, (think of T&& data)
    // compared to explicit below: no additional Ctor calls which would affect timings
    if (true) {
        // Note: whether data is moved or not depends on callsite of lambda!
        // One lambda is sufficient for loop which either copies or move data
        const auto timeit = [N, dims, time_ctor](const auto& string, auto callback) {
            auto timer = af::timer::start();
            for (auto n = 0; n < N; ++n) {
                auto data = af::constant(42, dims, type);
                // auto data_lambda = callback(data);
                auto data_lambda = callback(std::move(data));
                // Note: using std::move() here only important if
                // callback requires rval ref of expicit type (e.g. double&&)
                // For this timing we could skip it without creating extra ops
            }
            double time = af::timer::stop(timer);
            std::cout << "Lambda: Time " << string << " = " << time << " [s], t-ref= " << time - time_ctor << " [s]\n";
            return time;
        };

        // Minimal impact lambdas, passed lambda forward-moves arg by rval ref:
        std::cout << std::endl;
        // timeit("Ctor only (warmup)", [](auto&&) { return 0; }); // causes var not used warning
        timeit("Ctor warm", [](auto&& data) { return data; });
        timeit("Ctor copy", [](auto&& data) { return data; });
        timeit("Ctor move", [](auto&& data) { return std::move(data); });
        timeit("Ctor Copy Copy", [](auto&& data) { return CpCp_or_MvCp{data}; });
        timeit("Ctor Move Copy", [](auto&& data) { return CpCp_or_MvCp{std::move(data)}; });
        timeit("Ctor Copy Move", [](auto&& data) { return CpMv_or_MvMv{data}; });
        timeit("Ctor Move Move", [](auto&& data) { return CpMv_or_MvMv{std::move(data)}; });
        timeit("Ctor Only Copy", [](auto&& data) { return Cp_or_Mv{data}; });
        timeit("Ctor Only Move", [](auto&& data) { return Cp_or_Mv{std::move(data)}; });
    }

    // Explicit Timings
    if (false) {
        {
            auto timer = af::timer::start();
            for (auto n = 0; n < N; ++n) {
                auto data = af::constant(42, dims, type);
                auto data_copy_copy = CpCp_or_MvCp(data);
            }
            double time = af::timer::stop(timer);
            // std::cout << "Time copy copy = " << time << " [s]\n";
            std::cout << "Time copy copy = " << time << " [s], t-ref= " << time - time_ctor << " [s]\n";
        }
        {
            auto timer = af::timer::start();
            for (auto n = 0; n < N; ++n) {
                auto data = af::constant(42, dims, type);
                auto data_move_copy = CpCp_or_MvCp(std::move(data));
            }
            double time = af::timer::stop(timer);
            std::cout << "Time move copy = " << time << " [s], t-ref= " << time - time_ctor << " [s]\n";
        }
        {
            auto timer = af::timer::start();
            for (auto n = 0; n < N; ++n) {
                auto data = af::constant(42, dims, type);
                auto data_copy_move = CpMv_or_MvMv(data);
            }
            double time = af::timer::stop(timer);
            std::cout << "Time copy move = " << time << " [s], t-ref= " << time - time_ctor << " [s]\n";
        }
        {
            auto timer = af::timer::start();
            for (auto n = 0; n < N; ++n) {
                auto data = af::constant(42, dims, type);
                auto data_move_move = CpMv_or_MvMv(std::move(data));
            }
            double time = af::timer::stop(timer);
            std::cout << "Time move move = " << time << " [s], t-ref= " << time - time_ctor << " [s]\n";
        }
        {
            auto timer = af::timer::start();
            for (auto n = 0; n < N; ++n) {
                auto data = af::constant(42, dims, type);
                auto data_copy_copy = Cp_or_Mv(data);
            }
            double time = af::timer::stop(timer);
            std::cout << "Time single copy = " << time << " [s], t-ref= " << time - time_ctor << " [s]\n";
        }
        {
            auto timer = af::timer::start();
            for (auto n = 0; n < N; ++n) {
                auto data = af::constant(42, dims, type);
                auto data_copy_copy = Cp_or_Mv(std::move(data));
            }
            double time = af::timer::stop(timer);
            std::cout << "Time single move = " << time << " [s], t-ref= " << time - time_ctor << " [s]\n";
        }
    }
    if (false) {
        // Separation of lambdas by cp/mv not neccessary/useless:

        // lambda passing data as lval
        const auto lambda_cp = [N, &dims, time_ctor](auto callback, const auto& string) {
            auto timer = af::timer::start();
            for (auto n = 0; n < N; ++n) {
                auto data = af::constant(42, dims, type);
                // af::print("before calback cp", data);
                auto data_lambda = callback(data);
                // af::print("after calback cp", data);
            }
            double time = af::timer::stop(timer);
            std::cout << "Lambda_callsite_cp: Time " << string << " = " << time << " [s], t-ref= " << time - time_ctor
                      << " [s]\n";
        };

        // lambda passing data as rval
        const auto lambda_mv = [N, &dims, time_ctor](auto callback, const auto& string) {
            auto timer = af::timer::start();
            for (auto n = 0; n < N; ++n) {
                auto data = af::constant(42, dims, type);
                auto data_lambda = callback(std::move(data));
            }
            double time = af::timer::stop(timer);
            std::cout << "Lambda_callsite_mv: Time " << string << " = " << time << " [s], t-ref= " << time - time_ctor
                      << " [s]\n";
        };

        // Minimal impact lambdas, passed lambda forward-moves arg by rval ref:
        std::cout << std::endl;
        lambda_cp([](auto&& data) { return CpCp_or_MvCp{data}; }, "Copy Copy");
        // lambda_cp([](auto&& data) { return CpCp_or_MvCp{std::move(data)}; }, "Copy Copy"); // NOTE This moves data
        // anyway, even though not stated in body
        lambda_mv([](auto&& data) { return CpCp_or_MvCp{std::move(data)}; }, "Move Copy");
        lambda_cp([](auto&& data) { return CpMv_or_MvMv{data}; }, "Copy Move");
        // lambda_cp([](auto&& data) { return CpMv_or_MvMv{std::move(data)}; }, "Copy Move");
        lambda_mv([](auto&& data) { return CpMv_or_MvMv{std::move(data)}; }, "Move Move");
        lambda_cp([](auto&& data) { return Cp_or_Mv{data}; }, "Single Copy");
        lambda_mv([](auto&& data) { return Cp_or_Mv{std::move(data)}; }, "Single Move");

        // // Minimal impact lambdas, passed lambda forward-moves arg by rval ref:
        // // This way of lambda arguments causes an additional copy, bat for timing:
        // std::cout << std::endl;
        // lambdacp([](auto const& data) { return CpCp_or_MvCp{std::move(data)}; }, "Copy Copy");
        // lambdamv([](auto const& data) { return CpCp_or_MvCp{std::move(data)}; }, "Move Copy");
        // lambdacp([](auto const& data) { return CpMv_or_MvMv{std::move(data)}; }, "Copy Move");
        // lambdamv([](auto const& data) { return CpMv_or_MvMv{std::move(data)}; }, "Move Move");
        // lambdacp([](auto const& data) { return Cp_or_Mv{std::move(data)}; }, "Single Copy");
        // lambdamv([](auto const& data) { return Cp_or_Mv{std::move(data)}; }, "Single Move");
    }

    return 0;
}
