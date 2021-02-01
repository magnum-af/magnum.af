#include "cache_manager.hpp"
#include <filesystem>
#include <gtest/gtest.h>
using namespace magnumafcpp::util;
namespace fs = std::filesystem;

const auto cm_path = fs::temp_directory_path() / "test_magnumaf_cache";

TEST(util_cache_manager, ctor_dtor) {
    {
        auto cm = CacheManager{false, cm_path, 0, 0};
        EXPECT_TRUE(fs::exists(cm_path));
    }
    EXPECT_FALSE(fs::exists(cm_path));
}

const af::array a = af::iota(10, f64);
const af::array b = af::iota(5, f64);

TEST(util_cache_manager, write_array) {
    {
        auto cm = CacheManager{false, cm_path};
        cm.write_array(a, "a");
        cm.write_array(b, "b");
        EXPECT_TRUE(fs::exists(cm_path / "a"));
        EXPECT_TRUE(fs::exists(cm_path / "b"));
    }
    EXPECT_TRUE(fs::exists(cm_path));
}

TEST(util_cache_manager, read_array) {
    {
        auto cm = CacheManager{false, cm_path};
        auto a_op = cm.get_array_if_existent("a");
        auto b_op = cm.get_array_if_existent("b");
        af::array a_in = a_op.value();
        af::array b_in = b_op.value();
        // check if read arrays are equal to initial arrays:
        EXPECT_EQ(af::sum(a != a_in).scalar<unsigned>(), 0);
        EXPECT_EQ(af::sum(b != b_in).scalar<unsigned>(), 0);
    }
    EXPECT_TRUE(fs::exists(cm_path));
}

// test is dtor cleans up each file if max size is 0:
TEST(util_cache_manager, dtor_cleanup_if_dir_populated) {
    auto filepath = cm_path / "testfile.txt";
    {
        EXPECT_TRUE(fs::exists(cm_path / "a"));
        EXPECT_TRUE(fs::exists(cm_path / "b"));
        auto cm = CacheManager{false, cm_path, 0, 0};
    }
    EXPECT_FALSE(fs::exists(cm_path));
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
