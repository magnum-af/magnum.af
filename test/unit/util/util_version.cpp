#include "util/version.hpp"
#include <gtest/gtest.h>

using namespace magnumafcpp;

TEST(Util, Version) {
    EXPECT_TRUE(Version("3.7.8.0") == Version("3.7.8.0"));
    EXPECT_TRUE(Version("3.7.8.0") == Version("3.7.8"));
    EXPECT_FALSE(Version("3.7.8.0") < Version("3.7.8"));
    EXPECT_FALSE(Version("3.7.9") < Version("3.7.8"));
    EXPECT_TRUE(Version("3") < Version("3.7.9"));
    EXPECT_TRUE(Version("1.7.9") < Version("3.1"));

    EXPECT_EQ(Version("3.6.4") < Version("3.7.0"), true);
    EXPECT_EQ(Version("3.7.3") < Version("3.7.0"), false);
    EXPECT_EQ(Version("3.7.0") < Version("3.7.0"), false);

    EXPECT_EQ(Version("3.6.4") <= Version("3.7.0"), true);
    EXPECT_EQ(Version("3.7.3") <= Version("3.7.0"), false);
    EXPECT_EQ(Version("3.7.0") <= Version("3.7.0"), true);

    EXPECT_EQ(Version("3.6.4") > Version("3.7.0"), false);
    EXPECT_EQ(Version("3.7.3") > Version("3.7.0"), true);
    EXPECT_EQ(Version("3.7.0") > Version("3.7.0"), false);

    EXPECT_EQ(Version("3.6.4") >= Version("3.7.0"), false);
    EXPECT_EQ(Version("3.7.3") >= Version("3.7.0"), true);
    EXPECT_EQ(Version("3.7.0") >= Version("3.7.0"), true);
}

TEST(Util, version_from_string) {
    EXPECT_EQ(Version(version_from_string("Test 1.2.3 test")) == Version("1.2.3"), true);
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
