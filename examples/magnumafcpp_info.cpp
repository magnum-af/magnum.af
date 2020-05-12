#include "arrayfire.h"
#include "magnum_af.hpp"
#include "magnumafConfig.hpp"
#include "magnumafConfig_git.hpp"

using namespace magnumafcpp;

int main() {
    af::info();

    std::cout << "magnumafcpp_VERSION_MAJOR = " << magnumafcpp_VERSION_MAJOR
              << std::endl;
    std::cout << "magnumafcpp_VERSION_MINOR = " << magnumafcpp_VERSION_MINOR
              << std::endl;
    std::cout << "magnumafcpp_VERSION_PATCH = " << magnumafcpp_VERSION_PATCH
              << std::endl;
    std::cout << "magnumafcpp_VERSION_TWEAK = " << magnumafcpp_VERSION_TWEAK
              << std::endl;
    std::cout << "magnumafcpp_VERSION       = " << magnumafcpp_VERSION
              << std::endl;

    std::cout << "GIT_RETRIEVED_STATE       = " << GIT_RETRIEVED_STATE
              << std::endl;
    std::cout << "GIT_HEAD_SHA1             = " << GIT_HEAD_SHA1 << std::endl;
    std::cout << "GIT_IS_DIRTY              = " << GIT_IS_DIRTY << std::endl;
    std::cout << "GIT_AUTHOR_NAME           = " << GIT_AUTHOR_NAME << std::endl;
    std::cout << "GIT_AUTHOR_EMAIL          = " << GIT_AUTHOR_EMAIL
              << std::endl;
    std::cout << "GIT_COMMIT_DATE_ISO8601   = " << GIT_COMMIT_DATE_ISO8601
              << std::endl;
    std::cout << "GIT_COMMIT_SUBJECT        = " << GIT_COMMIT_SUBJECT
              << std::endl;
    std::cout << "GIT_COMMIT_BODY           = " << GIT_COMMIT_BODY << std::endl;
    return 0;
}
