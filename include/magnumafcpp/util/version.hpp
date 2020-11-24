#pragma once
#include "arrayfire.h"
#include <array>
#include <iostream>
#include <regex>
#include <string>

namespace magnumafcpp {
///
/// Struct for version comparison.
/// Using: major_.minor_.revision_.build_, e.g: 3.7.2.0
///
struct Version {
    int major_ = 0, minor_ = 0, revision_ = 0, build_ = 0;

    Version(std::string version) { std::sscanf(version.c_str(), "%d.%d.%d.%d", &major_, &minor_, &revision_, &build_); }

    Version(int major_, int minor_ = 0, int revision_ = 0, int build_ = 0)
        : major_(major_), minor_(minor_), revision_(revision_), build_(build_) {}

    bool operator==(const Version& other) {
        return major_ == other.major_ && minor_ == other.minor_ && revision_ == other.revision_ &&
               build_ == other.build_;
    }

    bool operator<(const Version& other) {
        if (major_ < other.major_) {
            return true;
        } else if (minor_ < other.minor_) {
            return true;
        } else if (revision_ < other.revision_) {
            return true;
        } else if (build_ < other.build_) {
            return true;
        } else {
            return false;
        }
    }

    bool operator<=(const Version& other) {
        if (*this < other) {
            return true;
        } else {
            return *this == other;
        }
    }

    bool operator>(const Version& other) {
        return !(*this <= other);
    }

    bool operator>=(const Version& other) {
        return !(*this < other);
    }

    friend std::ostream& operator<<(std::ostream& stream, const Version& ver) {
        stream << ver.major_;
        stream << '.';
        stream << ver.minor_;
        stream << '.';
        stream << ver.revision_;
        stream << '.';
        stream << ver.build_;
        return stream;
    }
};

/// Returning arrayfire major_, minor_, patch version
std::array<int, 3> af_version() {
    const std::string s = af::infoString();
    int major_ = std::stoi(s.substr(11, 1));
    int minor_ = std::stoi(s.substr(13, 1));
    int patch = std::stoi(s.substr(15, 1));
    return {major_, minor_, patch};
}

/// Returning arrayfire version string
std::string af_version_string() {
    const std::string s = af::infoString();
    return {s.substr(11, 5)};
}

/// Search version number from string, returning as Version object
// NOTE: causes std::bad_alloc if af::infoString() is invoked, presumably due to raw pointer handling
Version version_from_string(std::string s,
                            std::regex pattern = std::regex("([0-9]+).([0-9]+).([0-9]+)", std::regex::awk),
                            bool verbose = false) {
    std::smatch match;
    if (std::regex_search(s, match, pattern)) {
        if (verbose) {
            for (unsigned i = 0; i < match.size(); ++i) {
                std::cout << "match[" << i << "]=" << match[i] << std::endl;
            }
        }
        return Version(match[0]);
    } else {
        std::cout << "Warning: No regex match found, returning Version('-1.-1.-1.-1') instead." << std::endl;
        return Version("-1.-1.-1.-1");
    }
}

} // namespace magnumafcpp
