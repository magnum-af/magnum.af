#include "arrayfire.h"
#include <array>
#include <iostream>
#include <string>

///
/// Struct for version comparison.
/// Using: major_.minor_.revision_.build__, e.g: 3.7.2.0
///
struct Version {
    int major_ = 0, minor_ = 0, revision_ = 0, build__ = 0;

    Version(std::string version) {
        std::sscanf(version.c_str(), "%d.%d.%d.%d", &major_, &minor_, &revision_, &build__);
    }

    Version(int major_, int minor_ = 0, int revision_ = 0, int build__ = 0)
        : major_(major_), minor_(minor_), revision_(revision_), build__(build__) {}

    bool operator<(const Version& other) {
        if (major_ < other.major_) {
            return true;
        } else if (minor_ < other.minor_) {
            return true;
        } else if (revision_ < other.revision_) {
            return true;
        } else if (build__ < other.build__) {
            return true;
        } else {
            return false;
        }
    }

    bool operator>(const Version& other) {
        if (major_ > other.major_) {
            return true;
        } else if (minor_ > other.minor_) {
            return true;
        } else if (revision_ > other.revision_) {
            return true;
        } else if (build__ > other.build__) {
            return true;
        } else {
            return false;
        }
    }

    bool operator>=(const Version& other) {
        if (major_ > other.major_) {
            return true;
        } else if (minor_ > other.minor_) {
            return true;
        } else if (revision_ > other.revision_) {
            return true;
        } else if (build__ > other.build__) {
            return true;
        } else {
            return major_ == other.major_ && minor_ == other.minor_ && revision_ == other.revision_ &&
                   build__ == other.build__;
        }
    }

    bool operator==(const Version& other) {
        return major_ == other.major_ && minor_ == other.minor_ && revision_ == other.revision_ &&
               build__ == other.build__;
    }

    std::ostream& operator<<(std::ostream& stream) {
        stream << this->major_;
        stream << '.';
        stream << this->minor_;
        stream << '.';
        stream << this->revision_;
        stream << '.';
        stream << this->build__;
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

