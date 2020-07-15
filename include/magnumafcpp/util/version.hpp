#include "arrayfire.h"
#include <array>
#include <iostream>
#include <string>

///
/// Struct for version comparison.
/// Using: major_.minor.revision.build, e.g: 3.7.2.0
///
struct Version {
    int major_ = 0, minor = 0, revision = 0, build = 0;

    Version(std::string version) { std::sscanf(version.c_str(), "%d.%d.%d.%d", &major_, &minor, &revision, &build); }

    Version(int major_, int minor = 0, int revision = 0, int build = 0)
        : major_(major_), minor(minor), revision(revision), build(build) {}

    bool operator<(const Version& other) {
        if (major_ < other.major_) {
            return true;
        } else if (minor < other.minor) {
            return true;
        } else if (revision < other.revision) {
            return true;
        } else if (build < other.build) {
            return true;
        } else {
            return false;
        }
    }

    bool operator>(const Version& other) {
        if (major_ > other.major_) {
            return true;
        } else if (minor > other.minor) {
            return true;
        } else if (revision > other.revision) {
            return true;
        } else if (build > other.build) {
            return true;
        } else {
            return false;
        }
    }

    bool operator>=(const Version& other) {
        if (major_ > other.major_) {
            return true;
        } else if (minor > other.minor) {
            return true;
        } else if (revision > other.revision) {
            return true;
        } else if (build > other.build) {
            return true;
        } else {
            return major_ == other.major_ && minor == other.minor && revision == other.revision && build == other.build;
        }
    }

    bool operator==(const Version& other) {
        return major_ == other.major_ && minor == other.minor && revision == other.revision && build == other.build;
    }

    std::ostream& operator<<(std::ostream& stream) {
        stream << this->major_;
        stream << '.';
        stream << this->minor;
        stream << '.';
        stream << this->revision;
        stream << '.';
        stream << this->build;
        return stream;
    }
};

/// Returning arrayfire major_, minor, patch version
std::array<int, 3> af_version() {
    const std::string s = af::infoString();
    int major_ = std::stoi(s.substr(11, 1));
    int minor = std::stoi(s.substr(13, 1));
    int patch = std::stoi(s.substr(15, 1));
    return {major_, minor, patch};
}

/// Returning arrayfire version string
std::string af_version_string() {
    const std::string s = af::infoString();
    return {s.substr(11, 5)};
}

