#include "cache_manager.hpp"
#include <algorithm>
#include <filesystem>
#include <iostream>
#include <numeric>
#include <utility>

namespace magnumaf::util {

namespace fs = std::filesystem;

// getting files in search path ignoring directories
auto get_files_in_dir(const fs::path& search_path) {
    std::vector<fs::directory_entry> container;
    auto dir_iter = fs::directory_iterator(search_path);
    std::copy_if(fs::begin(dir_iter), fs::end(dir_iter), std::back_inserter(container),
                 [](const auto& elem) { return !fs::is_directory(elem); });
    return container;
}

// accumulated size of files in list, ignoring dirs
std::uintmax_t accum_size_in_byte(const std::vector<fs::directory_entry>& vec) {
    return std::accumulate(std::begin(vec), std::end(vec), std::uintmax_t{0}, [](std::uintmax_t i, const auto& elem) {
        return fs::is_directory(elem) ? i : i + fs::file_size(elem);
    });
}

void sort_by_write_time_newest_first(std::vector<fs::directory_entry>& container) {
    std::sort(std::begin(container), std::end(container),
              [](const auto& lhs, const auto& rhs) { return fs::last_write_time(lhs) > fs::last_write_time(rhs); });
}

enum class RunPolicy { silent, verbose, test };
/// remove elements at back of list until accum size is below /param size_in_byte
/// /param filelist arbitrary-sorted list of files, will delete at the back
void remove_last_file_until_below(std::uintmax_t size_in_byte, std::vector<fs::directory_entry>& filelist,
                                  RunPolicy r = RunPolicy::silent) {
    while (accum_size_in_byte(filelist) > size_in_byte and not filelist.empty()) {
        if (r != RunPolicy::silent) {
            std::cout << (r == RunPolicy::test ? "would remove file: " : "removing file: ") << filelist.back() << '\n';
        }
        if (r != RunPolicy::test) {
            fs::remove(filelist.back()); // delete last (oldest) file in list:
        }
        filelist.pop_back(); // remove file-entry from list
    }
}

CacheManager::CacheManager(bool verbose, fs::path p, std::uintmax_t max_size_in_byte,
                           std::uintmax_t shrink_size_in_byte)
    : verbose(verbose), cache_folder(std::move(p)), max_size_in_byte(max_size_in_byte),
      shrink_size_in_byte(shrink_size_in_byte) {
    // create dir if not existing
    if (!fs::exists(cache_folder)) {
        fs::create_directories(cache_folder);
        // setting file permissions to linux '777'
        fs::permissions(cache_folder, fs::perms::all);
    }
}

CacheManager::~CacheManager() {
    try {
        // checking cleanup
        shrink_cache_if_gt_maxsize();
        // remove cache_folder if empty
        if (fs::is_empty(cache_folder)) {
            fs::remove(cache_folder);
        }
    } catch (const af::exception& e) {
        std::cout << e.what() << "\33[1;31mWarning:\33[0m clearing chache failed." << std::endl;
    }
}

std::optional<af::array> CacheManager::get_array_if_existent(const std::string& filename) const {
    const auto file_path = cache_folder / filename;
    if (fs::exists(file_path)) {
        if (verbose) {
            std::cout << "\33[0;32mInfo:\33[0m Reading array from " << file_path << std::endl;
        }
        try {
            return af::readArray(file_path.c_str(), "");
        } catch (const af::exception& e) {
            std::cout << e.what()
                      << "\33[1;31mWarning:\33[0m af::readArrayCheck failed. Omit reading demag tensor, calculating it "
                         "instead."
                      << std::endl;
            return {};
        }

    } else {
        return {};
    }
}
void CacheManager::write_array(const af::array& a, const std::string& filename, const std::string& key) const {
    if (verbose) {
        std::cout << "\33[0;32mInfo:\33[0m Saving array to " << cache_folder / filename << std::endl;
    }
    try {
        const auto array_filename = cache_folder / filename;
        af::saveArray(key.c_str(), a, (array_filename).c_str());
        fs::permissions(array_filename, fs::perms::others_write,
                        fs::perm_options::add); // changes 'rw-rw-r--' to 'rw-rw-rw-'
    } catch (const af::exception& e) {
        std::cout << e.what() << "\33[1;31mWarning:\33[0m af::saveArray failed, omit saving demag tensor." << std::endl;
    }
}

void CacheManager::shrink_cache_if_gt_maxsize() const {
    auto vec = get_files_in_dir(cache_folder);
    if (accum_size_in_byte(vec) > max_size_in_byte) {
        sort_by_write_time_newest_first(vec);
        remove_last_file_until_below(shrink_size_in_byte, vec, verbose ? RunPolicy::verbose : RunPolicy::silent);
    }
}

} // namespace magnumaf::util
