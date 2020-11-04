#pragma once
#include "arrayfire.h"
#include <filesystem>
#include <optional>

namespace magnumafcpp::util {
class CacheManager {
  public:
    CacheManager(std::filesystem::path p = std::filesystem::temp_directory_path() / "magnumaf.cache",
                 std::uintmax_t max_size_in_byte = 2'000'000'000, std::uintmax_t shrink_size_in_byte = 200'000'000);

    ~CacheManager();

    std::optional<af::array> get_array_if_existent(const std::string& filename) const;
    void write_array(const af::array& a, const std::string& filename, const std::string& key = "") const;

  private:
    const std::filesystem::path cache_folder;
    const bool verbose{true};
    const std::uintmax_t max_size_in_byte;
    const std::uintmax_t shrink_size_in_byte;

    void shrink_cache_if_gt_maxsize();
};

} // namespace magnumafcpp::util
