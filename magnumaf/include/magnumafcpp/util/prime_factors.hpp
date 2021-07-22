#pragma once
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace magnumafcpp::util {

/// Returns all prime factors of a given number n
/// Special cases: n == 0 returns {0}, n == 1 returns {1}
inline std::vector<unsigned> prime_factors(unsigned n) {
    if (n == 0) {
        return {0};
    } else if (n == 1) {
        return {1};
    } else {
        std::vector<unsigned> primes;
        // Print the number of 2s that divide n
        while (n % 2 == 0) {
            primes.push_back(2);
            n = n / 2;
        }

        // n is odd and will only be further dividable by odd numbers
        // thus we iterate over all odd numbers up to sqrt(n)
        for (int i = 3; i <= std::sqrt(n); i = i + 2) {
            // While i divides n, print i and divide n
            while (n % i == 0) {
                primes.push_back(i);
                n = n / i;
            }
        }

        // when we got here n is either 1 or the final prime
        if (n != 1) {
            primes.push_back(n);
        }
        return primes;
    }
}

inline unsigned max_of_prime_factors(unsigned n) {
    auto primes = prime_factors(n);
    return *std::max_element(primes.begin(), primes.end());
}

} // namespace magnumafcpp::util
