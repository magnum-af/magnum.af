#pragma once
#include <algorithm>
#include <bits/stdc++.h>
#include <vector>
namespace magnumafcpp {
namespace util {

// adaped from: https://www.geeksforgeeks.org/print-all-prime-factors-of-a-given-number/
/// Returns all prime factors of a given number n
/// Special cases: n == 0 returns {0}, n == 1 returns {1}
std::vector<unsigned> inline prime_factors(unsigned n) {
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

        // n must be odd at this point. So we can skip
        // one element (Note i = i +2)
        for (int i = 3; i <= sqrt(n); i = i + 2) {
            // While i divides n, print i and divide n
            while (n % i == 0) {
                primes.push_back(i);
                n = n / i;
            }
        }

        // This condition is to handle the case when n
        // is a prime number greater than 2
        if (n > 2) {
            primes.push_back(n);
        }
        return primes;
    }
}

unsigned inline max_of_prime_factors(unsigned n) {
    auto primes = prime_factors(n);
    return *std::max_element(primes.begin(), primes.end());
}

} // namespace util
} // namespace magnumafcpp
