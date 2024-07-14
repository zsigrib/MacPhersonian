#pragma once
#include<type_traits>

// ====================
// FORWARD DECLARATIONS
// ====================
// #learningcpp: as these functions are constexpr, and they
//    are intented to be run at compile time, their *definitions*
//    must be available in compile time as well, and as compiling
//    happens file-by-file, this means the definitions must be
//    contained in the same header.


template<typename T,
    std::enable_if_t<std::is_integral_v<T>, bool> = true>
constexpr T factorial(T N);

template<typename T,
    std::enable_if_t<std::is_integral_v<T>, bool> = true>
constexpr T binomial_coefficient(T N, T k);

template<typename T,
    std::enable_if_t<std::is_integral_v<T>, bool> = true>
constexpr T division_rounded_up(T total, T divisor);

// ===========
// DEFINITIONS
// ===========

template<typename T,
    std::enable_if_t<std::is_integral_v<T>, bool>>
constexpr T factorial(T N) {
    if (N == 0) { 
        return 1; 
    } else if (N < 0) { 
        throw std::invalid_argument("Factorial of negative values is undefined; "+std::to_string(N)+" is an invalid input.");
    }
    T total = 1;
    for (T i=N; i > 1; i--) {
        total *= i;
    }
    return total;
}

template<typename T,
    std::enable_if_t<std::is_integral_v<T>, bool>>
constexpr T binomial_coefficient(T N, T k) {
    if (N < 0 || k < 0 || N < k) { 
        return 0; 
    } else if (k == 0) { 
        return 1; 
    } else {
        return (N * binomial_coefficient(N-1, k-1)) / k;
    }
}

template<typename T,
    std::enable_if_t<std::is_integral_v<T>, bool>>
constexpr T division_rounded_up(T total, T divisor) {
    if (total % divisor == 0) {
        return total / divisor;
    } else {
        return total / divisor + 1;
    }
}

// NchooseK<N,k>::array lists all k-element subsets of {0..N-1}, 
// and is calculated compile time.
template<typename T, T N, T k,
    std::enable_if_t<std::is_integral_v<T>, bool> = true>
struct NchooseK {
    constexpr static const auto array{[]() constexpr{
        std::array<std::array<T, k>, binomial_coefficient(N,k)> a{}; 
        for (int i = 0; i < k; i++) {
            a[0][i] = i;
        }
        T r = k -1;
        T s = 1;
        while(true) {
            r = k - 1;
            while(r >= 0 && a[s - 1][r] >= N-k+r)
                r--;
            if (r == -1)
                break;
            for (auto i = 0; i < r; i++)
                a[s][i] = a[s-1][i];
            for (auto i = r; i < k; i++)
                a[s][i] = a[s-1][r] + i - r + 1;
            s++;
        }
        return a;
    }()};
};