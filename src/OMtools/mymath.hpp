#pragma once
#include<type_traits>
#include <cstdint>

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

// NchooseK<N,k>::array lists all k-element subsets of {0..N-1}.
// Auxulary constants are provided to make querying this list a
// smoother experience. All constants are computed compile-time.
template<typename T, T N, T k,
    std::enable_if_t<std::is_integral_v<T>, bool> = true>
struct NchooseK {
    // `(N choose R)=N!/(R!*(N-R)!)`, the number of distinct `R`-element 
    // subsets of `0..N-1`.
    constexpr static const auto NR_RTUPLES = binomial_coefficient(N,k);
    // Number of (unsigned) 32-bit integers needed to have 1 bit for each 
    // `R` element subset of `0..N-1`.
    constexpr static const auto NR_INT32 = division_rounded_up(NR_RTUPLES, 32);
    // Number of bits used in the last 32-bit unsigned integer used to 
    // store the characteristic vector.
    constexpr static const auto NR_REMAINING_BITS = NR_RTUPLES % 32;
    // `array[i][t]` is the `t`th element of the `i`th `k`-element 
    // subset of `0..N-1`. This is computed at compile time.
    constexpr static const auto array{[]() constexpr{
        std::array<std::array<T, k>, NR_RTUPLES> a{}; 
        for (T i = 0; i < k; i++) {
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
    // `contained[e][i]` is true if and only if the `i`th `k`-element
    // subset of `0..N-1` contains the element `e`.
    constexpr static const auto contained{[]() constexpr{
        std::array<std::array<bool, NR_RTUPLES>, N> c{}; 
        for (T e = 0; e < k; e++) {
            for (T i = 0; i < NR_RTUPLES; i++) {
                c[array[i][e]][i] = true;
            }
        }
        return c;
    }()};
    // The `r`th lowest bit of `contained_mask32[e][idx]` is 1
    // if and only if the `idx * 32 + r`th `k`-element subset
    // of `0..N-1` contains the element `e`.
    constexpr static const auto contained_mask32{[]() constexpr{
        std::array<std::array<uint32_t, NR_INT32>, N> m{};
        for (T e = 0; e < N; e++) {
            for (auto idx = 0; idx < NR_INT32 - 1; idx++) {
                for (auto r = 0; r < 32; r++) {
                    m[e][idx] |= ((uint32_t)contained[e][idx * 32 + r]
                        & (uint32_t)1) << r;
                }
            }
            for (auto r = 0; r < NR_REMAINING_BITS; r++) {
                m[e][NR_INT32 - 1] |= 
                    ((uint32_t)contained[e][(NR_INT32 - 1) * 32 + r]
                    & (uint32_t)1) << r;
            }
        }
        return m;
    }()};
};