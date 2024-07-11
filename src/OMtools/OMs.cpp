#include <iostream>
#include "OMs.hpp"

// ===============
// Chirotope<R, N>
// ===============

template<int R, int N>
std::ostream& operator<<(std::ostream& os, const Chirotope<R, N>& chi) {
    for (auto i = 0; i < chi.NR_INT32 - 1; i++) {
        for (auto c = 0; c < 32; c++) {
            if (chi.plus[i] & ((uint32_t)1 << c))
                os << '+';
            else if (chi.minus[i] & ((uint32_t)1 << c))
                os << '-';
            else
                os << '0';
        }
    }
    for (auto c = 0; c < chi.NR_REMAINING_BITS; c++) {
        if (chi.plus[chi.NR_INT32 - 1] & ((uint32_t)1 << c))
            os << '+';
        else if (chi.minus[chi.NR_INT32 - 1] & ((uint32_t)1 << c))
            os << '-';
        else
            os << '0';
    }
    return os;
}

template<int R, int N>
std::ofstream& operator<<(std::ofstream& of, const Chirotope<R, N>& chi) {
    for (auto i = 0; i < chi.NR_INT32 - 1; i++) {
        for (auto c = 0; c < 32; c++) {
            if (chi.plus[i] & ((uint32_t)1 << c))
                of << '+';
            else if (chi.minus[i] & ((uint32_t)1 << c))
                of << '-';
            else
                of << '0';
        }
    }
    for (auto c = 0; c < chi.NR_REMAINING_BITS; c++) {
        if (chi.plus[chi.NR_INT32 - 1] & ((uint32_t)1 << c))
            of << '+';
        else if (chi.minus[chi.NR_INT32 - 1] & ((uint32_t)1 << c))
            of << '-';
        else
            of << '0';
    }
    return of;
}
