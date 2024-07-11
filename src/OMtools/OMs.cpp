#include <iostream>
#include <string>
#include "OMs.hpp"

// ===============
// Chirotope<R, N>
// ===============

template<int R, int N>
constexpr Chirotope<R, N>::Chirotope(std::string from): plus{}, minus{} {
    read(from);
}

template<int R, int N>
constexpr Chirotope<R, N> Chirotope<R,N>::inverse() const {
    auto chi = Chirotope<R, N>();
    for (auto i = 0; i < NR_INT32; i++) {
        chi.plus[i] = minus[i];
        chi.minus[i] = plus[i];
    }
    return chi;
} 

// Here "idx >> 5" is just "idx / 32" and "idx & 31" is "idx % 32"

template<int R, int N>
constexpr bool Chirotope<R, N>::get_plus(int idx) const {
    return plus[idx >> 5] & ((uint32_t)1 << (idx & 31));
}

template<int R, int N>
constexpr bool Chirotope<R, N>::get_minus(int idx) const {
    return minus[idx >> 5] & ((uint32_t)1 << (idx & 31));
}

template<int R, int N>
constexpr bool Chirotope<R, N>::is_basis(int idx) const {
    return get_plus(idx) || get_minus(idx);
}

template<int R, int N>
constexpr char Chirotope<R, N>::get(int idx) const {
    return evaluate(idx);
}

template<int R, int N>
constexpr char Chirotope<R, N>::evaluate(int idx) const {
    if (plus[idx >> 5] & ((uint32_t)1 << (idx & 31))) {
        return '+';
    } else if (minus[idx >> 5] & ((uint32_t)1 << (idx & 31))) {
        return '-';
    } else {
        return '0';
    }
}

template<int R, int N>
constexpr void Chirotope<R, N>::set_plus(int idx, bool value) {
    plus[idx >> 5] = plus[idx >> 5] & ~((uint32_t)1 << (idx & 31)) | uint32_t(value) << (idx & 31);
}

template<int R, int N>
constexpr void Chirotope<R, N>::set_minus(int idx, bool value) {
    minus[idx >> 5] = minus[idx >> 5] & ~((uint32_t)1 << (idx & 31)) | uint32_t(value) << (idx & 31);
}

template<int R, int N>
constexpr void Chirotope<R, N>::set(int i, int r, char c) {
    switch (c)
    {
    case '+':
        plus[i] |= (uint32_t)1 << r;
        minus[i] &= ~((uint32_t)1 << r);
        break;
    case '-':
        minus[i] |= (uint32_t)1 << r;
        plus[i] &= ~((uint32_t)1 << r);
        break;
    case '0':
        plus[i] &= ~((uint32_t)1 << r);
        minus[i] &= ~((uint32_t)1 << r);
        break;
    default:
        throw std::invalid_argument("Invalid character! "
        "The only valid characters are '0', '-', and '+', but the inputted "
        "character was '"+std::to_string(c)+"'.");
        break;
    }
}

template<int R, int N>
constexpr void Chirotope<R, N>::set(int idx, char c) {
    set(idx >> 5, idx & 31, c);
}

template<int R, int N>
constexpr bool Chirotope<R, N>::is_zero() const {
    for (auto i = 0; i < NR_INT32; i++) {
        if (plus[i] != 0 || minus[i] != 0) return true;
    }
    return false;
}

template<int R, int N>
constexpr int Chirotope<R, N>::countbases() const {
    int count = 0;
    for (auto i = 0; i < NR_RTUPLES; i++) {
        count += int(is_basis(i));
    }
    return count;
}

template<int R, int N>
constexpr Chirotope<R,N>& Chirotope<R, N>::read(const std::string& str) {
    if (str.length() != NR_RTUPLES) throw std::invalid_argument("The input string is of incorrect length! "
        "Input length: "+std::to_string(str.length())+", expected length: ("+std::to_string(N)+" "
        "choose "+std::to_string(R)+") = "+std::to_string(NR_RTUPLES)+"."
        " Input string: <"+str+">");
    for (auto i = 0; i < NR_INT32 - 1; i++) {
        for (auto c = 0; c < 32; c++) {
            set(i, c, str[i * 32 + c]);
        }
    }
    for (auto c = 0; c < NR_REMAINING_BITS; c++) {
        set(NR_INT32 - 1, c, str[(NR_INT32 - 1) * 32 + c]);
    }
    return (*this);
}

template<int R, int N>
constexpr bool Chirotope<R, N>::weak_maps_to(const Chirotope& chi) const {
    for (auto i = 0; i < NR_INT32; i++) {
        if (~plus[i] & chi.plus[i] | ~minus[i] & chi.minus[i]) return false;
    }
    return true;
}

template<int R, int N>
constexpr bool Chirotope<R, N>::OM_weak_maps_to(const Chirotope& chi) const {
    bool is_wm = true;
    for (auto i = 0; i < NR_INT32; i++) {
        if (~plus[i] & chi.minus[i] | ~minus[i] & chi.plus[i]) {
            is_wm = false;
            break;
        }
    }
    if (is_wm) return true;
    return weak_maps_to(chi);
}

template<int R, int N>
constexpr bool Chirotope<R, N>::operator==(const Chirotope& chi) const {
    for (auto i = 0; i < NR_INT32; i++) {
        if (plus[i] != chi.plus[i] || minus[i] != chi.minus[i])
            return false;
    }
    return true;
}

template<int R, int N>
constexpr bool Chirotope<R, N>::operator!=(const Chirotope& chi) const {
    for (auto i = 0; i < NR_INT32; i++) {
        if (plus[i] != chi.plus[i] || minus[i] != chi.minus[i])
            return true;
    }
    return false;
}

template<int R, int N>
constexpr bool Chirotope<R, N>::is_same_OM_as(const Chirotope& chi) const {
    bool is_same = true;
    for (auto i = 0; i < NR_INT32; i++) {
        if (plus[i] != chi.minus[i] || minus[i] != chi.plus[i]) {
            is_same = false;
            break;
        }
    }
    if (is_same) return true;
    return (*this) == chi;
}

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
