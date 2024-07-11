#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <cstdint>
#include "mymath.hpp"

//writebits? -> print an integer as a 0-1 sequence

// ////////////////////
// ====================
// FORWARD DECLARATIONS
// ====================
// ////////////////////

// ===============
// Chirotope<R, N>
// ===============

template<int R, int N> struct Chirotope;

template<int R, int N>
std::ostream& operator<<(std::ostream&, const Chirotope<R, N>&);

template<int R, int N>
std::ofstream& operator<<(std::ofstream&, const Chirotope<R, N>&);

template<int R, int N>
struct Chirotope
{
    // `(N choose R)=N!/(R!*(N-R)!)`, the number of distinct `R`-element 
    // subsets of `0..N-1`.
    static constexpr int NR_RTUPLES = binomial_coefficient(N, R);
    // Number of (unsigned) 32-bit integers needed to have 1 bit for each 
    // `R` element subset of `0..N-1`.
    static constexpr int NR_INT32 = division_rounded_up(NR_RTUPLES, 32);
    // Number of bits used in the last 32-bit unsigned integer used to 
    // store the chirotope.
    static constexpr int NR_REMAINING_BITS = NR_RTUPLES % 32;
    // Use `RTUPLES_LIST().array[i][k]` to access the `k`th element of 
    // the `i`th `R`-element subset of `0..N-1`. `RTUPLES_LIST().array` 
    // evaluates at compile-time.
    using RTUPLES_LIST = NchooseK<int, N, R>;
    
    // The `c`th least significant bit of `plus[i]` is 1 if and only if 
    // the chirotope evaluates to `+` on the `R`-tuple stored in
    // `RTUPLES_LIST().array[i * 32 + c]`.
    uint32_t plus[NR_INT32];
    // The `c`th least significant bit of `minus[i]` is 1 if and only if 
    // the chirotope evaluates to `-` on the `R`-tuple stored in
    // `RTUPLES_LIST().array[i * 32 + c]`.
    uint32_t minus[NR_INT32];

    // Initializes `plus` and `minus` to arrays of 0s. Note that this
    // does not satisfy chirotope axiom (B0), so it is not a valid
    // chirotope.
    Chirotope(): plus{}, minus{} {}

    // Return the inverse of this chirotope, i.e. if this chirotope is
    // denoted by `c`, then `mc`, for which 
    // - `c[i] == '0'` implies `mc[i] == '0'`,
    // - `c[i] == '+'` implies `mc[i] == '-'`, and
    // - `c[i] == '-'` implies `mc[i] == '+'`. 
    constexpr Chirotope inverse() const;
    
    // Returns whether the chirotope evaluates to `+` on the `R`-tuple
    // stored at the given index within `RTUPLES_LIST().array`.
    constexpr bool get_plus(int) const;
    // Returns whether the chirotope evaluates to `-` on the `R`-tuple
    // stored at the given index within `RTUPLES_LIST().array`.
    constexpr bool get_minus(int) const;
    // Returns whether the chirotope is non-zero on the `R`-tuple
    // stored at the given index within `RTUPLES_LIST().array`.
    constexpr bool is_basis(int) const;
    // Evaluates the chirotope on the `R`-tuple in `RTUPLES_LIST().array`
    // at the given index. Possible return values are `'0'`, `'+'`, and
    // `'-'`. Identival to `evaluate(int)`.
    constexpr char get(int) const;
    // Evaluates the chirotope on the `R`-tuple in `RTUPLES_LIST().array`
    // at the given index. Possible return values are `'0'`, `'+'`, and
    // `'-'`. Identival to `get(int)`.
    constexpr char evaluate(int) const;
    // Sets the appropriate bit of the appropriate integer in `plus`.
    // See the documentation of `plus`.
    constexpr void set_plus(int, bool);
    // Sets the appropriate bit of the appropriate integer in `minus`.
    // See the documentation of `minus`.
    constexpr void set_minus(int, bool);
    // Makes the chirotope evaluate for the `R`-tuple in 
    // `RTUPLES_LIST().array` at the given index to `0`, `+`, or `-`,
    // if the given `char` is `'0'`, `'+'`, or `'-'` respectively.
    constexpr void set(int, char);
    // Makes the chirotope evaluate for the `R`-tuple corresponding
    // to the `r`th bit of the `i`th integer of `plus` and `minus`
    // to `0`, `+`, or `-`, if the given `char` is `'0'`, `'+'`, or 
    // `'-'` respectively.
    constexpr void set(int, int, char);

    // Returns `true` if this chirotope is constant zero.
    constexpr bool is_zero() const;
    // Returns the number of bases of the chirotope, i.e. the number
    // of distinct `R`-tuples on which it evaluates to a non-zero value.
    constexpr int countbases() const;
    // Given a string of `'+'`, `'-'`, and `'0'` characters of length
    // `NR_RTUPLES`, sets this chirotope to have the values prescribed
    // by it: the chirotope will evaluate to the value specified by the
    // `i`th character of the string at the `R`-tuple stored in
    // `RTUPLES_LIST().array[i]`.
    constexpr Chirotope& read(const std::string&);
    // Returns whether this chirotope weak maps to the other chirotope.
    // This happens if for all `R`-tuples the other chirotope being
    // non-zero implies that the two chirotopes are equal on that given
    // `R`-tuple.
    constexpr bool weak_maps_to(const Chirotope&) const;
    // Returns whether the oriented matroid defined by this chirotope 
    // weak maps to the oriented matroid defined by the other chirotope.
    // This happens if and only if this chirotope weak maps to the
    // other chirotope or this chirotope weak maps to the inverse of
    // the other chirotope (see `inverse` for details).
    constexpr bool OM_weak_maps_to(const Chirotope&) const;
    // Returns whether this chirotope is equal to the other one.
    constexpr bool operator==(const Chirotope&) const;
    // Returns whether this chirotope is not equal to the other one.
    constexpr bool operator!=(const Chirotope&) const;
    // Returns whether the oriented matroid defined by this chirotope
    // is the same as the oriented matroid defined by the other
    // chirotope. This happens if and only if they are equal, or inverses
    // of each other.
    constexpr bool is_same_OM_as(const Chirotope&) const;

    // Writes the chirotope to the given output stream as a sequence of
    // `'0'`, `'+'`, and `'-'` characters: the `i`th character is just
    // the chirotope evaluated on the `i`th `R`-tuple stored in
    // `RTUPLES_LIST().array`.
    friend std::ostream& operator<< <>(std::ostream&, const Chirotope&);

    // Writes the chirotope to the given file stream as a sequence of
    // `'0'`, `'+'`, and `'-'` characters: the `i`th character is just
    // the chirotope evaluated on the `i`th `R`-tuple stored in
    // `RTUPLES_LIST().array`.
    friend std::ofstream& operator<< <>(std::ofstream&, const Chirotope&);
};

// =============
// Matroid<R, N>
// =============

template<int R, int N>
struct Matroid;

// /////////////////////
// =====================
// CONSTEXPR DEFINITIONS
// =====================
// /////////////////////
// #learningcpp: these must be available in the same header,
//    as they will be executed at compile time

// ===============
// Chirotope<R, N>
// ===============

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
        break;
    case '-':
        minus[i] |= (uint32_t)1 << r;
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

// This file declares templates, so their implementations must
// be in this same header file as well.
#include "OMs.cpp"
