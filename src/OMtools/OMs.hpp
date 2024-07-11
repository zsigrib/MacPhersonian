#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <cstdint>
#include "mymath.hpp"

//writebits? -> print an integer as a 0-1 sequence

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
    constexpr Chirotope(): plus{}, minus{} {}

    // Constructs a chirotope from a string, analogously to `read()`.
    constexpr Chirotope(std::string);

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

// This file declares templates, so their implementations must
// be in this same header file as well.
#include "OMs.cpp"
