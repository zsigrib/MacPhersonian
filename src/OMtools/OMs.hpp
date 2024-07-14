#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <cstdint>
#include "mymath.hpp"

// =======
// HELPERS 
// =======

// Counts the number of 1s in the base 2 expansion of
// the given integer; works for 32-bit ints.
constexpr int count_1bits(uint32_t x);

// Sorts an array of 1-byte signed integers, and returns the 
// of the permutation applied.
template<int R, int N>
constexpr int sort_array(std::array<char, R>&);

// Returns the index of the given `R`-tuple inside
// `Chirotope<R, N>::RTUPLES_LIST::array`.
//
// TODO: move into `NchooseK`
template<int R, int N>
constexpr int index_of_Rtuple(const std::array<char, R>&);

//writebits? -> print an integer as a 0-1 sequence

// =============
// Matroid<R, N>
// =============

template<int R, int N> struct Matroid;

template<int R, int N>
std::ostream& operator<<(std::ostream&, const Matroid<R, N>&);

template<int R, int N>
std::ofstream& operator<<(std::ofstream&, const Matroid<R, N>&);

template<int R, int N>
std::istream& operator>>(std::istream&, Matroid<R, N>&);

template<int R, int N>
std::ifstream& operator>>(std::ifstream&, Matroid<R, N>&);

template<int R, int N>
struct Matroid {
    // `(N choose R)=N!/(R!*(N-R)!)`, the number of distinct `R`-element 
    // subsets of `0..N-1`.
    static constexpr int NR_RTUPLES = binomial_coefficient(N, R);
    // Number of (unsigned) 32-bit integers needed to have 1 bit for each 
    // `R` element subset of `0..N-1`.
    static constexpr int NR_INT32 = division_rounded_up(NR_RTUPLES, 32);
    // Number of bits used in the last 32-bit unsigned integer used to 
    // store the characteristic vector.
    static constexpr int NR_REMAINING_BITS = NR_RTUPLES % 32;
    // Use `RTUPLES_LIST::array[i][k]` to access the `k`th element of 
    // the `i`th `R`-element subset of `0..N-1`. `RTUPLES_LIST::array` 
    // evaluates at compile-time.
    using RTUPLES_LIST = NchooseK<int, N, R>;

    // The `c`th least significant bit of `char_vector[i]` is 1 if and 
    // only if the `R`-tuple stored in `RTUPLES_LIST::array[i * 32 + c]`
    // is a basis of the matroid.
    uint32_t char_vector[NR_INT32];

    // Initializes `char_vector` to an array of 0s. Note that this
    // does not satisfy the matroid axioms, so it is not a valid
    // matroid.
    constexpr Matroid(): char_vector{} {}

    // Constructs a matroid from a string, by reading a 0-1 characteristic
    // vector, analogously to `read()`.
    constexpr Matroid(const std::string&);

    // Returns whether the `R`-tuple stored at the given index within 
    // `RTUPLES_LIST::array` is a basis of this matroid.
    constexpr bool is_basis(int) const;
    // Change the matroid such that the `R`-tuple stored at the given
    // index withing `RTUPLES_LIST::array` is a basis if and only if
    // the `bool` argument is `true`.
    constexpr void set_basis(int, bool);
    // Set the `r`th bit of the `i`th integer of `char_vector` to the given
    // value (`1` if `true`, `0` if `false`).
    constexpr void set(int, int, bool);
    // Set the `r`th bit of the `i`th integer of `char_vector` to be 0 or 1
    // depending on the `char` argument being `'0'` or `'1'`.
    constexpr void set_using_char(int, int, char);

    // Returns `true` if this matroid has no bases.
    constexpr bool is_zero() const;
    // Return the number of bases of the matroid.
    constexpr int countbases() const;
    // Given a string of `'0'` and `'1'` characters of length `NR_RTUPLES`,
    // sets this matroid to be the matroid prescribed by it: this matroid
    // will have the `R`-tuple stored in RTUPLES_LIST::array[i]` as a basis
    // if and only if the `i`th character of the given string is `'1'`.
    constexpr Matroid& read(const std::string&);
    // Returns whether this matroid weak maps to the other matroid. This 
    // holds if and only if every basis of the other matroid is also a basis
    // of this matroid.
    constexpr bool weak_maps_to(const Matroid&) const;
    // Returns whether this matroid is equal to the other one.
    constexpr bool operator==(const Matroid&) const;
    // Returns whether this matroid is not equal to the other one.
    constexpr bool operator!=(const Matroid&) const;
    // Returns whether the other matroid weak maps to this one; see `weak_maps_to()`.
    constexpr bool operator<=(const Matroid&) const;
    // Returns whether this matroid weak maps to the other one; see `weak_maps_to()`.
    constexpr bool operator>=(const Matroid&) const;

    // Writes the matroid to the given output stream as a sequence of
    // `'0'` and `'1'` characters: the `i`th character is `'1'` if and
    // only if the `i`th `R`-tuple stored in `RTUPLES_LIST::array` is
    // a basis of this matroid.
    friend std::ostream& operator<< <>(std::ostream&, const Matroid&);

    // Writes the matroid to the given output stream as a sequence of
    // `'0'` and `'1'` characters: the `i`th character is `'1'` if and
    // only if the `i`th `R`-tuple stored in `RTUPLES_LIST::array` is
    // a basis of this matroid.
    friend std::ofstream& operator<< <>(std::ofstream&, const Matroid&);

    // After discarding all leading whitespace, read the next `NR_RTUPLES`
    // characters into the given matroid, as specified in `Matroid::read(...)`.
    friend std::istream& operator>> <>(std::istream&, Matroid&);

    // After discarding all leading whitespace, read the next `NR_RTUPLES`
    // characters into the given matroid, as specified in `Matroid::read(...)`.
    friend std::ifstream& operator>> <>(std::ifstream&, Matroid&);
};

// ===============
// Chirotope<R, N>
// ===============

template<int R, int N> struct Chirotope;

template<int R, int N>
std::ostream& operator<<(std::ostream&, const Chirotope<R, N>&);

template<int R, int N>
std::ofstream& operator<<(std::ofstream&, const Chirotope<R, N>&);

template<int R, int N>
std::istream& operator>>(std::istream&, Chirotope<R, N>&);

template<int R, int N>
std::ifstream& operator>>(std::ifstream&, Chirotope<R, N>&);

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
    // Use `RTUPLES_LIST::array[i][k]` to access the `k`th element of 
    // the `i`th `R`-element subset of `0..N-1`. `RTUPLES_LIST::array` 
    // evaluates at compile-time.
    using RTUPLES_LIST = NchooseK<int, N, R>;
    
    // The `c`th least significant bit of `plus[i]` is 1 if and only if 
    // the chirotope evaluates to `+` on the `R`-tuple stored in
    // `RTUPLES_LIST::array[i * 32 + c]`.
    uint32_t plus[NR_INT32];
    // The `c`th least significant bit of `minus[i]` is 1 if and only if 
    // the chirotope evaluates to `-` on the `R`-tuple stored in
    // `RTUPLES_LIST::array[i * 32 + c]`.
    uint32_t minus[NR_INT32];

    // Initializes `plus` and `minus` to arrays of 0s. Note that this
    // does not satisfy chirotope axiom (B0), so it is not a valid
    // chirotope.
    constexpr Chirotope(): plus{}, minus{} {}

    // Constructs a chirotope from a string, analogously to `read()`.
    constexpr Chirotope(const std::string&);

    // Return the inverse of this chirotope, i.e. if this chirotope is
    // denoted by `c`, then `mc`, for which 
    // - `c[i] == '0'` implies `mc[i] == '0'`,
    // - `c[i] == '+'` implies `mc[i] == '-'`, and
    // - `c[i] == '-'` implies `mc[i] == '+'`. 
    constexpr Chirotope inverse() const;

    // Return the underlying matroid of this chirotope, i.e. the set of
    // bases of this chirotope.
    constexpr Matroid<R, N> underlying_matroid() const;
    
    // Returns whether the chirotope evaluates to `+` on the `R`-tuple
    // stored at the given index within `RTUPLES_LIST::array`.
    constexpr bool get_plus(int) const;
    // Returns whether the chirotope evaluates to `-` on the `R`-tuple
    // stored at the given index within `RTUPLES_LIST::array`.
    constexpr bool get_minus(int) const;
    // Returns whether the chirotope is non-zero on the `R`-tuple
    // stored at the given index within `RTUPLES_LIST::array`.
    constexpr bool is_basis(int) const;
    // Evaluates the chirotope on the `R`-tuple in `RTUPLES_LIST::array`
    // at the given index. Possible return values are `'0'`, `'+'`, and
    // `'-'`. Identival to `evaluate(int)`.
    constexpr char get(int) const;
    // Evaluates the chirotope on the `R`-tuple in `RTUPLES_LIST::array`
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
    // `RTUPLES_LIST::array` at the given index to `0`, `+`, or `-`,
    // if the given `char` is `'0'`, `'+'`, or `'-'` respectively.
    constexpr void set(int, char);
    // Makes the chirotope evaluate for the `R`-tuple corresponding
    // to the `r`th bit of the `i`th integer of `plus` and `minus`
    // to `0`, `+`, or `-`, if the given `char` is `'0'`, `'+'`, or 
    // `'-'` respectively.
    constexpr void set(int, int, char);
    // Set this chirotope to `0` on all non-bases of the given matroid.
    // Note that the resulting object might not satisfy the chirotope
    // axioms. 
    constexpr Chirotope& restrict_to_matroid(const Matroid<R, N>&) const;

    // Returns `true` if this chirotope is constant zero.
    constexpr bool is_zero() const;
    // Returns the number of bases of the chirotope, i.e. the number
    // of distinct `R`-tuples on which it evaluates to a non-zero value.
    constexpr int countbases() const;
    // Given a string of `'+'`, `'-'`, and `'0'` characters of length
    // `NR_RTUPLES`, sets this chirotope to have the values prescribed
    // by it: the chirotope will evaluate to the value specified by the
    // `i`th character of the string at the `R`-tuple stored in
    // `RTUPLES_LIST::array[i]`.
    constexpr Chirotope& read(const std::string&);
    // Returns whether this chirotope weak maps to the other chirotope.
    // This happens if for all `R`-tuples the other chirotope being
    // non-zero implies that the two chirotopes are equal on that given
    // `R`-tuple.
    constexpr bool weak_maps_to(const Chirotope&) const;
    // Return whether the underlying matroid of this chirotope weak maps
    // to the given matroid; this happens if every basis of the matroid
    // is also a basis of this chirotope.
    constexpr bool weak_maps_to(const Matroid<R, N>&) const;
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

    // Checks whether this chirotope satisfies the chirotope axioms.
    // The implementation is ported legacy code from Nevena.
    constexpr bool is_chirotope() const;
    // Ported legacy code from Nevena.
    constexpr bool b2prime(char sign, const std::array<char,R>&, const std::array<char,R>&) const;
    // Ported legacy code from Nevena.
    constexpr bool axB2(char sign, char s1, char s2, int in1, int in2) const;

    // Writes the chirotope to the given output stream as a sequence of
    // `'0'`, `'+'`, and `'-'` characters: the `i`th character is just
    // the chirotope evaluated on the `i`th `R`-tuple stored in
    // `RTUPLES_LIST::array`.
    friend std::ostream& operator<< <>(std::ostream&, const Chirotope&);

    // Writes the chirotope to the given file stream as a sequence of
    // `'0'`, `'+'`, and `'-'` characters: the `i`th character is just
    // the chirotope evaluated on the `i`th `R`-tuple stored in
    // `RTUPLES_LIST::array`.
    friend std::ofstream& operator<< <>(std::ofstream&, const Chirotope&);

    // After discarding all leading whitespace, read the next `NR_RTUPLES`
    // characters into the given chirotope, as specified in `Chirotope::read(...)`.
    friend std::istream& operator>> <>(std::istream&, Chirotope&);

    // After discarding all leading whitespace, read the next `NR_RTUPLES`
    // characters into the given chirotope, as specified in `Chirotope::read(...)`.
    friend std::ifstream& operator>> <>(std::ifstream&, Chirotope&);
};

// This file declares templates, so their implementations must
// be in this same header file as well.
#include "OMs.cpp"
