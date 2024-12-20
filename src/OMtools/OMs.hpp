#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <array>
#include <vector>
#include <utility>
#include <cstdint>
#include "mymath.hpp"
#include "NchooseK.hpp"
#include "signvectors.hpp"

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

// Represents a matroid of rank `R` on `N` elements.
template<int R, int N>
struct Matroid: public bit_vector<binomial_coefficient(N, R)> {
    // =============
    //   CONSTANTS
    // =============

    // Use `RTUPLES::LIST::array[i][k]` to access the `k`th element of 
    // the `i`th `R`-element subset of `0..N-1`. `RTUPLES::LIST::array` 
    // evaluates at compile-time.
    using RTUPLES = Rtuples::RTUPLES<char, R, N, int>;
    // The base class of `Matroid<R,N>` is `bit_vector<binomial_coefficient(N, R)>`
    using BASE = bit_vector<RTUPLES::NR>;

    // ================
    //   CONSTRUCTORS
    // ================

    // Initializes `char_vector` to an array of 0s. Note that this
    // does not satisfy the matroid axioms, so it is not a valid
    // matroid.
    constexpr Matroid(): BASE{} {}
    // Constructs a matroid from a string, by reading a 0-1 characteristic
    // vector, analogously to `read()`.
    constexpr Matroid(const std::string& from): BASE{} 
    { read(from); }

    constexpr Matroid(const BASE& other): BASE(other) {} 

    // ===============================
    //   WRAPPED ACCESS TO VARIABLES
    // ===============================
    
    // Returns whether the `R`-tuple stored at the given index within 
    // `RTUPLES::LIST::array` is a basis of this matroid.
    constexpr bool is_basis(int idx) const 
    { return BASE::get_bit(idx); }
    // Change the matroid such that the `R`-tuple stored at the given
    // index withing `RTUPLES::LIST::array` is a basis if and only if
    // the `bool` argument is `true`.
    constexpr Matroid& set_basis(int idx, bool value) 
    { BASE::set_bit(idx, value); return *this; }
    // Given a string of `'0'` and `'1'` characters of length `RTUPLES::NR`,
    // sets this matroid to be the matroid prescribed by it: this matroid
    // will have the `R`-tuple stored in RTUPLES::LIST::array[i]` as a basis
    // if and only if the `i`th character of the given string is `'1'`.
    constexpr Matroid& read(const std::string& from) 
    { BASE::read(from); return *this; }

    // ===================
    //   COMPLEX QUERIES
    // ===================

    // Returns `true` if this matroid has no bases.
    constexpr bool is_zero() const 
    { return BASE::is_zero(); }
    // Return the number of bases of the matroid.
    constexpr int countbases() const 
    { return BASE::count_ones(); }
    // Returns whether the given element is a loop in this matroid.
    constexpr bool is_loop(int) const;
    // Returns whether the given element is a coloop in this matroid.
    constexpr bool is_coloop(int) const;
    // Returns whether this matroid is loopfree.
    constexpr bool is_loopfree() const;
    // Decides whether a subset of elements is independent. The empty
    // set is considered independent.
    //
    // The type `Iterable` should be iterable with a range-based for 
    // loop, any number of times, yielding `char`s.
    template<typename Iterable>
    constexpr bool is_independent(const Iterable&) const;
    // Compute the rank of a subset of elements.
    //
    // The type `Iterable` should be iterable with a range-based for 
    // loop, any number of times, yielding `char`s. If the variable
    // `max_rank` is specified, than this value will be returned
    // immediately once it is clear that this subset has at least this
    // rank.
    template<typename Iterable>
    constexpr int rank(const Iterable&, int max_rank = R) const;
    // Compute the lexicographically first maximal independent subset 
    // of the given list of elements.
    //
    // The type `Iterable` should be iterable with a range-based for 
    // loop, any number of times, yielding `char`s. If the variable
    // `max_rank` is specified, than the lexicographically first
    // independent subset of this size will be returned instead, if 
    // this exists.
    template<typename Iterable>
    constexpr std::vector<char> maximal_independent_subset(const Iterable&, int max_rank=R) const;
    // Same as `maximal_independent_subset(iterable, rank)`, but returns
    // a partially filled array of the given size, if it fails to find
    // an independent subset of such a size. The remaining entries will
    // be filled with 0.
    template<int rank, typename Iterable>
    constexpr std::array<char, rank> maximal_independent_subset_of_rank(const Iterable&) const;
    // Returns whether this matroid weak maps to the other matroid. This 
    // holds if and only if every basis of the other matroid is also a basis
    // of this matroid.
    constexpr bool weak_maps_to(const Matroid& other) const
    { return BASE::bitwise_greater_than(other); }
    // Returns whether this matroid is equal to the other one.
    constexpr bool operator==(const Matroid& other) const 
    { return BASE::operator==(other); }
    // Returns whether this matroid is not equal to the other one.
    constexpr bool operator!=(const Matroid& other) const
    { return BASE::operator!=(other); }
    // Returns whether the other matroid weak maps to this one; see `weak_maps_to()`.
    constexpr bool operator<=(const Matroid& other) const 
    { return BASE::operator<=(other); }
    // Returns whether this matroid weak maps to the other one; see `weak_maps_to()`.
    constexpr bool operator>=(const Matroid& other) const
    { return BASE::operator>=(other); }

    // ========================
    //   WRITING AND PRINTING
    // ========================

    // Writes the matroid to the given output stream as a sequence of
    // `'0'` and `'1'` characters: the `i`th character is `'1'` if and
    // only if the `i`th `R`-tuple stored in `RTUPLES::LIST::array` is
    // a basis of this matroid.
    friend std::ostream& operator<< <>(std::ostream&, const Matroid&);
    // Writes the matroid to the given output stream as a sequence of
    // `'0'` and `'1'` characters: the `i`th character is `'1'` if and
    // only if the `i`th `R`-tuple stored in `RTUPLES::LIST::array` is
    // a basis of this matroid.
    friend std::ofstream& operator<< <>(std::ofstream&, const Matroid&);
    // After discarding all leading whitespace, read the next `RTUPLES::NR`
    // characters into the given matroid, as specified in `Matroid::read(...)`.
    friend std::istream& operator>> <>(std::istream&, Matroid&);
    // After discarding all leading whitespace, read the next `RTUPLES::NR`
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

// Represents a chirotope of rank `R` on `N` elements.
template<int R, int N>
struct Chirotope: public sign_vector<binomial_coefficient(N, R)> {
    // =============
    //   CONSTANTS
    // =============

    // Use `RTUPLES::LIST::array[i][k]` to access the `k`th element of 
    // the `i`th `R`-element subset of `0..N-1`. `RTUPLES::LIST::array` 
    // evaluates at compile-time.
    using RTUPLES = Rtuples::RTUPLES<char, R, N, int>;
    // The base class of `Chirotope<R,N>` is `sign_vector<binomial_coefficient(N, R)>`
    using BASE = sign_vector<RTUPLES::NR>;

    // ================
    //   CONSTRUCTORS
    // ================

    // Initializes `plus` and `minus` to arrays of 0s. Note that this
    // does not satisfy chirotope axiom (B0), so it is not a valid
    // chirotope.
    constexpr Chirotope(): BASE{} {}
    // Constructs a chirotope from a string, analogously to `read()`.
    constexpr Chirotope(const std::string& from): BASE{}
    { read(from); }

    constexpr Chirotope(const BASE& from): BASE(from) {}

    // =====================================
    //   CONSTRUCT NEW (ORIENTED) MATROIDS
    // =====================================

    // Return the inverse of this chirotope, i.e. if this chirotope is
    // denoted by `c`, then `mc`, for which 
    // - `c[i] == '0'` implies `mc[i] == '0'`,
    // - `c[i] == '+'` implies `mc[i] == '-'`, and
    // - `c[i] == '-'` implies `mc[i] == '+'`. 
    //
    // See also `invert()` for an in-place version of this operation.
    constexpr Chirotope inverse() const
    { return Chirotope(BASE::inverse()); }
    // Changes this chirotope to be its inverse, i.e. if the original
    // version of this chirotope is denoted by `c` and the post-operation
    // state by `mc` then
    // - `c[i] == '0'` implies `mc[i] == '0'`,
    // - `c[i] == '+'` implies `mc[i] == '-'`, and
    // - `c[i] == '-'` implies `mc[i] == '+'`.
    //
    // See also `inverse()` for a non-mutating version of this operation.
    constexpr Chirotope& invert()
    { BASE::invert(); return *this;}
    // Return the underlying matroid of this chirotope, i.e. the set of
    // bases of this chirotope.
    constexpr Matroid<R, N> underlying_matroid() const;

    // ===============================
    //   WRAPPED ACCESS TO VARIABLES
    // ===============================

    // Returns whether the chirotope evaluates to `+` on the `R`-tuple
    // stored at the given index within `RTUPLES::LIST::array`.
    constexpr bool get_plus(int idx) const 
    { return BASE::plus.get_bit(idx); }
    // Returns whether the chirotope evaluates to `-` on the `R`-tuple
    // stored at the given index within `RTUPLES::LIST::array`.
    constexpr bool get_minus(int idx) const
    { return BASE::minus.get_bit(idx); }
    // Returns whether the chirotope is non-zero on the `R`-tuple
    // stored at the given index within `RTUPLES::LIST::array`.
    constexpr bool is_basis(int idx) const
    { return BASE::is_nonzero(idx); }
    // Evaluates the chirotope on the `R`-tuple in `RTUPLES::LIST::array`
    // at the given index. Possible return values are `'0'`, `'+'`, and
    // `'-'`. Identival to `get(int)`.
    constexpr char evaluate(int idx) const
    { return BASE::get_char(idx); }
    // Sets the appropriate bit of the appropriate integer in `plus`.
    // See the documentation of `plus`.
    constexpr Chirotope& set_plus(int idx, bool value)
    { BASE::plus.set_bit(idx, value); return *this; }
    // Sets the appropriate bit of the appropriate integer in `minus`.
    // See the documentation of `minus`.
    constexpr Chirotope& set_minus(int idx, bool value)
    { BASE::minus.set_bit(idx, value); return *this; }
    // Makes the chirotope evaluate for the `R`-tuple in 
    // `RTUPLES::LIST::array` at the given index to `0`, `+`, or `-`,
    // if the given `char` is `'0'`, `'+'`, or `'-'` respectively.
    constexpr Chirotope& set(int idx, char c)
    { BASE::set_sign(idx, c); return *this;}
    // Makes the chirotope evaluate for the `R`-tuple corresponding
    // to the `r`th bit of the `i`th integer of `plus` and `minus`
    // to `0`, `+`, or `-`, if the given `char` is `'0'`, `'+'`, or 
    // `'-'` respectively.
    constexpr Chirotope& set(int i, int r, char c)
    { BASE::set_sign(i, r, c); return *this; }
    // Given a string of `'+'`, `'-'`, and `'0'` characters of length
    // `RTUPLES::NR`, sets this chirotope to have the values prescribed
    // by it: the chirotope will evaluate to the value specified by the
    // `i`th character of the string at the `R`-tuple stored in
    // `RTUPLES::LIST::array[i]`.
    constexpr Chirotope& read(const std::string& from)
    { BASE::read(from); return *this; }
    
    // ===================
    //   COMPLEX QUERIES
    // ===================

    // Returns `true` if this chirotope is constant zero.
    constexpr bool is_zero() const
    { return BASE::is_zero(); }
    // Returns the number of bases of the chirotope, i.e. the number
    // of distinct `R`-tuples on which it evaluates to a non-zero value.
    constexpr int countbases() const
    { return BASE::count_nonzero(); }
    // Returns whether the given element is a loop in this chirotope.
    constexpr bool is_loop(int) const;
    // Returns whether the given element is a coloop in this chirotope.
    constexpr bool is_coloop(int) const;
    // Returns whether this chirotope is loopfree.
    constexpr bool is_loopfree() const;
    // Decides whether a subset of elements is independent.
    //
    // The type `Iterable` should be iterable with a range-based for 
    // loop, any number of times, yielding `char`s.
    template<typename Iterable>
    constexpr bool is_independent(const Iterable&) const;
    // Compute the rank of a subset of elements.
    //
    // The type `Iterable` should be iterable with a range-based for 
    // loop, any number of times, yielding `char`s. If the variable
    // `max_rank` is specified, than this value will be returned
    // immediately once it is clear that this subset has at least this
    // rank.
    template<typename Iterable>
    constexpr int rank(const Iterable&, int max_rank=R) const;
    // Compute the lexicographically first maximal independent subset 
    // of the given list of elements.
    //
    // The type `Iterable` should be iterable with a range-based for 
    // loop, any number of times, yielding `char`s. If the variable
    // `max_rank` is specified, than the lexicographically first
    // independent subset of this size will be returned instead, if 
    // this exists.
    template<typename Iterable>
    constexpr std::vector<char> maximal_independent_subset(const Iterable&, int max_rank=R) const;
    // Same as `maximal_independent_subset(iterable, rank)`, but returns
    // a partially filled array of the given size, if it fails to find
    // an independent subset of such a size. The remaining entries will
    // be filled with 0.
    template<int rank, typename Iterable>
    constexpr std::array<char, rank> maximal_independent_subset_of_rank(const Iterable&) const;
    // Returns whether this chirotope weak maps to the other chirotope.
    // This happens if for all `R`-tuples the other chirotope being
    // non-zero implies that the two chirotopes are equal on that given
    // `R`-tuple.
    constexpr bool weak_maps_to(const Chirotope& other) const
    { return BASE::signwise_greater_than(other); }
    // Return whether the underlying matroid of this chirotope weak maps
    // to the given matroid; this happens if every basis of the matroid
    // is also a basis of this chirotope.
    constexpr bool weak_maps_to(const Matroid<R, N>& matroid) const;
    // Returns whether the oriented matroid defined by this chirotope 
    // weak maps to the oriented matroid defined by the other chirotope.
    // This happens if and only if this chirotope weak maps to the
    // other chirotope or this chirotope weak maps to the inverse of
    // the other chirotope (see `inverse` for details).
    constexpr bool OM_weak_maps_to(const Chirotope& other) const;
    // Returns whether this chirotope is equal to the other one.
    constexpr bool operator==(const Chirotope& other) const
    { return BASE::operator==(other); }
    // Returns whether this chirotope is not equal to the other one.
    constexpr bool operator!=(const Chirotope& other) const
    { return BASE::operator!=(other); }
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

    // ========================
    //   WRITING AND PRINTING
    // ========================

    // Writes the chirotope to the given output stream as a sequence of
    // `'0'`, `'+'`, and `'-'` characters: the `i`th character is just
    // the chirotope evaluated on the `i`th `R`-tuple stored in
    // `RTUPLES::LIST::array`.
    friend std::ostream& operator<< <>(std::ostream&, const Chirotope&);
    // Writes the chirotope to the given file stream as a sequence of
    // `'0'`, `'+'`, and `'-'` characters: the `i`th character is just
    // the chirotope evaluated on the `i`th `R`-tuple stored in
    // `RTUPLES::LIST::array`.
    friend std::ofstream& operator<< <>(std::ofstream&, const Chirotope&);
    // After discarding all leading whitespace, read the next `RTUPLES::NR`
    // characters into the given chirotope, as specified in `Chirotope::read(...)`.
    friend std::istream& operator>> <>(std::istream&, Chirotope&);
    // After discarding all leading whitespace, read the next `RTUPLES::NR`
    // characters into the given chirotope, as specified in `Chirotope::read(...)`.
    friend std::ifstream& operator>> <>(std::ifstream&, Chirotope&);
};

// This file declares templates, so their implementations must
// be in this same header file as well.
#include "OMs_impl.hpp"
