#pragma once
#include <cstdint>
#include <string>
#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <bit>
#include "mymath.hpp"

// ==============
//   bit_vector
// ==============

// Represents a 0-1 sequence of length `L`.
template<int L>
struct bit_vector;

template<int L>
std::ostream& operator<<(std::ostream&, const bit_vector<L>&);
template<int L>
std::ofstream& operator<<(std::ofstream&, const bit_vector<L>&);
template<int L>
std::istream& operator>>(std::istream&, bit_vector<L>&);
template<int L>
std::ifstream& operator>>(std::ifstream&, bit_vector<L>&);

template<int L>
struct bit_vector {
    // =============
    //   CONSTANTS
    // =============

    // The number of bits stored in a bitvector.
    constexpr static const int LENGTH = L;
    // Number of (unsigned) 32-bit integers needed to store all
    // the bits in a bitvector.
    constexpr static const int NR_INT32 = division_rounded_up(L, 32);
    // Number of bits used in the last 32-bit unsigned integer used to 
    // store the bitvector.
    constexpr static const int NR_REMAINING_BITS = L - 32 * (NR_INT32 - 1);

    // =============
    //   VARIABLES
    // =============

    // The `r`th least significant bit of `bits[i]` is 1 if and
    // only if the `i * 32 + r`th bit of this bitvector is 1.
    uint32_t bits[NR_INT32];

    // ================
    //   CONSTRUCTORS
    // ================

    // Creates a bitvector that only contains 0s.
    constexpr bit_vector(): bits{} {}
    // Create a bitvector by copying the contents of an array.
    constexpr bit_vector(const std::array<uint32_t,NR_INT32>& from);
    // Reads in a string of length `L` which only contains the
    // characters `'0'` and `'1'`, and constructs a bitvector
    // corresponding to it: the `idx`th bit of the bitvector will
    // be `1` if and only if the `idx`th character in the string
    // was `'1'`. See also `read(...)`.
    constexpr bit_vector(const std::string& from);

    // ===============================
    //   WRAPPED ACCESS TO VARIABLES
    // ===============================

    // Gives mutable access to `bits[idx]`.
    constexpr       uint32_t& operator[](int idx);
    // Gives immutable access to `bits[idx]`.
    constexpr const uint32_t& operator[](int idx) const;
    // Returns `true` if and only if the `i * 32 + r`th 
    // bit of this bitvector is `1`.
    constexpr bool get_bit(int i, int r) const;
    // Returns `true` if and only if the `idx`th bit of this
    // bitvector is `1`.
    constexpr bool get_bit(int idx) const;
    // Set the `i * 32 + r`th bit of this bitvector to `1`
    // if `value` is true and `0` if `value` is false.
    constexpr bit_vector& set_bit(int i, int r, bool value);
    // Set the `idx`th bit of this bitvector to `1` if `value`
    // is true and `0` if `value` is false.
    constexpr bit_vector& set_bit(int idx, bool value);
    // Set the `i * 32 + r`th bit of this bitvector to `1` if
    // `value` is `'1'` and `0` if `value` is `'0'`.
    constexpr bit_vector& set_using_char(int i, int r, char value);
    // Reads in a string of length `L` which only contains the
    // characters `'0'` and `'1'`, and stores the resulting 
    // bitvector: the `idx`th bit of this bitvector will
    // be `1` if and only if the `idx`th character in the string
    // was `'1'`.
    constexpr bit_vector& read(const std::string&);

    // ============================
    //   CONSTRUCT NEW BITVECTORS
    // ============================

    // Flips all bits in this bit-vector from 0 to 
    // 1 and from 1 to 0. In-place.
    constexpr bit_vector& invert();
    // Flips all bits in this bit-vector from 0 to 
    // 1 and from 1 to 0, and returns the result.
    constexpr bit_vector inverse() const;
    // Flips all bits in this bit-vector from 0 to 
    // 1 and from 1 to 0, and returns the result.
    constexpr bit_vector operator~() const;

    // Construct the bitwise AND of the two bitvectors.
    constexpr bit_vector operator&(const bit_vector&) const;
    // Construct the bitwise OR of the two bitvectors.
    constexpr bit_vector operator|(const bit_vector&) const;
    // Construct the bitwise XOR of the two bitvectors.
    constexpr bit_vector operator^(const bit_vector&) const;

    // ===================
    //   COMPLEX QUERIES
    // ===================

    // Returns whether all bits in this bitvector are zero.
    constexpr bool is_zero() const;
    // Returns the number of 1 bits in this bitvector.
    constexpr int count_ones() const;
    // Returns the list of indices, in increasing order, at which
    // a `1` is present in the bit-vector.
    std::vector<int> indices_with_ones() const;
    // Returns true if this bitvector contains a 1 wherever
    // the other bitvector also contains a 1.
    constexpr bool bitwise_greater_than(const bit_vector&) const;
    // Returns true if this bitvector is the same as the other
    // one as 0-1 sequences.
    constexpr bool operator==(const bit_vector&) const;
    // Returns true if this bitvector is NOT the same as the
    // other one as 0-1 sequences.
    constexpr bool operator!=(const bit_vector&) const;
    // Returns true if the other bitvector contains a 1
    // wherever this bitvector also contains a 1.
    constexpr bool operator<=(const bit_vector& other) const;
    // Returns true if this bitvector contains a 1 wherever
    // the other bitvector also contains a 1.
    constexpr bool operator>=(const bit_vector& other) const;

    // ========================
    //   WRITING AND PRINTING
    // ========================

    // Writes the bitvector to the given output stream as a sequence of
    // `'0'` and `'1'` characters: the `i`th character is `'1'` if and
    // only if the `i`th bit of this bitvector is `1`.
    friend std::ostream& operator<< <>(std::ostream&, const bit_vector&);
    // Writes the bitvector to the given output stream as a sequence of
    // `'0'` and `'1'` characters: the `i`th character is `'1'` if and
    // only if the `i`th bit of this bitvector is `1`.
    friend std::ofstream& operator<< <>(std::ofstream&, const bit_vector&);
    // After discarding all leading whitespace, read the next `L`
    // characters into the given bitvector, as specified in 
    // `bit_vector::read(...)`.
    friend std::istream& operator>> <>(std::istream&, bit_vector&);
    // After discarding all leading whitespace, read the next `L`
    // characters into the given bitvector, as specified in 
    // `bit_vector::read(...)`.
    friend std::ifstream& operator>> <>(std::ifstream&, bit_vector&);
};

// ===============
//   sign_vector
// ===============

// Represents a sequence with elements 0, +, and -.
template<int L>
struct sign_vector;

template<int L>
std::ostream& operator<<(std::ostream&, const sign_vector<L>&);
template<int L>
std::ofstream& operator<<(std::ofstream&, const sign_vector<L>&);
template<int L>
std::istream& operator>>(std::istream&, sign_vector<L>&);
template<int L>
std::ifstream& operator>>(std::ifstream&, sign_vector<L>&);

template<int L>
struct sign_vector {
    // =============
    //   CONSTANTS
    // =============

    // The number of signs stored in a signvector.
    constexpr static const int LENGTH = L;
    // Number of (unsigned) 32-bit integers needed to store as many
    // bits as there are signs in a signvector.
    constexpr static const int NR_INT32 = division_rounded_up(L, 32);
    // Grouping the signs in 32-long chunks, this many signs are left
    // for the last group (0 < n <= 32).
    constexpr static const int NR_REMAINING_BITS = L - 32 * (NR_INT32 - 1);

    // =============
    //   VARIABLES
    // =============

    // The characteristic (bit)vector of the `+`-es in the signvector.
    // This stores a `1` wherever this signvector has a `+`, and a `0`
    // everywhere else.
    bit_vector<L> plus;
    // The characteristic (bit)vector of the `-`-es in the signvector.
    // This stores a `1` wherever this signvector has a `-`, and a `0`
    // everywhere else.
    bit_vector<L> minus;

    // ================
    //   CONSTRUCTORS
    // ================

    // Creates a signvector that only contains 0s.
    constexpr sign_vector(): plus(), minus() {}
    // Construct a new signvector given the characteristic vectors
    // of `+`s and `-`s in it.
    constexpr sign_vector(
        const bit_vector<L>& p, 
        const bit_vector<L>& m
    ): plus(p), minus(m) {}
    // Reads in a string of length `L` which only contains the
    // characters `'-'`, `'0'`, and `'1'`, and constructs a signvector
    // corresponding to it: the `idx`th sign in the signvector will
    // be `+` if and only if the `idx`th character in the string
    // was `'+'`, and similarly for `-` and `0`. See also `read(...)`.
    constexpr sign_vector(const std::string& from);

    // ===============================
    //   WRAPPED ACCESS TO VARIABLES
    // ===============================

    // Returns `true` if and only if the `idx`th sign of this
    // signvector is `0`.
    constexpr bool is_zero(int idx) const;
    // Returns `true` if and only if the `idx`th sign of this
    // signvector is `-` or `+`.
    constexpr bool is_nonzero(int idx) const;
    // Return the `i * 32 + r`th sign of this signvector as a 
    // character: `'-'` for `-`, `'0'` for `0` and `'+'` for `+`.
    constexpr char get_char(int i, int r) const;
    // Return the `idx`th sign of this signvector as a character:
    // `'-'` for `-`, `'0'` for `0` and `'+'` for `+`.
    constexpr char get_char(int idx) const;
    // Set the `i * 32 + r`th sign of this signvector to the
    // sign represented by the given character `value`.
    constexpr sign_vector& set_sign(int i, int r, char value);
    // Set the `idx`th sign of this signvector to the
    // sign represented by the given character `value`.
    constexpr sign_vector& set_sign(int idx, char value);
    // Reads in a string of length `L` which only contains the
    // characters `'-'`, `'0'`, and `'+'`, and stores the resulting 
    // signvector: the `idx`th sign of this signvector will
    // be `+` if and only if the `idx`th character in the string
    // was `'+'`, and similarly for `0` and `-`.
    constexpr sign_vector& read(const std::string&);

    // =============================
    //   CONSTRUCT NEW SIGNVECTORS
    // =============================

    // Turn all `+`s to `-`s and `-`s to `+`s in this
    // signvector. In-place operation.
    constexpr sign_vector& invert();
    // Turn all `+`s to `-`s and `-`s to `+`s in this
    // signvector, and return the result.
    constexpr sign_vector inverse() const;
    // Turn all `+`s to `-`s and `-`s to `+`s in this
    // signvector, and return the result.
    constexpr sign_vector operator~() const;
    
    // Replace this signvector with its composition
    // with the other one: the composition `c` 
    // equals this signvector wherever this one is non-zero, 
    // and the other signvector at all other places. 
    constexpr sign_vector& compose(const sign_vector&);
    // Return the composition of this sign vector with the
    // other one as a new signvector: the composition `c` 
    // equals this signvector wherever this one is non-zero, 
    // and the other signvector at all other places. 
    constexpr sign_vector composed(const sign_vector&) const;
    // Return the composition of this sign vector with the
    // other one as a new signvector: the composition `c` 
    // equals this signvector wherever this one is non-zero, 
    // and the other signvector at all other places. 
    constexpr sign_vector operator*(const sign_vector&) const;

    // ===================
    //   COMPLEX QUERIES
    // ===================

    // Returns whether all signs in this signvector are zero.
    constexpr bool is_zero() const;
    // Returns the number of nonzero signs in this signvector.
    constexpr int count_nonzero() const;
    // Returns true if wherever the other signvector is nonzero,
    // the two signvectors agree.
    constexpr bool signwise_greater_than(const sign_vector&) const;
    // Returns true if this signvector is the same as the other
    // one as sign-sequences.
    constexpr bool operator==(const sign_vector&) const;
    // Returns true if this signvector is NOT the same as the other
    // one as sign-sequences.
    constexpr bool operator!=(const sign_vector&) const;
    // Returns true if wherever this signvector is nonzero,
    // the two signvectors agree.
    constexpr bool operator<=(const sign_vector&) const;
    // Returns true if wherever the other signvector is nonzero,
    // the two signvectors agree.
    constexpr bool operator>=(const sign_vector&) const;

    // ========================
    //   WRITING AND PRINTING
    // ========================

    // Writes the signvector to the given output stream as a sequence of
    // `'-'`, `'0'`, and `'+'` characters: the `idx`th character is `'+'` 
    // if and only if the `idx`th sign of this signvector is `+`, and
    // similarly for `-` and `0`.
    friend std::ostream& operator<< <>(std::ostream&, const sign_vector&);
    // Writes the signvector to the given output stream as a sequence of
    // `'-'`, `'0'`, and `'+'` characters: the `idx`th character is `'+'` 
    // if and only if the `idx`th sign of this signvector is `+`, and
    // similarly for `-` and `0`.
    friend std::ofstream& operator<< <>(std::ofstream&, const sign_vector&);
    // After discarding all leading whitespace, read the next `L`
    // characters into the given signvector, as specified in 
    // `sign_vector::read(...)`.
    friend std::istream& operator>> <>(std::istream&, sign_vector&);
    // After discarding all leading whitespace, read the next `L`
    // characters into the given signvector, as specified in 
    // `sign_vector::read(...)`.
    friend std::ifstream& operator>> <>(std::ifstream&, sign_vector&);
};

#include "signvectors.cpp"