#include <cstdint>
#include <iostream>
#include <fstream>
#include "signvectors.hpp"

// ==============
//   bit_vector 
// ==============

template<int L>
constexpr bit_vector<L>::bit_vector(const std::string& from): bits {} {
    read(from);
}

template<int L>
constexpr       uint32_t& bit_vector<L>::operator[](int idx) { 
    return bits[idx]; 
}

template<int L>
constexpr const uint32_t& bit_vector<L>::operator[](int idx) const { 
    return bits[idx]; 
}

template<int L>
constexpr bool bit_vector<L>::get_bit(int i, int r) const {
    return (bits[i] >> r) & 1; 
}

template<int L>
constexpr bool bit_vector<L>::get_bit(int idx) const {
    return get_bit(idx >> 5, idx & 31);
}

template<int L>
constexpr bit_vector<L>& bit_vector<L>::set_bit(int i, int r, bool value) {
    bits[i] = bits[i] & ~((uint32_t)1 << r) | (uint32_t)value << r;
    return *this;
}

template<int L>
constexpr bit_vector<L>& bit_vector<L>::set_bit(int idx, bool value) {
    return set_bit(idx >> 5, idx & 31, value);
}

template<int L>
constexpr bit_vector<L>& bit_vector<L>::set_using_char(int i, int r, char value) {
    switch (value)
    {
    case '0':
        bits[i] &= ~((uint32_t)1 << r);
        break;
    case '1':
        bits[i] |= (uint32_t)1 << r;
        break;
    
    default:
        throw std::invalid_argument("An individual bit of a bitvector "
        "can only be set using the characters '0' or '1'. A bit was attempted"
        " to be set using the character '"+std::to_string(value)+"'.");
        break;
    }
    return *this;
}

template<int L>
constexpr bit_vector<L>& bit_vector<L>::read(const std::string& from) {
    if (from.length() != L) throw std::invalid_argument("The input string is of incorrect length! "
        "Input length: "+std::to_string(from.length())+", expected length: "+std::to_string(L)+"."
        " Input string: <"+from+">");
    for (auto i = 0; i < NR_INT32 - 1; i++) {
        for (auto r = 0; r < 32; r++) {
            set_using_char(i, r, from[i * 32 + r]);
        }
    }
    for (auto r = 0; r < NR_REMAINING_BITS; r++) {
        set_using_char(NR_INT32 - 1, r, from[(NR_INT32 - 1) * 32 + r]);
    }
    return *this;
}

template<int L>
constexpr bit_vector<L>& bit_vector<L>::invert() {
    for (auto i = 0; i < NR_INT32 - 1; i++) {
        bits[i] = ~bits[i];
    } 
    bits[NR_INT32 - 1] = ~bits[NR_INT32 - 1] 
        & (((uint32_t)1 << NR_REMAINING_BITS) - 1);
    return *this;
}

template<int L>
constexpr bit_vector<L> bit_vector<L>::inverse() const {
    bit_vector<L> ret;
    for (auto i = 0; i < NR_INT32 - 1; i++) {
        ret.bits[i] = ~bits[i];
    } 
    ret.bits[NR_INT32 - 1] = ~bits[NR_INT32 - 1] 
        & (((uint32_t)1 << NR_REMAINING_BITS) - 1);
    return ret;
}

template<int L>
constexpr bit_vector<L> bit_vector<L>::operator~() const {
    return inverse();
}

template<int L>
constexpr bit_vector<L> bit_vector<L>::operator&(const bit_vector<L>& other) const {
    bit_vector<L> ret;
    for (auto i = 0; i < NR_INT32; i++) {
        ret.bits[i] = bits[i] & other.bits[i];
    }
    return ret;
}

template<int L>
constexpr bit_vector<L> bit_vector<L>::operator|(const bit_vector<L>& other) const {
    bit_vector<L> ret;
    for (auto i = 0; i < NR_INT32; i++) {
        ret.bits[i] = bits[i] | other.bits[i];
    }
    return ret;
}

template<int L>
constexpr bit_vector<L> bit_vector<L>::operator^(const bit_vector<L>& other) const {
    bit_vector<L> ret;
    for (auto i = 0; i < NR_INT32; i++) {
        ret.bits[i] = bits[i] ^ other.bits[i];
    }
    return ret;
}

template<int L>
constexpr bool bit_vector<L>::is_zero() const {
    for (auto i = 0; i < NR_INT32; i++) {
        if (bits[i]) return false;
    }
    return true;
}

template<int L>
constexpr int bit_vector<L>::count_ones() const {
    int sum = 0;
    for (auto i = 0; i < NR_INT32; i++) {
        uint32_t n = bits[i];
        // Bit magic; n = number of 1s in bits[i] by the end
        n = (((uint32_t)0xaaaaaaaa & n) >> 1) + ((uint32_t)0x55555555 & n);
        n = (((uint32_t)0xcccccccc & n) >> 2) + ((uint32_t)0x33333333 & n);
        n = (((uint32_t)0xf0f0f0f0 & n) >> 4) + ((uint32_t)0x0f0f0f0f & n);
        n = (((uint32_t)0xff00ff00 & n) >> 8) + ((uint32_t)0x00ff00ff & n);
        n = (((uint32_t)0xffff0000 & n) >> 16) + ((uint32_t)0x0000ffff & n);
        sum += n;
    }
    return sum;
}

template<int L>
constexpr bool bit_vector<L>::bitwise_greater_than(const bit_vector<L>& other) const {
    for (auto i = 0; i < NR_INT32; i++) {
        if (~bits[i] & other.bits[i]) return false;
    }
    return true;
}

template<int L>
constexpr bool bit_vector<L>::operator==(const bit_vector<L>& other) const {
    for (auto i = 0; i < NR_INT32; i++) {
        if (bits[i] != other.bits[i]) return false;
    }
    return true;
}

template<int L>
constexpr bool bit_vector<L>::operator!=(const bit_vector<L>& other) const {
    for (auto i = 0; i < NR_INT32; i++) {
        if (bits[i] != other.bits[i]) return true;
    }
    return false;
}

template<int L>
constexpr bool bit_vector<L>::operator<=(const bit_vector<L>& other) const {
    return other.bitwise_greater_than(*this);
}

template<int L>
constexpr bool bit_vector<L>::operator>=(const bit_vector<L>& other) const {
    return bitwise_greater_than(other);
}

template<int L>
std::ostream& operator<<(std::ostream& os, const bit_vector<L>& v) {
    for (auto i = 0; i < v.NR_INT32 - 1; i++) {
        for (auto r = 0; r < 32; r++) {
            if (v.bits[i] & ((uint32_t)1 << r))
                os << '1';
            else
                os << '0';
        }
    }
    for (auto r = 0; r < v.NR_REMAINING_BITS; r++) {
        if (v.bits[v.NR_INT32 - 1] & ((uint32_t)1 << r))
            os << '1';
        else
            os << '0';
    }
    return os;
}

template<int L>
std::ofstream& operator<<(std::ofstream& of, const bit_vector<L>& v) {
    for (auto i = 0; i < v.NR_INT32 - 1; i++) {
        for (auto r = 0; r < 32; r++) {
            if (v.bits[i] & ((uint32_t)1 << r))
                of << '1';
            else
                of << '0';
        }
    }
    for (auto r = 0; r < v.NR_REMAINING_BITS; r++) {
        if (v.bits[v.NR_INT32 - 1] & ((uint32_t)1 << r))
            of << '1';
        else
            of << '0';
    }
    return of;
}

template<int L>
std::istream& operator>>(std::istream& is, bit_vector<L>& v) {
    std::string str;
    is >> str;
    v.read(str);
    return is;
}

template<int L>
std::ifstream& operator>>(std::ifstream& ifs, bit_vector<L>& v) {
    std::string str;
    ifs >> str;
    v.read(str);
    return ifs;
}

// ===============
//   sign_vector
// ===============

template<int L>
constexpr sign_vector<L>::sign_vector(const std::string& from): plus(), minus() {
    read(from);
}

template<int L>
constexpr bool sign_vector<L>::is_zero(int idx) const {
    return !is_nonzero(idx);
}

template<int L>
constexpr bool sign_vector<L>::is_nonzero(int idx) const {
    return (plus[idx >> 5] | minus[idx >> 5]) & ((uint32_t)1 << (idx & 31));
}

template<int L>
constexpr char sign_vector<L>::get_char(int i, int r) const {
    if (plus.get_bit(i, r)) return '+';
    else if (minus.get_bit(i, r)) return '-';
    else return '0';
}

template<int L>
constexpr char sign_vector<L>::get_char(int idx) const {
    return get_char(idx >> 5, idx & 31);
}

template<int L>
constexpr sign_vector<L>& sign_vector<L>::set_sign(int i, int r, char value) {
    switch (value)
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
        "character was '"+std::to_string(value)+"'.");
        break;
    }
	return *this;
}

template<int L>
constexpr sign_vector<L>& sign_vector<L>::set_sign(int idx, char value) {
    set_sign(idx >> 5, idx & 31, value);
    return *this;
}

template<int L>
constexpr sign_vector<L>& sign_vector<L>::read(const std::string& str) {
    if (str.length() != L) throw std::invalid_argument("The input string is of incorrect length! "
        "Input length: "+std::to_string(str.length())+", expected length: "+std::to_string(L)+"."
        " Input string: <"+str+">");
    for (auto i = 0; i < NR_INT32 - 1; i++) {
        for (auto c = 0; c < 32; c++) {
            set_sign(i, c, str[i * 32 + c]);
        }
    }
    for (auto c = 0; c < NR_REMAINING_BITS; c++) {
        set_sign(NR_INT32 - 1, c, str[(NR_INT32 - 1) * 32 + c]);
    }
    return *this;
}

template<int L>
constexpr sign_vector<L>& sign_vector<L>::invert() {
    for (auto i = 0; i < NR_INT32; i++) {
        uint32_t temp = plus[i];
        plus[i] = minus[i];
        minus[i] = temp;
    }
    return *this;
}

template<int L>
constexpr sign_vector<L> sign_vector<L>::inverse() const {
    sign_vector<L> ret;
    for (auto i = 0; i < NR_INT32; i++) {
        ret.plus[i] = minus[i];
        ret.minus[i] = plus[i];
    }
    return ret;
}

template<int L>
constexpr sign_vector<L> sign_vector<L>::operator~() const {
    return inverse();
}

template<int L>
constexpr sign_vector<L>& sign_vector<L>::compose(
    const sign_vector<L>& other
) {
    for (auto i = 0; i < NR_INT32; i++) {
        plus[i] |= ~(plus[i] | minus[i]) & other.plus[i];
        minus[i] |= ~(plus[i] | minus[i]) & other.minus[i];
    }
    return *this;
}

template<int L>
constexpr sign_vector<L> sign_vector<L>::composed(
    const sign_vector<L>& other
) const {
    sign_vector<L> ret;
    for(auto i = 0; i < NR_INT32; i++) {
        ret.plus[i] = plus[i] |
            (~(plus[i] | minus[i]) & other.plus[i]);
        ret.minus[i] = minus[i] |
            (~(plus[i] | minus[i]) & other.minus[i]);
    }
    return ret;
}

template<int L>
constexpr sign_vector<L> sign_vector<L>::operator*(
    const sign_vector<L>& other
) const {
    return composed(other);
}

template<int L>
constexpr bool sign_vector<L>::is_zero() const {
    return plus.is_zero() && minus.is_zero();
}

template<int L>
constexpr bool sign_vector<L>::signwise_greater_than(
    const sign_vector<L>& other
) const {
    return plus.bitwise_greater_than(other.plus) 
        && minus.bitwise_greater_than(other.minus);
}

template<int L>
constexpr bool sign_vector<L>::operator==(
    const sign_vector<L>& other
) const {
    return (plus == other.plus) && (minus == other.minus);
}

template<int L>
constexpr bool sign_vector<L>::operator!=(
    const sign_vector<L>& other
) const {
    return (plus != other.plus) || (minus != other.minus);
}

template<int L>
constexpr bool sign_vector<L>::operator<=(
    const sign_vector<L>& other
) const {
    return other.signwise_greater_than(*this);
}

template<int L>
constexpr bool sign_vector<L>::operator>=(
    const sign_vector<L>& other
) const {
    return signwise_greater_than(other);
}

template<int L>
std::ostream& operator<<(std::ostream& os, const sign_vector<L>& v) {
    for (auto i = 0; i < v.NR_INT32 - 1; i++) {
        for (auto c = 0; c < 32; c++) {
            os << v.get_char(i, c);
        }
    }
    for (auto c = 0; c < v.NR_REMAINING_BITS; c++) {
        os << v.get_char(v.NR_INT32 - 1, c);
    }
    return os;
}

template<int L>
std::ofstream& operator<<(std::ofstream& of, const sign_vector<L>& v) {
    for (auto i = 0; i < v.NR_INT32 - 1; i++) {
        for (auto c = 0; c < 32; c++) {
            of << v.get_char(i, c);
        }
    }
    for (auto c = 0; c < v.NR_REMAINING_BITS; c++) {
        of << v.get_char(v.NR_INT32 - 1, c);
    }
    return of;
}

template<int L>
std::istream& operator>>(std::istream& is, sign_vector<L>& v) {
    std::string str;
    is >> str;
    v.read(str);
    return is;
}

template<int L>
std::ifstream& operator>>(std::ifstream& ifs, sign_vector<L>& v) {
    std::string str;
    ifs >> str;
    v.read(str);
    return ifs;
}
