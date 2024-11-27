#include <cstdint>
#include <array>
#include "signvectors.hpp"
#include "signvectoroperations.hpp"

namespace sign_vector_operations {

// ==================
//   index_function
// ==================

template<int L0, int L1>
constexpr int index_function<L0, L1>::operator()(int idx) const {
    return i[idx] * 32 + r[idx];
}

template<int L0, int L1>
template<int L>
constexpr index_function<L,L1> index_function<L0,L1>::compose(
    const index_function<L,L0>& other
) const {
    index_function<L,L1> ret;
    for (auto idx = 0; idx < L; idx++) {
        ret.i[idx] = i[other.i[idx] * 32 + other.r[idx]];
        ret.r[idx] = r[other.i[idx] * 32 + other.r[idx]];
    }
    return ret;
}

template<int L0, int L1>
template<int L>
constexpr index_function<L,L1> index_function<L0,L1>::operator*(
    const index_function<L,L0>& other
) const {
    return compose(other);
}

template<int L0, int L1>
bit_vector<L1> push_forward(
    const index_function<L0,L1>& f, const bit_vector<L0>& v
) {
    auto ret = bit_vector<L1>();
    for (auto idx = 0; idx < L0; idx++) {
        ret[f.i[idx]] |= (uint32_t)v.get_bit(idx) << f.r[idx];
    }
    return ret;
}

template<int L0, int L1>
bit_vector<L0> pull_back(
    const index_function<L0,L1>& f, const bit_vector<L1>& v
) {
    bit_vector<L0> ret;
    for (auto idx = 0; idx < L0; idx++) {
        ret.set_bit(idx, v.get_bit(f.i[idx], f.r[idx]));
    }
    return ret;
}

template<int L0, int L1>
sign_vector<L1> push_forward(
    const index_function<L0,L1>& f, const sign_vector<L0>& v
) {
    return sign_vector<L1>{
        push_forward(f, v.plus), 
        push_forward(f, v.minus)
    };
}

template<int L0, int L1>
sign_vector<L0> pull_back(
    const index_function<L0,L1>& f, const sign_vector<L1>& v
) {
    return sign_vector<L0>{
        pull_back(f, v.plus),
        pull_back(f, v.minus)
    };
}

// ============================
//   Specific operation types
// ============================

template<typename T>
constexpr T Multiply<T>::applied(const T& to) const {
    T ret;
    for (auto i = 0; i < T::NR_INT32; i++) {
        ret.minus[i] = (to.minus[i] & ~minus[i]) | (to.plus[i] & minus[i]);
        ret.plus[i] = (to.plus[i] & ~minus[i]) | (to.minus[i] & minus[i]);
    }
    return ret;
}

template<typename T>
constexpr T Multiply_P0<T>::applied(const T& to) const {
    T ret;
    for (auto i = 0; i < T::NR_INT32; i++) {
        ret.minus[i] = to.minus[i] & nonzeros[i];
        ret.plus[i] = to.plus[i] & nonzeros[i];
    }
    return ret;
}

template<typename T>
constexpr T Multiply_PM0<T>::applied(const T& to) const {
    T ret;
    for (auto i = 0; i < T::NR_INT32; i++) {
        ret.minus[i] = ((to.minus[i] & ~minus[i]) | (to.plus[i] & minus[i])) & nonzeros[i];
        ret.plus[i] = ((to.plus[i] & ~minus[i]) | (to.minus[i] & minus[i])) & nonzeros[i];
    }
    return ret;
}

template<typename T0, typename T1>
constexpr T1 PushForward<T0,T1>::applied(const T0& to) const {
    return {push_forward(function, to)};
}

template<typename T0, typename T1>
constexpr T1 PullBack<T0,T1>::applied(const T0& to) const {
    return {pull_back(function, to)};
}

// ==================================
//   Composition of operation types
// ==================================

template<typename T>
constexpr Multiply<T> Compose<Multiply<T>, Multiply<T>>::operator()(
    const Multiply<T>& m1, const Multiply<T>& m0
) {
    return {m1.minus ^ m0.minus};
}

template<typename T>
constexpr Multiply_PM0<T> Compose<Multiply<T>, Multiply_P0<T>>::operator()(
    const Multiply<T>& m1, const Multiply_P0<T>& m0
) {
    return {m1.minus & m0.nonzeros, m0.nonzeros};
}

template<typename T>
constexpr Multiply_PM0<T> Compose<Multiply_P0<T>, Multiply<T>>::operator()(
    const Multiply_P0<T>& m1, const Multiply<T>& m0
) {
    return {m1.nonzeros & m0.minus, m1.nonzeros};
}

template<typename T>
constexpr Multiply_P0<T> Compose<Multiply_P0<T>, Multiply_P0<T>>::operator()(
    const Multiply_P0<T>& m1, const Multiply_P0<T>& m0
) {
    return {m1.nonzeros & m0.nonzeros};
}

template<typename T>
constexpr Multiply_PM0<T> Compose<Multiply<T>, Multiply_PM0<T>>::operator()(
    const Multiply<T>& m1, const Multiply_PM0<T>& m0
) {
    return {m1.minus ^ m0.minus, m0.nonzeros};
}

template<typename T>
constexpr Multiply_PM0<T> Compose<Multiply_PM0<T>, Multiply<T>>::operator()(
    const Multiply_PM0<T>& m1, const Multiply<T>& m0
) {
    return {m1.minus ^ m0.minus, m1.nonzeros};
}

template<typename T>
constexpr Multiply_PM0<T> Compose<Multiply_P0<T>, Multiply_PM0<T>>::operator()(
    const Multiply_P0<T>& m1, const Multiply_PM0<T>& m0
) {
    return {m1.nonzeros & m0.minus, m1.nonzeros & m0.nonzeros};
}

template<typename T>
constexpr Multiply_PM0<T> Compose<Multiply_PM0<T>, Multiply_P0<T>>::operator()(
    const Multiply_PM0<T>& m1, const Multiply_P0<T>& m0
) {
    return {m1.minus & m0.nonzeros, m1.nonzeros & m0.nonzeros};
}

template<typename T>
constexpr Multiply_PM0<T> Compose<Multiply_PM0<T>, Multiply_PM0<T>>::operator()(
    const Multiply_PM0<T>& m1, const Multiply_PM0<T>& m0
) {
    return {(m1.minus & m0.nonzeros) ^ (m0.minus & m1.nonzeros), 
        m1.nonzeros & m0.nonzeros};
}

template<typename T0, typename T1>
constexpr composes_to<Multiply<T1>, PushForward<T0, T1>>
Compose<Multiply<T1>, PushForward<T0, T1>>::operator()(
    const Multiply<T1>& m, const PushForward<T0, T1>& pf
) {
    return {pf, {pull_back(pf.function, m.minus)}};
}

template<typename T0, typename T1>
constexpr composes_to<Multiply_P0<T1>, PushForward<T0, T1>>
Compose<Multiply_P0<T1>, PushForward<T0, T1>>::operator()(
    const Multiply_P0<T1>& m, const PushForward<T0, T1>& pf
) {
    return {pf, {pull_back(pf.function, m.nonzeros)}};
}

template<typename T0, typename T1>
constexpr composes_to<Multiply_PM0<T1>, PushForward<T0, T1>>
Compose<Multiply_PM0<T1>, PushForward<T0, T1>>::operator()(
    const Multiply_PM0<T1>& m, const PushForward<T0, T1>& pf
) {
    return {pf, {
        pull_back(pf.function, m.minus), 
        pull_back(pf.function, m.nonzeros)
    }};
}

template<typename T0, typename T1>
constexpr composes_to<PullBack<T0, T1>, Multiply<T0>>
Compose<PullBack<T0, T1>, Multiply<T0>>::operator()(
    const PullBack<T0, T1>& pb, const Multiply<T0>& m
) {
    return {{pull_back(pb.function, m.minus)}, pb};
}

template<typename T0, typename T1>
constexpr composes_to<PullBack<T0, T1>, Multiply_P0<T0>>
Compose<PullBack<T0, T1>, Multiply_P0<T0>>::operator()(
    const PullBack<T0, T1>& pb, const Multiply_P0<T0>& m
) {
    return {{pull_back(pb.function, m.nonzeros)}, pb};
}

template<typename T0, typename T1>
constexpr composes_to<PullBack<T0, T1>, Multiply_PM0<T0>>
Compose<PullBack<T0, T1>, Multiply_PM0<T0>>::operator()(
    const PullBack<T0, T1>& pb, const Multiply_PM0<T0>& m
) {
    return {{
        pull_back(pb.function, m.minus), 
        pull_back(pb.function, m.nonzeros)
    }, pb};
}

template<typename T0, typename T1, typename T2>
constexpr PushForward<T0, T2> 
Compose<PushForward<T1, T2>, PushForward<T0, T1>>::operator()(
    const PushForward<T1, T2>& pf12,
    const PushForward<T0, T1>& pf01
) {
    return {pf12.function * pf01.function};
}

template<typename T0, typename T1, typename T2>
constexpr PullBack<T0, T2> 
Compose<PullBack<T1, T2>, PullBack<T0, T1>>::operator()(
    const PullBack<T1, T2>& pb12,
    const PullBack<T0, T1>& pb01
) {
    return {pb01.function * pb12.function};
}

template<typename T0, typename T1, typename T2>
constexpr composes_to<PullBack<T1, T2>, PushForward<T0, T1>>
Compose<PullBack<T1, T2>, PushForward<T0, T1>>::operator()(
    const PullBack<T1, T2>& pb, const PushForward<T0, T1>& pf
) {
    int idx_pushed_into[T1::LENGTH];
    for (auto idx = 0; idx < T1::LENGTH; idx++) {
        idx_pushed_into[idx] = -1;
    }
    for (auto idx = 0; idx < T0::LENGTH; idx++) {
        idx_pushed_into[pf.function(idx)] = idx;
    }
    bit_vector<T2::LENGTH> nonzeros;
    std::array<int, T2::LENGTH> f = {};
    for (auto idx = 0; idx < T2::LENGTH; idx++) {
        int x = idx_pushed_into[pb.function(idx)];
        nonzeros.set_bit(idx, x != -1);
        f[idx] = (x == -1)? 0 : x;
    }
    return {{nonzeros}, {{f}}};
}


}