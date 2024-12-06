#pragma once
#include <cstdint>
#include <concepts>
#include <array>
#include "mymath.hpp"
#include "signvectors.hpp"

namespace sign_vector_operations {

// ===============================================
//                  index_function
//   and elementary push-forwards and pull-backs
// ===============================================

// Represents a function `f: {0..L0-1} -> {0..L1-1}` by
// storing the composite functions `f/32` and `f%32`. 
template<int L0, int L1>
struct index_function {
    // =============
    //   CONSTANTS
    // =============

    // The number of elements in the domain.
    constexpr static const int LENGTH = L0;
    // Number of (unsigned) 32-bit integers needed to store as many
    // bits as there are elements in the domain.
    constexpr static const int NR_INT32 = division_rounded_up(LENGTH, 32);
    // Grouping the elements in the domain in 32-long chunks, this 
    // many signs are left for the last group (0 < n <= 32).
    constexpr static const int NR_REMAINING_BITS = LENGTH - 32 * (NR_INT32 - 1);

    // =============
    //   VARIABLES
    // =============

    // `i[idx] == f(idx)/32` if `f` is the function being represented.
    int i[L0];
    // `r[idx] == f(idx)%32` if `f` is the function being represented.
    char r[L0];

    // ================
    //   CONSTRUCTORS
    // ================

    // Initializes the constant 0 function.
    constexpr index_function(): i{}, r{} {}
    // Given an `Iterable` of length `L0` and of entries in
    // `0..L1`, consider `f` to be represented by it, and
    // initialize this `f`.
    template<typename Iterable>
    constexpr index_function(const Iterable& iterable) {
        int idx = 0;
        for (auto value: iterable) {
            i[idx] = value >> 5;
            r[idx] = value & 31;
            idx++;
        }    
    }

    // ===============================
    //   WRAPPED ACCESS TO VARIABLES
    // ===============================

    // Evaluates the function on a given input.
    constexpr int operator()(int idx) const;

    // ================================
    //   CONSTRUCT NEW INDEX FUNTIONS
    // ================================
    
    template<int L>
    constexpr index_function<L,L1> compose(const index_function<L,L0>&) const;
    template<int L>
    constexpr index_function<L,L1> operator*(const index_function<L,L0>&) const; 
};

// Given an `index_function` representing an injective 
// function `f: {0..L-1} -> {0..L1-1}`, push 
// forward this bitvector along it: the resulting
// bitvector at index `idx` should be equal to the input
// bitvector at index obtained as the preimage of `idx` 
// under `f` if this makes sense, and otherwise the
// bitvector at this index should be `0`.
// 
// Behaviour undefined for non-injective functions.
template<int L0, int L1>
bit_vector<L1> push_forward(
    const index_function<L0,L1>&, const bit_vector<L0>&
);
// Given an `index_function` representing a
// function `f: {0..L0-1} -> {0..L-1}` pull back this
// bitvector along it: viewing the input bitvector as
// a function `b: {0..L-1} -> {0, 1}`, the result is
// `b` composed with `f`. 
template<int L0, int L1>
bit_vector<L0> pull_back(
    const index_function<L0,L1>&, const bit_vector<L1>&
);
// Given an `index_function` representing an injective 
// function `f: {0..L-1} -> {0..L1-1}`, push 
// forward this signvector along it: the resulting
// signvector at index `idx` should be equal to the input
// signvector at index obtained as the preimage of `idx` 
// under `f` if this makes sense, and otherwise the
// signvector at this index should be `0`.
// 
// Behaviour undefined for non-injective functions.
template<int L0, int L1>
sign_vector<L1> push_forward(
    const index_function<L0,L1>&, const sign_vector<L0>&
);
// Given an `index_function` representing a
// function `f: {0..L0-1} -> {0..L-1}` pull back this
// signvector along it: viewing the input signvector as
// a function `b: {0..L-1} -> {-,0,+}`, the result is
// `b` composed with `f`. 
template<int L0, int L1>
sign_vector<L0> pull_back(
    const index_function<L0,L1>&, const sign_vector<L1>&
);

// =======================
//        Operations
//   and formal products
// =======================

// Represents an operation that transforms a term of type
// `T0` into a term of type `T1`. Both `T0` and `T1` are assumed
// to be wrappers of `bit_vector` or `sign_vector`.
template<typename T0, typename T1>
struct Operation {
    constexpr static const int DOMAIN_LENGTH = T0::LENGTH;
    constexpr static const int CODOMAIN_LENGTH = T1::LENGTH;
    using DOMAIN_TYPE = T0;
    using CODOMAIN_TYPE = T1;
};

// This lets us include a type description in a template argument
// without actually doing anything for it. This is useful for partial
// specialization.
template<typename First, typename Second>
using first_type = First;

// `IsOperation<O>::value` holds if and only if `O` is a derived class
// of `Operation` with some template parameters. In case `O` does not
// even have `DOMAIN_TYPE` and `CODOMAIN_TYPE` typedefs, it can not be
// derived from `Operation`.
template<typename O, typename DOMAIN_TYPE=void, typename CODOMAIN_TYPE=void>
struct IsOperation {
    constexpr static const bool value = false;
};

// `IsOperation<O>::value` holds if and only if `O` is a derived class
// of `Operation` with some template parameters. In case `O` has
// `DOMAIN_TYPE` and `CODOMAIN_TYPE` typedefs, whether it is derived
// from `O` can be determined using `std::derived_from`.
template<typename O>
struct IsOperation<O, 
first_type<void, typename O::DOMAIN_TYPE>, 
first_type<void, typename O::CODOMAIN_TYPE>> {
    constexpr static const bool value = std::derived_from<O, Operation<
        typename O::DOMAIN_TYPE,
        typename O::CODOMAIN_TYPE
    >>;
};

// `is_operation<O>` holds if and only if `O` is a derived class
//  of `Operation` with some template parameters.
template<typename O>
constexpr bool is_operation = IsOperation<O>::value;

// `AreComposableOperations<Op12,Op01>::value` holds if `Op12` and
// `Op01` are both derived classes of `Operation` and the
// codomain type of `Op01` coincides with the domain type of
// `Op12`.
template<typename Op12, typename Op01, 
    bool AreOperations = is_operation<Op12> && is_operation<Op01>>
struct AreComposableOperations {
    constexpr static const bool value = false;
};

template<typename Op12, typename Op01>
struct AreComposableOperations<Op12, Op01, true> {
    constexpr static const bool value = std::is_same_v<
        typename Op01::CODOMAIN_TYPE,
        typename Op12::DOMAIN_TYPE
    >;
};

// `are_composable_operations<Op12,Op01>` holds if `Op12` and
// `Op01` are both derived classes of `Operation` and the
// codomain type of `Op01` coincides with the domain type of
// `Op12`.
template<typename Op12, typename Op01>
constexpr bool are_composable_operations = AreComposableOperations<Op12,Op01>::value;

// Represents a formal product (composition) of two `Operation`s
// of types `Op12` and `Op01`, where `Op01::CODOMAIN_TYPE` and
// `Op12::DOMAIN_TYPE` coincide.
template<typename Op12, typename Op01>
struct FormalProductOfOperations: 
public Operation<
    typename Op01::DOMAIN_TYPE,
    typename Op12::CODOMAIN_TYPE
> {
    constexpr static const int L0 = Op01::DOMAIN_LENGTH;
    constexpr static const int L1 = Op01::CODOMAIN_LENGTH;
    constexpr static const int L2 = Op12::CODOMAIN_LENGTH;
    using T0 = typename Op01::DOMAIN_TYPE;
    using T1 = typename Op01::CODOMAIN_TYPE;
    using T2 = typename Op12::CODOMAIN_TYPE;

    static_assert(are_composable_operations<Op12,Op01>);

    Op12 op12;
    Op01 op01;

    constexpr FormalProductOfOperations(): op12{}, op01{} {}
    constexpr FormalProductOfOperations(
        const Op12& o12, const Op01& o01
    ): op12(o12), op01(o01) {}

    constexpr T2 applied(const T0& to) const { 
        return op12.applied(op01.applied(to));
    }
    constexpr T2 operator*(const T0& arg) const { return applied(arg); }
    constexpr T2 operator()(const T0& arg) const { return applied(arg); }
};

// Has a single typedef `type` as a member that specifies what is
// the type of the composition of two `Operation`s of types `Op12`
// and `Op01`. Defaults to a formal product type, but specializations
// can specify "smarter" compositions.
template<typename Op12, typename Op01, typename Condition=
    std::enable_if_t<are_composable_operations<Op12,Op01>>>
struct ComposesTo {
    using type = FormalProductOfOperations<Op12,Op01>;
};

// Specifies what is the type of the composition of two `Operation`s 
// of types `Op12` and `Op01`.
template<typename Op12, typename Op01>
using composes_to = typename ComposesTo<Op12,Op01>::type;

// True if and only if there is a "smart" way to compose `Operation`s
// of types `Op12` and `Op01`, i.e. if their composition is not just
// their formal product.
template<typename Op12, typename Op01>
constexpr bool composes_nontrivially = !std::is_same_v<
    composes_to<Op12, Op01>,
    FormalProductOfOperations<Op12,Op01>
>;

// A static functor which performs the composition of two `Operation`s
// of types `Op12` and `Op01`. Can be specialized to provide
// "smarter" ways to combine specific operator types.
template<typename Op12, typename Op01, typename Condition=void>
struct Compose {
    constexpr static FormalProductOfOperations<Op12,Op01> operator()(
        const Op12& op12, const Op01& op01
    ) {
        return {op12, op01};
    }
};

// A function (an instance of a static functor really) which performs
// the composition of two `Operation`s of types `Op12` and `Op01`.
template<typename Op12, typename Op01>
constexpr Compose<Op12,Op01> compose;

// Compose two operators, possibly only formally.
template<typename Op12, typename Op01, std::enable_if_t<
    are_composable_operations<Op12,Op01>, bool> = true>
composes_to<Op12, Op01> operator*(
    const Op12& op12,
    const Op01& op01
) {
    return compose<Op12, Op01>(op12, op01);
}

// Composition should be parenthesised to the left; hence
// formal products on the right hand side of a composition
// must decompose.
template<typename Op23, typename Op12, typename Op01>
struct ComposesTo<Op23, FormalProductOfOperations<Op12, Op01>> {
    using type = composes_to<composes_to<Op23, Op12>, Op01>;
};

// Composition should be parenthesised to the left; hence
// formal products on the right hand side of a composition
// must decompose.
template<typename Op23, typename Op12, typename Op01>
struct Compose<Op23, FormalProductOfOperations<Op12, Op01>> {
    constexpr static composes_to<composes_to<Op23, Op12>, Op01> operator()(
        const Op23& op23,
        const FormalProductOfOperations<Op12, Op01>& op02
    ) {
        return compose<composes_to<Op23, Op12>, Op01>(
            compose<Op23, Op12>(op23, op02.op12),
            op02.op01
        );
    }
};

// If breaking left-parenthesisation would result in a "smarter"
// composition, do so: then first compose `Op12` and `Op01`, and
// then postcompose `Op23`.
template<typename Op23, typename Op12, typename Op01>
struct ComposesTo<FormalProductOfOperations<Op23, Op12>, Op01,
    std::enable_if_t<composes_nontrivially<Op12, Op01>>> {
    using type = composes_to<Op23, composes_to<Op12, Op01>>;
};

// If breaking left-parenthesisation would result in a "smarter"
// composition, do so: then first compose `Op12` and `Op01`, and
// then postcompose `Op23`.
template<typename Op23, typename Op12, typename Op01>
struct Compose<FormalProductOfOperations<Op23, Op12>, Op01,
    std::enable_if_t<composes_nontrivially<Op12, Op01>>> {
    constexpr static composes_to<Op23, composes_to<Op12, Op01>> operator()(
        const FormalProductOfOperations<Op23, Op12>& op13,
        const Op01& op01
    ) {
        return compose<Op23, composes_to<Op12, Op01>>(
            op13.op12,
            compose<Op12, Op01>(op13.op01, op01)
        );
    }
};

// ============================
//   Specific operation types
// ============================

template<typename T>
struct Multiply: public Operation<T,T> {
    constexpr static const int L = T::LENGTH;

    bit_vector<L> minus;

    constexpr Multiply(): minus() {}
    constexpr Multiply(const bit_vector<L>& m): minus(m) {}

    constexpr T applied(const T&) const;
    constexpr T operator*(const T& arg) const { return applied(arg); }
    constexpr T operator()(const T& arg) const { return applied(arg); }
};

template<typename T>
struct Multiply_P0: public Operation<T,T> {
    constexpr static const int L = T::LENGTH;

    bit_vector<L> nonzeros;

    constexpr Multiply_P0(): nonzeros{} {}
    constexpr Multiply_P0(
        const bit_vector<L>& nz
    ): nonzeros(nz) {}
    constexpr Multiply_P0( // Do we need this?
        const Multiply<T>& m
    ): nonzeros(m.minus) {}

    constexpr T applied(const T&) const;
    constexpr T operator*(const T& arg) const { return applied(arg); }
    constexpr T operator()(const T& arg) const { return applied(arg); }
};

template<typename T>
struct Multiply_PM0: public Operation<T,T> {
    constexpr static const int L = T::LENGTH;

    bit_vector<L> minus;
    bit_vector<L> nonzeros;

    constexpr Multiply_PM0(): minus{}, nonzeros{} {}
    constexpr Multiply_PM0(
        const bit_vector<L>& m, 
        const bit_vector<L>& nz
    ): minus(m), nonzeros(nz) {}
    constexpr Multiply_PM0(
        const Multiply<T>& m
    ): minus(m.minus), nonzeros(m.minus) {}

    constexpr T applied(const T&) const;
    constexpr T operator*(const T& arg) const { return applied(arg); }
    constexpr T operator()(const T& arg) const { return applied(arg); }
};

template<typename T0, typename T1>
struct PushForward: Operation<T0,T1> {
    constexpr static const int L0 = T0::LENGTH;
    constexpr static const int L1 = T1::LENGTH;

    index_function<L0, L1> function;

    constexpr PushForward(): function() {}
    constexpr PushForward(const index_function<L0, L1>& f): function(f) {}

    constexpr T1 applied(const T0&) const;
    constexpr T1 operator*(const T0& arg) const { return applied(arg); }
    constexpr T1 operator()(const T0& arg) const { return applied(arg); }
};

template<typename T0, typename T1>
struct PullBack: Operation<T0,T1> {
    constexpr static const int L0 = T0::LENGTH;
    constexpr static const int L1 = T1::LENGTH;

    index_function<L1, L0> function;

    constexpr PullBack(): function() {}
    constexpr PullBack(const index_function<L1, L0>& f): function(f) {}

    constexpr T1 applied(const T0&) const;
    constexpr T1 operator*(const T0& arg) const { return applied(arg); }
    constexpr T1 operator()(const T0& arg) const { return applied(arg); }
};

// ==================================
//   Composition of operation types
// ==================================

template<typename T>
struct ComposesTo<Multiply<T>, Multiply<T>> {
    using type = Multiply<T>;
};
template<typename T>
struct Compose<Multiply<T>, Multiply<T>> {
    constexpr static Multiply<T> operator()(
        const Multiply<T>&, const Multiply<T>&
    );
};

template<typename T>
struct ComposesTo<Multiply<T>, Multiply_P0<T>> {
    using type = Multiply_PM0<T>;
};
template<typename T>
struct Compose<Multiply<T>, Multiply_P0<T>> {
    constexpr static Multiply_PM0<T> operator()(
        const Multiply<T>&, const Multiply_P0<T>&
    );
};

template<typename T>
struct ComposesTo<Multiply_P0<T>, Multiply<T>> {
    using type = Multiply_PM0<T>;
};
template<typename T>
struct Compose<Multiply_P0<T>, Multiply<T>> {
    constexpr static Multiply_PM0<T> operator()(
        const Multiply_P0<T>&, const Multiply<T>&
    );
};

template<typename T>
struct ComposesTo<Multiply_P0<T>, Multiply_P0<T>> {
    using type = Multiply_P0<T>;
};
template<typename T>
struct Compose<Multiply_P0<T>, Multiply_P0<T>> {
    constexpr static Multiply_P0<T> operator()(
        const Multiply_P0<T>&, const Multiply_P0<T>&
    );
};

template<typename T>
struct ComposesTo<Multiply<T>, Multiply_PM0<T>> {
    using type = Multiply_PM0<T>;
};
template<typename T>
struct Compose<Multiply<T>, Multiply_PM0<T>> {
    constexpr static Multiply_PM0<T> operator()(
        const Multiply<T>&, const Multiply_PM0<T>&
    );
};

template<typename T>
struct ComposesTo<Multiply_PM0<T>, Multiply<T>> {
    using type = Multiply_PM0<T>;
};
template<typename T>
struct Compose<Multiply_PM0<T>, Multiply<T>> {
    constexpr static Multiply_PM0<T> operator()(
        const Multiply_PM0<T>&, const Multiply<T>&
    );
};

template<typename T>
struct ComposesTo<Multiply_P0<T>, Multiply_PM0<T>> {
    using type = Multiply_PM0<T>;
};
template<typename T>
struct Compose<Multiply_P0<T>, Multiply_PM0<T>> {
    constexpr static Multiply_PM0<T> operator()(
        const Multiply_P0<T>&, const Multiply_PM0<T>&
    );
};

template<typename T>
struct ComposesTo<Multiply_PM0<T>, Multiply_P0<T>> {
    using type = Multiply_PM0<T>;
};
template<typename T>
struct Compose<Multiply_PM0<T>, Multiply_P0<T>> {
    constexpr static Multiply_PM0<T> operator()(
        const Multiply_PM0<T>&, const Multiply_P0<T>&
    );
};

template<typename T>
struct ComposesTo<Multiply_PM0<T>, Multiply_PM0<T>> {
    using type = Multiply_PM0<T>;
};
template<typename T>
struct Compose<Multiply_PM0<T>, Multiply_PM0<T>> {
    constexpr static Multiply_PM0<T> operator()(
        const Multiply_PM0<T>&, const Multiply_PM0<T>&
    );
};

template<typename T0, typename T1>
struct ComposesTo<Multiply<T1>, PushForward<T0, T1>> {
    using type = FormalProductOfOperations<PushForward<T0, T1>, Multiply<T0>>;
};
template<typename T0, typename T1>
struct Compose<Multiply<T1>, PushForward<T0, T1>> {
    constexpr static composes_to<Multiply<T1>, PushForward<T0, T1>> operator()(
        const Multiply<T1>&, const PushForward<T0, T1>&
    );
};

template<typename T0, typename T1>
struct ComposesTo<Multiply_P0<T1>, PushForward<T0, T1>> {
    using type = FormalProductOfOperations<PushForward<T0, T1>, Multiply_P0<T0>>;
};
template<typename T0, typename T1>
struct Compose<Multiply_P0<T1>, PushForward<T0, T1>> {
    constexpr static composes_to<Multiply_P0<T1>, PushForward<T0, T1>> operator()(
        const Multiply_P0<T1>&, const PushForward<T0, T1>&
    );
};

template<typename T0, typename T1>
struct ComposesTo<Multiply_PM0<T1>, PushForward<T0, T1>> {
    using type = FormalProductOfOperations<PushForward<T0, T1>, Multiply_PM0<T0>>;
};
template<typename T0, typename T1>
struct Compose<Multiply_PM0<T1>, PushForward<T0, T1>> {
    constexpr static composes_to<Multiply_PM0<T1>, PushForward<T0, T1>> operator()(
        const Multiply_PM0<T1>&, const PushForward<T0, T1>&
    );
};

template<typename T0, typename T1>
struct ComposesTo<PullBack<T0,T1>, Multiply<T0>> {
    using type = FormalProductOfOperations<Multiply<T1>, PullBack<T0, T1>>;
};
template<typename T0, typename T1>
struct Compose<PullBack<T0, T1>, Multiply<T0>> {
    constexpr static composes_to<PullBack<T0, T1>, Multiply<T0>> operator()(
        const PullBack<T0, T1>&, const Multiply<T0>&
    );
};

template<typename T0, typename T1>
struct ComposesTo<PullBack<T0,T1>, Multiply_P0<T0>> {
    using type = FormalProductOfOperations<Multiply_P0<T1>, PullBack<T0, T1>>;
};
template<typename T0, typename T1>
struct Compose<PullBack<T0, T1>, Multiply_P0<T0>> {
    constexpr static composes_to<PullBack<T0, T1>, Multiply_P0<T0>> operator()(
        const PullBack<T0, T1>&, const Multiply_P0<T0>&
    );
};

template<typename T0, typename T1>
struct ComposesTo<PullBack<T0,T1>, Multiply_PM0<T0>> {
    using type = FormalProductOfOperations<Multiply_PM0<T1>, PullBack<T0, T1>>;
};
template<typename T0, typename T1>
struct Compose<PullBack<T0, T1>, Multiply_PM0<T0>> {
    constexpr static composes_to<PullBack<T0, T1>, Multiply_PM0<T0>> operator()(
        const PullBack<T0, T1>&, const Multiply_PM0<T0>&
    );
};

template<typename T0, typename T1, typename T2>
struct ComposesTo<PushForward<T1, T2>, PushForward<T0, T1>> {
    using type = PushForward<T0, T2>;
};
template<typename T0, typename T1, typename T2>
struct Compose<PushForward<T1, T2>, PushForward<T0, T1>> {
    constexpr static PushForward<T0, T2> operator()(
        const PushForward<T1, T2>&, const PushForward<T0, T1>&
    );
};

template<typename T0, typename T1, typename T2>
struct ComposesTo<PullBack<T1, T2>, PullBack<T0, T1>> {
    using type = PullBack<T0, T2>;
};
template<typename T0, typename T1, typename T2>
struct Compose<PullBack<T1, T2>, PullBack<T0, T1>> {
    constexpr static PullBack<T0, T2> operator()(
        const PullBack<T1, T2>&, const PullBack<T0, T1>&
    );
};

template<typename T0, typename T1, typename T2>
struct ComposesTo<PullBack<T1, T2>, PushForward<T0, T1>> {
    using type = FormalProductOfOperations<Multiply_P0<T2>, PullBack<T0, T2>>;
};
template<typename T0, typename T1, typename T2>
struct Compose<PullBack<T1, T2>, PushForward<T0, T1>> {
    constexpr static composes_to<PullBack<T1, T2>, PushForward<T0, T1>> operator()(
        const PullBack<T1, T2>&, const PushForward<T0, T1>&
    );
};

}

#include "signvectoroperations.cpp"
