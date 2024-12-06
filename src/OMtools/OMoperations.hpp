#pragma once
#include <vector>
#include <array>
#include "OMs.hpp"
#include "mymath.hpp"
#include "NchooseK.hpp"
#include "signvectors.hpp"
#include "signvectoroperations.hpp"

// Conventions:
// - when not specified, the template parameters are the template
//      parameters of the input chirotopes.

// Standard functors ("chirotope-operations" and "matroid-operations") 
// constructing new matroids, chirotopes from old ones.
//
// Unless specified otherwise, the rank (`R` or `R0`) and number
// of elements (`N` or `N0`) template parameters parametrize the
// input of a functor defined in this namespace.
namespace OM_operations {

namespace svo=sign_vector_operations;

// =========
// UTILITIES
// =========

// Constructs an array `f` where `f[i]` increases one by
// one as `i` increases, starting from 0, but skips the
// specified integer.
//
// Example: 
// `monotone_function_skipping_number<int, 5>(2)`
// is `std::array<int, 5>{0,1,3,4,5}`.
template<typename NumberType, NumberType N>
std::array<NumberType, N> monotone_function_skipping_number(NumberType);

// Constructs an array `f` where `f[i]` increases one by
// one as `i` increases, starting from 0, but skips all
// integers which are listed in `numbers`.
//
// Example: 
// `monotone_function_skipping_numbers<char, 5>(std::vector<char>{1,4,5})`
// is `std::array<char, 5>{0,2,3,6,7}`.
//
// The type `Iterable` must be range-based-for-loop iterable
// where it yields `NumberType`s.
template<typename NumberType, NumberType N, typename Iterable>
std::array<NumberType, N> monotone_function_skipping_numbers(const Iterable&);


// Given an injective function `f` from `0..N0` to `0..N1`, 
// there is an induced function from the set of `R`-element 
// subsets of `0..N0` to the set of `R`-element subsets of
// `0..N1`. Both these sets of subsets have canonical 
// orderings defined in `NchooseK`, and so can be identified with
// `0..binomial_coefficient(N0,R)` and `0..binomial_coefficient(N1,R)`
// respectively. This method computes the function induced
// by `f` between these two sets of integer.
//
// Both `f` and the output are encoded as a iterables, so 
// that the `i`th item in the iterable is `f(i)`, and
// similarly for the output.
template<typename SmallLabelType, typename BigLabelType, 
SmallLabelType R, SmallLabelType N0, BigLabelType N1, typename Iterable, 
typename IndexType=BigLabelType, typename FlatIndexType=IndexType>
constexpr std::array<IndexType, binomial_coefficient(N0, R)>
Rtuple_function_induced_by_relabeling(const Iterable&);

// This is identical to `Rtuple_function_induced_by_relabeling`
// but assumes that the input is strictly monotone. In particular:
//
// Given a strictly monotone function `f` from `0..N0` to `0..N1`, 
// there is an induced function from the set of `R`-element 
// subsets of `0..N0` to the set of `R`-element subsets of
// `0..N1`. Both these sets of subsets have canonical 
// orderings defined in `NchooseK`, and so can be identified with
// `0..binomial_coefficient(N0,R)` and `0..binomial_coefficient(N1,R)`
// respectively. This method computes the function induced
// by `f` between these two sets of integer.
//
// Both `f` and the output are encoded as a iterables, so 
// that the `i`th item in the iterable is `f(i)`, and
// similarly for the output.
template<typename SmallLabelType, typename BigLabelType, 
SmallLabelType R, SmallLabelType N0, BigLabelType N1, typename Iterable, 
typename IndexType=BigLabelType, typename FlatIndexType=IndexType>
constexpr std::array<IndexType, binomial_coefficient(N0, R)>
Rtuple_function_induced_by_monotone_relabeling(const Iterable&);


// =============
// REORIENTATION
// =============

// Constructs a chirotope-operator, which reorients the given 
// element of any rank `R` chirotope on set of elements `0..N-1`. 
// In other words, this chirotope-operator changes the sign of 
// all `R`-tuples which contain the given element.
template<int R, int N>
constexpr svo::Multiply<Chirotope<R, N>> reorient_element(int);

// Constructs a chirotope-operator, which reorients the given 
// elements of any rank `R` chirotope on set of elements `0..N-1`. 
// In other words, this chirotope-operator changes the sign of 
// all `R`-tuples which contain an even number of the given set 
// of elements.
//
// The set of elements must be provided in a data structure which
// can be iterated using a range-based for loop, such that it
// yields `int`s.
template<int R, int N, typename Iterable>
constexpr svo::Multiply<Chirotope<R, N>> 
reorient_elements(const Iterable&);

// ==========
// RELABELING
// ==========

// Constructs a chirotope-operator given an injective function 
// `f` from `0..N0-1` to `0..N1-1`, which computes the relabeling
// `chi1` of any chirotope `chi0` of rank `R` on set of elements 
// `0..N0-1`. In other words, this chirotope-operator computes
// a chirotope `chi1` of rank `R` on set of elements `0..N1-1`
// such that `chi0` evaluated on the `R`-tuple `(x1,...,xR)` is
// always equal to `chi1` evaluated on `(f(x1),...,f(xR))`, and
// the latter is `0` on all `R`-tuples which are not of this form.
//
// The function `f` is represented by a data structure which can
// be iterated using a range-based for loop yielding `int`s. The
// structure yields `f(i)` on the `i`th iteration, starting from
// `i=0`.
template<int R, int N0, int N1=N0, typename Iterable=std::vector<int>>
constexpr svo::composes_to<svo::PushForward<Chirotope<R, N0>, Chirotope<R, N1>>,svo::Multiply<Chirotope<R,N0>>>
relabel_elements(const Iterable&);

// This behaves identically to `relabel_elements`, except that
// it has undefined behaviour whenever the specified function
// is not monotonely increasing. In exchange it produces a simpler
// and faster chirotope-operation.
template<int R, int N0, int N1, typename Iterable>
constexpr svo::PushForward<Chirotope<R, N0>, Chirotope<R, N1>>
relabel_monotonely(const Iterable&);

// Constructs a chirotope-operator given an injective function 
// `f` from `0..N1-1` to `0..N0-1`, which computes to a
// chirotope `chi0` of rank `R0` on set of elements `0..N0-1`
// the chirotope `chi1` of rank `R1` on set of elements `0..N1-1`
// given by first contracting all elements not in the image of
// `f`, and then relabeling by the inverse of `f` as described
// in `relabel_elements`.
//
// The function `f` is represented by a data structure which can
// be iterated using a range-based for loop yielding `int`s. The
// structure yields `f(i)` on the `i`th iteration, starting from
// `i=0`.
template<int R0, int N0, int R1, int N1, typename Iterable>
constexpr svo::composes_to<svo::Multiply<Chirotope<R1,N1>>,svo::PullBack<Chirotope<R0,N0>,Chirotope<R1,N1>>>
general_minor(const Iterable&);

// This behaves identically to `general_minor`, except that
// it has undefined behaviour whenever the specified function
// is not monotonely increasing. In exchange it produces a simpler
// and faster chirotope-operation.
template<int R0, int N0, int R1, int N1, typename Iterable>
constexpr svo::PullBack<Chirotope<R0, N0>, Chirotope<R1, N1>>
monotone_minor(const Iterable&);

// Constructs a chirotope-operator, which inserts a loop with
// the given label into any rank `R` chirotope on set of
// elements `0..N-1`. All elements whose labels were originally
// larger than or equal to the given integer will be shifted up
// by one.
//
// The left inverse of this operation is `delete_label`:
// `delete_label(i)*insert_loop_at(i)` is the identity.
template<int R, int N>
constexpr svo::PushForward<Chirotope<R,N>,Chirotope<R,N+1>>
insert_loop_at(int);

// Constructs a chirotope-operator, which inserts a coloop with
// the given label into any rank `R` Chirotope on set of
// elements `0..N-1`, thereby increasing the rank to `R+1` and
// the set of elements to `0..N`. All elements whose labels were
// originally larger than or equal to the given integer will be
// shifted up by one.
//
// The left inverse of this operation is `contract_label`:
// `contract_label(i)*insert_coloop_at(i)` is the identity.
template<int R, int N>
constexpr svo::PushForward<Chirotope<R,N>,Chirotope<R+1,N+1>>
insert_coloop_at(int);

// =====================
// DELETION, CONTRACTION
// =====================

// Constructs a chirotope-operator, which deletes the given
// element of any chirotope of rank `R` on set of elements `0..N-1`.
// In other words, it sets the chirotope to `0` on all
// `R`-tuples which contain the given element, without changing
// the rank or the set of elements of said chirotope.
//
// See also the chirotope-operator `delete_label`, which
// produces chirotopes on the set of elements `0..N-2`, by
// forgetting the given element as well and relabeling while
// preserving the order of elements. See also `delete_elements`
// for deleting multiple elements simultaneously.
template<int R, int N>
constexpr svo::Multiply_P0<Chirotope<R, N>> delete_element(int);

// Constructs a chirotope-operator, which deletes the given
// label of any chirotope of rank `R` on set of elements `0..N-1`.
// In other words, it forgets the given element, and relabels
// the others in an order preserving way, to get a chirotope
// of rank `R` on set of elements `0..N-2`.
//
// See also the chirotope-operator `delete_element`, which
// produces chirotopes on the set of elements `0..N-1`, by
// setting the chirotope to `0` on all `R`-tuples which
// contain the given element. See also `delete_labels` for 
// deleting multiple labels simultaneously.
template<int R, int N>
constexpr svo::PullBack<Chirotope<R,N>,Chirotope<R,N-1>> 
delete_label(int);

// Constructs a chirotope-operator, which deletes the given
// set of elements of any chirotope of rank `R` on set of 
// elements `0..N-1`. In other words, it sets the chirotope 
// to `0` on all `R`-tuples which contain at least one given 
// element, without changing the rank or the set of elements 
// of said chirotope.
//
// See also the chirotope-operator `delete_labels`, which does
// the same, but afterwards forgets all elements in the given 
// set, and relabels while preserving the order of elements,
// thus changing the set of elements of the resulting chirotope.
// See also `delete_element` for deleting only a single element
// at a time.
template<int R, int N, typename Iterable>
constexpr svo::Multiply_P0<Chirotope<R, N>> 
delete_elements(const Iterable&);

// Constructs a chirotope-operator, which deletes the given
// set of labels of any chirotope of rank `R` on set of 
// elements `0..N0-1`. In other words, it forgets about the
// labels in the given set, and relabels the remaining elements
// in an order preserving way to get a new chirotope of rank
// `R` on set of elements `0..N1-1`.
//
// See also the chirotope-operator `delete_elements`, produces
// a chirotope on the set of elements `0..N0-1` by setting it
// to `0` on all `R`-tuples which contain at least one given
// element. See also `delete_label` for deleting only a single 
// label at a time.
template<int R, int N0, int N1, typename Iterable>
constexpr svo::PullBack<Chirotope<R,N0>,Chirotope<R,N1>> 
delete_labels(const Iterable&);

// Constructs a chirotope-operator, which contracts the given
// element of any chirotope of rank `R` on set of elements `0..N-1`.
// In other words, it sets the chirotope to `0` on all
// `R`-tuples which do not contain the given element, without 
// changing the rank or the set of elements of said chirotope.
//
// See also the chirotope-operator `contract_label`, which
// produces chirotopes on the set of elements `0..N-2`, by
// forgetting the given element as well and relabeling while
// preserving the order of elements. See also `contract_elements`
// for contracting multiple elements simultaneously.
template<int R, int N>
constexpr svo::Multiply_P0<Chirotope<R, N>> contract_element(int);

// Constructs a chirotope-operator, which contracts the given
// label of any chirotope of rank `R` on set of elements `0..N-1`,
// producing a chirotope of rank `R-1` on set of elements 
// `0..N-2`. In other words, it forgets about the given
// label, keeping only those `R`-tuples which contain it, 
// removing it from those `R`-tuples, and relabeling the
// remaining elements in an order preserving way.
//
// See also the chirotope-operator `contract_element`, which
// produces chirotopes on the set of elements `0..N-1`, by
// setting the chirotope to `0` on all `R`-tuples which do 
// not contain the specified element. See also `contract_labels`
// for contracting multiple labels simultaneously.
template<int R, int N>
constexpr svo::PullBack<Chirotope<R,N>,Chirotope<R-1,N-1>>
contract_label(int);

// Constructs a chirotope-operator, which contracts the given
// set of elements of any chirotope of rank `R` on set of 
// elements `0..N-1`. In other words, it sets the chirotope 
// to `0` on all `R`-tuples which do not contain all given 
// elemenst, without changing the rank or the set of elements 
// of said chirotope. Note that if the given set of elements 
// is not independent, then the resulting chirotope will be `0`.
//
// See also the chirotope-operator `contract_labels`, which does
// the same, but afterwards forgets all elements in the given 
// set, and relabels while preserving the order of elements,
// thus changing the set of elements of the resulting chirotope.
// See also `contract_element` for contracting only a single
// element at a time.
template<int R, int N, typename Iterable>
constexpr svo::Multiply_P0<Chirotope<R ,N>> 
contract_elements(const Iterable&);

// Constructs a chirotope-operator, which contracts the given
// set of labels of any chirotope of rank `R` on set of 
// elements `0..N-1`. In other words, it forgets all specified
// labels and all `R`-tuples which do not contain all of them,
// removes these labels from the remaining `R`-tuples, and
// relabels the remaining elements in an order preserving manner.
// Note that if the given set of elements is not independent, then 
// the resulting chirotope (which is of smaller rank) will be `0`.
//
// See also the chirotope-operator `contract_elements`, which
// produces chirotopes of rank `R0` on set of elements `0..N0-1`
// by setting the chirotope to `0` on all `R`-tuples which
// do not contain all specified elements. See also 
// `contract_label` for contracting only a single element 
// at a time.
template<int R0, int N0, int R1, int N1, typename Iterable>
constexpr svo::PullBack<Chirotope<R0 ,N0>, Chirotope<R1, N1>> 
contract_labels(const Iterable&);

// =======
// DUALITY
// =======

// Constructs a chirotope-operator, which computes the dual
// of any chirotope of rank `R` on `N` elements. Be warned that
// `dualize<N-R,N>*dualize<R,N>` is the identity if `R*(N-R)`
// is even, and `.inverted()` if it is odd.
template<int R, int N>
constexpr svo::composes_to<svo::PushForward<Chirotope<R,N>,Chirotope<N-R,N>>,svo::Multiply<Chirotope<R,N>>>
dualize{[]() constexpr{
    const auto NcR = binomial_coefficient(N, R);
    std::array<int, NcR> induced{};
    for (int idx = 0; idx < NcR; ++idx) {
        induced[idx] = NcR - 1 - idx;
    }
    return svo::composes_to<svo::PushForward<Chirotope<R,N>,Chirotope<N-R,N>>,svo::Multiply<Chirotope<R,N>>>{
        {induced},
        {bit_vector<NcR>{Chirotope<R,N>::RTUPLES::LIST::is_dual_inverted_mask32}}
    };
}()};

// ===============
// CRYPTOMORPHISMS
// ===============

// An array of signvector-operators, where the `i`th entry
// extracts from a `Chirotope<R,N>` one of its two cocircuits 
// which vanish on the `i`th `R-1`-subset of `0..N-1`. In case
// the specified elements are not independent, the produced
// signvector is constantly zero instead.
template<int R, int N>
constexpr auto cocircuit_extractors_from_chirotope{[]() constexpr {
    using RTUPLES = Chirotope<R, N>::RTUPLES; // R-tuples
    using R1TUPLES = Rtuples::RTUPLES<char,R-1,N,int>; // (R-1)-tuples
    std::array<
        svo::composes_to<
            svo::Multiply_PM0<sign_vector<N>>,
            svo::PullBack<sign_vector<RTUPLES::NR>, sign_vector<N>>
        >,
        R1TUPLES::NR
    > cocircuit_constructors{};
    for (int cocirc_idx = 0; cocirc_idx < R1TUPLES::NR; ++cocirc_idx) {
        std::array<int, N> induced{};
        std::array<char, R-1> R1tuple = R1TUPLES::LIST::array[cocirc_idx];
        bit_vector<N> to_reorient{};
        bit_vector<N> zero{};
        for (int element = 0; element < N; ++element) {
            std::array<char, R> full_Rtuple{};
            char t = 0;
            while (t < R-1 && R1tuple[t] < element) {
                full_Rtuple[t] = R1tuple[t];
                ++t;
            }
            to_reorient.set_bit(element, t % 2);
            if (t < R-1 && R1tuple[t] == element) {
                zero.set_bit(element, true);
                continue;
            }
            full_Rtuple[t] = element;
            for (++t; t < R; ++t) {
                full_Rtuple[t] = R1tuple[t - 1];
            }
            induced[element] = RTUPLES::index_of_ordered(full_Rtuple);
        }
        cocircuit_constructors[cocirc_idx] = {
            {to_reorient, ~zero}, 
            {{induced}}
        };
    }
    return cocircuit_constructors;
}()};

// Computes the list of cocircuits of a chirotope.
// From each pair of cocircuits `C`, `-C` only one is
// included in the list. Those included are in lexicographically
// sorted by the sets on which vanish.
template<int R, int N>
constexpr std::vector<sign_vector<N>> cocircuits(const Chirotope<R, N>&);

}

#include "OMoperations.cpp"