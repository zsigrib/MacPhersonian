#pragma once

#include <vector>
#include <concepts>
#include "OMtools.hpp"
#include "abstractly_solvable.hpp"
#include "research_file_template.hpp"

namespace research {

// Find an element `e` of the oriented matroid `chi` satisfying the
// property below; return `-1` if no such element exists.
//
// Property: all weak insertion problems of the form
// `(?,e,chi)` are both non-isolated and abstractly solvable.
//
// `check_abstractly_solvability` should be a function which tells
// if all weak insertion problems like above are abstractly solvable.
// See other overloads of this functions, where this argument is replaced
// by more concrete parameters.
template<int R, int N,
std::invocable<const Chirotope<R, N>&, int, enum verboseness> Function>
int _weakly_reducible_by(
    const Chirotope<R, N>& chi,
    Function check_abstractly_solvability,
    enum verboseness verbose = verboseness::result
);

// Find an element `e` of the oriented matroid `chi` satisfying the
// property below; return `-1` if no such element exists.
//
// Property: all weak insertion problems of the form
// `(?,e,chi)` are both non-isolated and abstractly solvable.
//
// It is required that `lower_cone_of_chi` contains precisely 
// those weak map images of `chi` which have at least 
// `minimum_N_for_not_abstractly_solvable<R>` many nonloops.
// `matroids_to_filter_deletion_with` must contain all matroids with
// at least `minimum_N_for_not_abstractly_solvable<R> - 1`
// nonloops, sorted into subvectors by basecount (shifted
// down by 1), but may contain more matroids.
template<int R, int N>
int weakly_reducible_by(
    const Chirotope<R, N>& chi,
    const std::vector<Chirotope<R, N>>& lower_cone_of_chi,
    const std::vector<std::vector<Matroid<R, N>>>& matroids_to_filter_deletion_with,
    enum verboseness verbose = verboseness::result
) {
    return _weakly_reducible_by(
        chi,
        [&] (const Chirotope<R, N>& c, int e, enum verboseness v) {
            return is_always_abstractly_solvabe(
                c,
                e,
                lower_cone_of_chi,
                matroids_to_filter_deletion_with,
                v
            );
        },
        verbose
    );
}

// This is the same as a different function of the same name,
// but some input variables are computed automatically.
// `matroids` must contain all matroids with at least
// `minimum_N_for_not_abstractly_solvable<R>` many nonloops.
template<int R, int N>
int weakly_reducible_by(
    const Chirotope<R, N>& chi,
    const std::vector<std::vector<Matroid<R, N>>>& matroids,
    const std::vector<std::vector<Matroid<R, N>>>& matroids_to_filter_deletion_with,
    enum verboseness verbose = verboseness::result
) {
    return _weakly_reducible_by(
        chi,
        [&] (const Chirotope<R, N>& c, int e, enum verboseness v) {
            return is_always_abstractly_solvabe(
                c,
                e,
                matroids,
                matroids_to_filter_deletion_with,
                v
            );
        },
        verbose
    );
}

// This is the same as a different function of the same name,
// but some input variables are computed automatically.
template<int R, int N>
int weakly_reducible_by(
    const Chirotope<R, N>& chi,
    const std::vector<std::vector<Matroid<R, N>>>& matroids_to_filter_deletion_with,
    enum verboseness verbose = verboseness::result
) {
    return _weakly_reducible_by(
        chi,
        [&] (const Chirotope<R, N>& c, int e, enum verboseness v) {
            return is_always_abstractly_solvabe(
                c,
                e,
                matroids_to_filter_deletion_with,
                v
            );
        },
        verbose
    );
}

template<int R, int N>
int weakly_reducible_by(
    const Chirotope<R, N>& chi,
    enum verboseness verbose = verboseness::result
) {
    return _weakly_reducible_by(
        chi,
        is_always_abstractly_solvabe<R, N>,
        verbose
    );
}

}

#include "weakly_reducible_impl.hpp"