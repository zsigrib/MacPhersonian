#pragma once

#include <vector>
#include "OMtools.hpp"
#include "lowercones.hpp"
#include "research_file_template.hpp"
#include "verboseness.hpp"

namespace research {

template<int R>
constexpr int minimum_N_for_not_abstractly_solvable = R + 3;

template<>
constexpr int minimum_N_for_not_abstractly_solvable<3> = 7; //NOLINT

// Determine if all weak insertion problems of the form
// `(?,element,chi)` are abstractly solvable. It is
// required that `lower_cone_of_chi` contains precisely 
// those weak map images of `chi` which have at least 
// `minimum_N_for_not_abstractly_solvable<R>` many nonloops.
// `matroids_to_filter_deletion_with` must contain all matroids with
// at least `minimum_N_for_not_abstractly_solvable<R> - 1`
// nonloops, sorted into subvectors by basecount (shifted
// down by 1), but may contain more matroids.
template<int R, int N>
bool is_always_abstractly_solvabe(
    const Chirotope<R, N>& chi,
    int element,
    const std::vector<Chirotope<R, N>>& lower_cone_of_chi,
    const std::vector<std::vector<Matroid<R, N>>>& matroids_to_filter_deletion_with,
    enum verboseness verbose = verboseness::result
);

// This is the same as a different function of the same name,
// but some input variables are computed automatically.
// `matroids` must contain all matroids with at least
// `minimum_N_for_not_abstractly_solvable<R>` many nonloops.
template<int R, int N>
bool is_always_abstractly_solvabe(
    const Chirotope<R, N>& chi,
    int element,
    const std::vector<std::vector<Matroid<R, N>>>& matroids,
    const std::vector<std::vector<Matroid<R, N>>>& matroids_to_filter_deletion_with,
    enum verboseness verbose = verboseness::result
) {
    return is_always_abstractly_solvabe(
        chi,
        element,
        lower_cone_with_few_loops(
            chi, 
            matroids,
            N - minimum_N_for_not_abstractly_solvable<R>,
            verbose
        ),
        matroids_to_filter_deletion_with,
        verbose
    );
}

// This is the same as a different function of the same name,
// but some input variables are computed automatically.
template<int R, int N>
bool is_always_abstractly_solvabe(
    const Chirotope<R, N>& chi,
    int element,
    const std::vector<std::vector<Matroid<R, N>>>& matroids_to_filter_deletion_with,
    enum verboseness verbose = verboseness::result
) {
    return is_always_abstractly_solvabe(
        chi,
        element,
        matroids_to_filter_deletion_with,
        matroids_to_filter_deletion_with,
        verbose
    );
}

// Determine if all weak insertion problems of the form
// `(?,element,chi)` are abstractly solvable.
template<int R, int N>
bool is_always_abstractly_solvabe(
    const Chirotope<R, N>& chi,
    int element,
    enum verboseness verbose = verboseness::result
) {
    auto matroids_to_filter_deletion_with = matroids_with_few_loops<R, N>(
        N + 1 - minimum_N_for_not_abstractly_solvable<R>,
        verboseness::silent
    );
    return is_always_abstractly_solvabe(
        chi,
        element,
        matroids_to_filter_deletion_with,
        verbose
    );
}

}

#include "abstractly_solvable_impl.hpp"
