#pragma once

#include <vector>
#include "OMtools.hpp"
#include "research_file_template.hpp"

namespace research {

// ===============================
// MAIN DATABASE READING FUNCTIONS
// ===============================

// Returns the list of all matroids, grouped by basecount.
// The index is shifted by 1 from the actual basecount.
//
// To do so, a database of all oriented matroids is read.
template<int R, int N>
std::vector<std::vector<Matroid<R, N>>> read_matroids(
    enum verboseness verbose = verboseness::checkpoints
);

// Return the set of all weak map images of `top`, grouped by basecount.
// The index is shifted from the basecount by 1.
//
// This uses a database of all matroids.
template<int R, int N>
std::vector<std::vector<Chirotope<R,N>>> generate_lower_cone(
    const Chirotope<R, N>& top, 
    enum verboseness verbose = verboseness::checkpoints
);

// Calculate all weak map images of `top` using a precomputed
// set of matroids - only those weak map images will be returned
// whose underlying matroid is in the given set. Both `matroids`
// and the output are grouped by base count, and indices are 
// shifted by 1 from the actual base count.
template<int R, int N>
std::vector<std::vector<Chirotope<R, N>>> generate_lower_cone(
    const Chirotope<R, N>& top, 
    const std::vector<std::vector<Matroid<R, N>>>& matroids,
    enum verboseness verbose = verboseness::checkpoints
);

// Reads all OMs of the given rank and number of elements, and selects
// only those which fall inside the given lower cone.
//
// This uses a database of all oriented matroids.
template<int R, int N>
std::vector<std::vector<Chirotope<R, N>>> filter_lower_cone(
    const Chirotope<R, N>& top,
    enum verboseness verbose = verboseness::checkpoints
);



// ====================
// FILTERED LOWER CONES
// ====================

template<int R, int N>
bool loopcount_at_most(Chirotope<R, N> chi, int bound) {
    return chi.loopcount(bound + 1) <= bound;
}

template<int R, int N>
bool loopcount_at_most(Matroid<R, N> matroid, int bound) {
    return matroid.loopcount(bound + 1) <= bound;
}

// Compute the list of all chirotopes in the lower cone
// of a given chirotope (i.e. all weak map images) which 
// have at most `max_nr_of_loops` many loops in them.
template<int R, int N>
std::vector<Chirotope<R, N>> lower_cone_with_few_loops(
    const Chirotope<R, N>& top,
    int max_nr_of_loops = 0,
    enum verboseness verbose = verboseness::result
);

// Compute the list of all chirotopes in the lower cone
// of a given chirotope (i.e. all weak map images) which 
// have at most `max_nr_of_loops` many loops in them.
// `matroids` must contain all matroids with at most
// `max_nr_of_loops` many loops, sorted into subvectors
// by basecount (shifted down by 1), but may contain more
// matroids.
template<int R, int N>
std::vector<Chirotope<R, N>> lower_cone_with_few_loops(
    const Chirotope<R, N>& top,
    const std::vector<std::vector<Matroid<R, N>>>& matroids,
    int max_nr_of_loops = 0,
    enum verboseness verbose = verboseness::result
);

// Read in the list of all matroids which have at most 
// `max_nr_of_loops` many loops in them.
template<int R, int N>
std::vector<std::vector<Matroid<R, N>>> matroids_with_few_loops(
    int max_nr_of_loops = 0,
    enum verboseness verbose = verboseness::result
);

}

#include "lowercones_impl.hpp"
