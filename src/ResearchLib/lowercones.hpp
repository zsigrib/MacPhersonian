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
// set of matroids. Both `matroids` and the output are grouped
// by base count, and indices are shifted by 1 from the actual
// base count.
//
// This uses a database of all matroids.
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



// ===================================
// EASY-TO-USE WRAPPERS FOR OM READING
// ===================================

// Compute the lower cone of a single chirotope. If a database of
// all oriented matroids is available, the non-weak-map-images will
// be filtered out of it, otherwise this lower cone will be computed
// using a database of all matroids.
template<int R, int N>
std::vector<std::vector<Chirotope<R, N>>> lower_cone(
    const Chirotope<R, N>& top,
    enum verboseness verbose = verboseness::result
);

// Compute the lower cone of a single chirotope. If a database of
// all oriented matroids is available, the non-weak-map-images will
// be filtered out of it, otherwise this lower cone will be computed
// using a database of all matroids, which should be provided as an
// argument, with `matroids[b]` containing all matroids of rank `R`
// and number of elements `N` which have `b+1` many bases.
template<int R, int N>
std::vector<std::vector<Chirotope<R, N>>> lower_cone(
    const Chirotope<R, N>& top,
    const std::vector<std::vector<Matroid<R, N>>>& matroids,
    enum verboseness verbose = verboseness::result
);



// ====================
// FILTERED LOWER CONES
// ====================

// Compute a list of all loopfree chirotopes in the lower cone
// of a given chirotope. See also `lower_cone(...)`.
template<int R, int N>
std::vector<Chirotope<R, N>> loopfree_lower_cone(
    const Chirotope<R, N>& top,
    enum verboseness verbose = verboseness::result
);

// Compute a list of all loopfree chirotopes in the lower cone
// of a given chirotope. See also `lower_cone(...)`.
template<int R, int N>
std::vector<Chirotope<R, N>> loopfree_lower_cone(
    const Chirotope<R, N>& top,
    const std::vector<std::vector<Matroid<R, N>>>& matroids,
    enum verboseness verbose = verboseness::result
);

}

#include "lowercones_impl.hpp"
