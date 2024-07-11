#pragma once

#include <vector>
#include <array>
#include "OMs.hpp"

// Computes the Euler-characteristic associated to an f-vector.
template<int dim_bound>
long long euler_characteristic(const std::array<size_t, dim_bound>& f_vector) {
    long long c = 0;
    int sgn = 1;
    for (size_t f: f_vector) {
        c += sgn * f;
        sgn *= -1;
    }
    return c;
}

// Returns the (increasing) list of indices `i` for which
// `less_than(vector[i], bound) == true`.
template<typename Comparable, bool (*less_than)(const Comparable&, const Comparable&)>
std::vector<size_t> smaller_elements(
    const std::vector<Comparable>& vector, 
    const Comparable& bound
) {
    std::vector<size_t> indices;
    size_t i = 0;
    for (auto e: vector) {
        if (less_than(e, bound)) indices.push_back(i);
        i++;
    }
    return indices;
}

// This is a specialized version of `smaller_elements` for
// a list of chirotopes, where the base counts of the chirotopes
// is also known, and the chirotopes are supposed to be ordered
// by them. This means the `for` loop can terminate earlier.
template<int R, int N>
std::vector<size_t> smaller_OMs(
    const std::vector<Chirotope<R, N>>& previous_OMs,
    const std::vector<int>& previous_base_counts,
    const Chirotope<R, N>& bound,
    int base_count
) {
    std::vector<size_t> indices;
    for (size_t id = 0; id < previous_OMs.size(); id++) {
        if (bound.OM_weak_maps_to(previous_OMs[id])) indices.push_back(id);
        else if (previous_base_counts[id] == base_count) break;
    }
    return indices;
}

// Given a partially ordered set `P`, for each element the face 
// vector of the order complex of its strict lower cone, and a 
// lower set `S` of this poset, this function computes the face 
// vector of the order complex of `S`.
// 
// The elements of the poset are labelled `0..N-1`, and
// - `face_vectors_of_lower_cones[i]` is the face vector of the
//   lower cone of the `i`th element of `P`,
// - `indices_of_elements` lists the elements of `S`.
//
// All face vectors are stored as `std::array<size_t, dim_bound>`,
// where `f[i]` is the count of `i`-dimensional simplices. Simplices
// which are higher dimensional than `dim_bound - 1` are ignored. 
//
// Note: `face_vectors_of_lower_cones` can be computed recursively
// using this function.
template<int dim_bound>
std::array<size_t, dim_bound> face_vector(
    const std::vector<std::array<size_t, dim_bound>>& face_vectors_of_lower_cones,
    const std::vector<size_t>& indices_of_elements
) {
    static_assert(dim_bound > 0, "The number of entries allocated to store the face vectors must be positive!");
    std::array<size_t, dim_bound> f; f.fill(0);
    f[0] = indices_of_elements.size();
    for (auto idx : indices_of_elements) {
        for (auto d = 0; d < dim_bound - 1; d++) {
            f[d+1] += face_vectors_of_lower_cones[idx][d];
        }
    }
    return f;
}

// Identical to `face_vector(face_vectors_of_lower_cones, /* a vector storing 0..N-1 */)`.
template<int dim_bound>
std::array<size_t, dim_bound> face_vector(
    const std::vector<std::array<size_t, dim_bound>>& face_vectors_of_lower_cones
) {
static_assert(dim_bound > 0, "The number of entries allocated to store the face vectors must be positive!");
    std::array<size_t, dim_bound> f; f.fill(0);
    f[0] = face_vectors_of_lower_cones.size();
    for (auto idx = 0; idx < face_vectors_of_lower_cones.size(); idx++) {
        for (auto d = 0; d < dim_bound - 1; d++) {
            f[d+1] += face_vectors_of_lower_cones[idx][d];
        }
    }
    return f;
}
