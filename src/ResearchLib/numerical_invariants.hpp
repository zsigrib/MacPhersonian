#pragma once

#include "OMtools.hpp"
#include "research_file_template.hpp"
#include "signvectors.hpp"

// Numerical invariants of oriented matroids.

namespace research {

inline int euler_char_of_sphere_of_dimension(int dimension) {
    return 2 - 2 * (dimension % 2);
}

// The dimension of the sphere which is homotopy equivalent to 
// the (strict) upper cone of a rank 2 oriented matroid. The formula
// is well-defined for any oriented matroid.
template<int R, int N>
int dim_of_upper_cone_in_rk_2(const Chirotope<R, N>& chi) {
    auto cocircs = OM_operations::cocircuits(chi);
    int total_size_of_kernels = 0;
    for (sign_vector<N> cocirc: cocircs) {
        total_size_of_kernels += N - cocirc.count_nonzero();
    }
    return total_size_of_kernels - cocircs.size() * (chi.loopcount() + 1) - 1 + 2 * chi.loopcount();
}

// The dimension of the sphere which is homotopy equivalent to 
// the (strict) upper cone of a rank 1 oriented matroid. The formula
// is well-defined for any oriented matroid.
template<int R, int N>
int dim_of_upper_cone_in_rk_1(const Chirotope<R, N>& chi) {
    return chi.loopcount() - 1;
}

}