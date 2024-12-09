#pragma once
#include <iostream>
#include <vector>
#include "OMtools.hpp"
#include "program_template.hpp"

namespace programs {

// Computes the f-vector of the lower cone of a given oriented
// matroid inside `MacP(R,N)`. In doing so, the f-vectors of all
// weak map images are computed, and those corresponding to
// oriented matroids in `targets` will be printed.
//
// The function accepts as optional input a vector of all weak map
// images of `top` (`wmis_precalculated` signifies if this is given).
template<int R, int N>
int compute_fvector_of_lowercone(
    const Chirotope<R, N>& top,
    const std::vector<std::vector<Chirotope<R, N>>>& targets,
    const std::vector<std::vector<Chirotope<R,N>>>& wmis_by_bases_={},
    bool wmis_precalculated = true,
    enum verboseness verbose = verboseness::checkpoints
) {

for (auto target_list: targets) {
    for (auto target : target_list) {
        if (!top.OM_weak_maps_to(target)) {
            std::cerr << "The list of targets included        " << target
            << ",\nbut this is not a weak map image of " << top;
            return 1;
        }
    }
}

if (!wmis_precalculated && verbose >= verboseness::info) {
    std::cout << "The set of weak map images is not precalculated. "
    "Calculating it...\n\n";
}
const auto& wmis_by_bases(wmis_precalculated? wmis_by_bases_ : research::generate_lower_cone(top));

std::vector<Chirotope<R, N>> all_wmis;
std::vector<int> base_counts;
//std::vector<std::vector<size_t>> smaller_OM_indices;
std::vector<std::array<size_t, binomial_coefficient(N,R)>> lower_cone_face_vectors;
int top_basecount = top.countbases();
for (int b = 0; b < top_basecount; b++) {
    if (verbose >= verboseness::checkpoints) {
        std::cout << "Computing face vectors of OMs with " << b+1 << " bases...\n";
    }
    size_t non_contr = 0;
    for (auto c : wmis_by_bases[b]) {
        // PARSE
        std::vector<size_t> weak_images = smaller_OMs(
            all_wmis,
            base_counts,
            c,
            b+1
        );
        auto f_vector = face_vector<binomial_coefficient(N,R)>(
            lower_cone_face_vectors, weak_images
        );
        auto ec = euler_characteristic<binomial_coefficient(N,R)>(f_vector);
        if (ec != 1) non_contr++;
        // PRINT
        for (auto target : targets[b]) {
            if (target.is_same_OM_as(c) && verbose >= verboseness::result) {
                std::cout << "Found target OM " << target << "! "
                "Euler-characteristic: " << ec << ", face vector: ("
                << f_vector[0];
                for (auto d = 1; d < binomial_coefficient(N,R); d++) {
                    std::cout << ", " << f_vector[d];
                }
                std::cout << ")\n";
            }
        }
        // SAVE
        all_wmis.push_back(c);
        base_counts.push_back(b+1);
        //smaller_OM_indices.push_back(weak_images);
        lower_cone_face_vectors.push_back(f_vector);
    }
    if (verbose >= verboseness::checkpoints) {
        std::cout << "# of non 1 Euler-characteristic lower cones: "
        << non_contr << "/" << wmis_by_bases[b].size() << " (# of bases: "
        << b+1 << ").\n";
    }
}
if (verbose >= verboseness::info) {
    std::cout << "All face vectors computed successfully. Terminating.\n";
}
return 0;

}



// Same as `compute_fvector_of_lowercone`, but it is assumed
// that the only target is `chi`; in particular computes the
// f-vector of the lower cone of `chi` inside `MacP(R,N)`.
template<int R, int N>
int compute_f_vector_of_single_element(
    const Chirotope<R,N>& chi,
    enum verboseness verbose = verboseness::checkpoints
) {

std::vector<std::vector<Chirotope<R,N>>> targets(
    binomial_coefficient(N,R), 
    std::vector<Chirotope<R,N>>()
);
targets[chi.countbases() - 1].push_back(chi);
return compute_fvector_of_lowercone<R,N>(
    chi,
    targets,
    std::vector<std::vector<Chirotope<R,N>>>{},
    false,
    verbose
);

}



// Computes the f-vector of the lower cones of the given
// target oriented matroids inside `MacP(R,N)`. Faster than
// calling `compute_f_vector_of_single_element` multiple times
// as this one only reads in the set of matroids once. No
// relation is assumed between the given oriented matroids.
template<int R, int N>
int compute_f_vector_of_independent_elements(
    const std::vector<Chirotope<R, N>>& targets,
    enum verboseness verbose = verboseness::checkpoints
) {

auto matroids = research::read_matroids<R, N>(verbose);

int target_idx = 0;
for (auto chi: targets) {
    std::cout << "\n[TARGET " << target_idx+1 << "/" 
    << targets.size() << "]\n";
    std::vector<std::vector<Chirotope<R,N>>> targets_(
        binomial_coefficient(N,R), 
        std::vector<Chirotope<R,N>>()
    );
    targets_[chi.countbases() - 1].push_back(chi);
    auto lc = research::generate_lower_cone(
        chi, 
        matroids, 
        std::min(verboseness::result, verbose)
    );
    std::cout << "\n";
    auto ret = compute_fvector_of_lowercone(
        chi,
        targets_,
        lc,
        true,
        std::min(verboseness::result, verbose)
    );
    if (ret) return ret;
    target_idx++;
}
std::cout << "\nFinished with all targets.";
return 0;

}

}