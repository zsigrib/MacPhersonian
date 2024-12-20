#pragma once

#include <iostream>
#include <vector>
#include "OMtools.hpp"
#include "program_template.hpp"

namespace programs {

// Computes the euler characteristic of the lower
// cone of the chirotope `RIN(9)` within `MacP(3,9)`
// modulo 3. For this, a database of all oriented
// matroids in `MacP(3,9)` is used which contains
// precisely those which are fixed under an action of
// `Z_3`. The total face vector of the complex of such
// oriented matroids is also computed.
inline int compute_euler_char_of_JRG_mod_3()
{
constexpr int R0 = 3;
constexpr int N0 = 9;
std::vector<std::vector<Chirotope<R0,N0>>> fixed_OMs_by_bases
(binomial_coefficient(N0, R0), std::vector<Chirotope<R0,N0>>());
auto input = ReadOMDataFromFile<Chirotope<R0,N0>>(
    "../../../resources/fixed_oriented_matroid_sets/fixed_OMs_rank3_9elements_group_Z3.txt",
    6
);
int count = 0;
for (auto p: input) {
    int base_count = p.countbases();
    fixed_OMs_by_bases[base_count - 1].push_back(p);
    count++;
}
std::vector<Chirotope<R0,N0>> all_fixed_OMs;
std::vector<int> base_counts;
std::vector<std::vector<size_t>> smaller_OM_indices;
std::vector<std::array<size_t, binomial_coefficient(N0,R0)>> lower_cone_face_vectors;
for (int b = 0; b < binomial_coefficient(N0, R0); b++) {
    std::cout << "Scanning OMs with " << b+1 << " bases...\n";
    size_t non_contr = 0;
    for (auto c : fixed_OMs_by_bases[b]) {
        std::vector<size_t> weak_images = smaller_OMs(
            all_fixed_OMs,
            base_counts,
            c,
            b+1
        );
        auto f_vector = face_vector<binomial_coefficient(N0,R0)>(
            lower_cone_face_vectors, weak_images
        );
        auto ec = euler_characteristic<binomial_coefficient(N0,R0)>(f_vector);
        if (ec % 3 != 1 && ec % 3 != -2) non_contr++;
        // Save OM to list of all OMS:
        all_fixed_OMs.push_back(c);
        base_counts.push_back(b+1);
        smaller_OM_indices.push_back(weak_images);
        lower_cone_face_vectors.push_back(f_vector);
        if (c == OMexamples::RIN9) {
            std::cout << "Euler characteristic of RIN9: " << euler_characteristic<binomial_coefficient(N0,R0)>(f_vector) 
            << ", mod 3:" <<  euler_characteristic<binomial_coefficient(N0,R0)>(f_vector) % 3 << "\n";
        }
    }
    std::cout << "# of non 1 (mod 3) Euler-characteristic lower cones: "
    << non_contr << "/" << fixed_OMs_by_bases[b].size() << " (# of bases: "
    << b+1 <<").\n";
}
auto f = face_vector<binomial_coefficient(N0,R0)>(lower_cone_face_vectors);
std::cout << "Total f-vector: (" << f[0];
for (auto d = 1; d < binomial_coefficient(N0,R0); d++) {
    std::cout << ", " << f[d];
}
std::cout << ")";
return 0;

}

}