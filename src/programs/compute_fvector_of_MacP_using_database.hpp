#pragma once
#include <iostream>
#include <vector>
#include <array>
#include "OMtools.hpp"
#include "program_template.hpp"

namespace programs {

// Using a database of all oriented matroids of a given
// rank and number of elements, compute the f-vector of
// `MacP(R,N)`, as well as the euler characteristic of
// all subcomplexes defined as lower cones of certain
// oriented matroids.
template<int R, int N>
int compute_fvector_of_MacP_using_database()
{
static_assert((R == 3 && N == 6) || (R == 3 && N == 7),
"This program must be compiled with parameters (3,6) or (3,7)!");
auto input = ReadOMDataFromFiles<Chirotope<R,N>>(
    &database_names::OM_set<R,N>,
    6
);

std::vector<Chirotope<R,N>> all_OMs;
std::vector<int> base_counts;
std::vector<std::vector<size_t>> smaller_OM_indices;
std::vector<std::array<size_t, binomial_coefficient(N,R)>> lower_cone_face_vectors;
size_t id = 0;
int current_basecount = 1;
size_t OMs_with_fixed_basecount = 0;
size_t OMs_with_fixed_basecount_and_good_ec = 0;
size_t total_good_ec = 0;
for (auto p: input) {
    // Parse new OM:
    auto weak_images = smaller_OMs(
        all_OMs,
        base_counts,
        p.second,
        p.first
    );
    auto fvector = face_vector<binomial_coefficient(N,R)>(
        lower_cone_face_vectors,
        weak_images
    );
    // Print message:
    if (p.first > current_basecount) {
        std::cout << "Finished parsing OMs with " << current_basecount << " bases. There were "
        << OMs_with_fixed_basecount_and_good_ec << "/" << OMs_with_fixed_basecount 
        << " many of them with non-1 Euler-characteristic.\n";
        total_good_ec += OMs_with_fixed_basecount_and_good_ec;
        OMs_with_fixed_basecount = 0;
        OMs_with_fixed_basecount_and_good_ec = 0;
        current_basecount = p.first;
    }
    // Continue ec counting:
    auto ec = euler_characteristic<binomial_coefficient(N,R)>(fvector);
    if (ec != 1) OMs_with_fixed_basecount_and_good_ec++;
    // Save results:
    all_OMs.push_back(p.second);
    base_counts.push_back(p.first);
    smaller_OM_indices.push_back(weak_images);
    lower_cone_face_vectors.push_back(fvector);
    // Increment
    OMs_with_fixed_basecount++;
    id++;
}
std::cout << "Finished parsing OMs with " << current_basecount 
<< " bases. There were " << OMs_with_fixed_basecount_and_good_ec << "/" 
<< OMs_with_fixed_basecount << " many of them with non-1 Euler-characteristic.\n";
total_good_ec += OMs_with_fixed_basecount_and_good_ec;

std::cout << "There were " << total_good_ec << "/" << all_OMs.size()
<< " many non-1 Euler characteristic lower cones.\n";

auto f = face_vector<binomial_coefficient(N,R)>(lower_cone_face_vectors);
std::cout << "Total f-vector: (" << f[0];
for (auto d = 1; d < binomial_coefficient(N,R); d++) {
    std::cout << ", " << f[d];
}
std::cout << ")";
return 0;
}

}
