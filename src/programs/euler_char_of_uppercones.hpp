#pragma once

#include <cstddef>
#include <iostream>
#include "OMtools.hpp"
#include "mymath.hpp"
#include "ordercomplexes.hpp"
#include "program_template.hpp"
#include "euler_char_of_lowercones.hpp"

namespace programs {

// Using a database of all oriented matroids of a given rank and
// number of elements, compute the Euler characteristic of the
// upper cone of each of them. Compare this to the Euler
// characteristic predicted by the function `predictor`.
// This must be a callable with two parameters, a chirotope and
// the resulting Euler characteristic, which returns a bool
// describing whether this Euler characteristic was expected or not.
template<int R, int N, typename Predictor>
inline int compute_euler_char_of_upper_cones_using_database(Predictor predictor) {
static_assert((R == 3 && N == 6) || (R == 3 && N == 7),
"This program must be compiled with parameters (3,6) or (3,7)!");

std::cout << "This program reads in the database of rank " << R << " oriented "
"matroids with " << N << " elements, and for each oriented matroid computes the "
"Euler characteristic of its upper cone, i.e. of the poset of oriented matroids "
"which are strictly greater than it.\n\n";

auto input = ReadOMDataFromFiles<Chirotope<R,N>>(
    &database_names::OM_set<R,N>,
    6
);

std::vector<Chirotope<R,N>> all_OMs;
std::vector<int> base_counts;
int current_basecount = 1;
for (auto p: input) {
    if (p.first > current_basecount) {
        std::cout << "[" << current_basecount << "] Finished reading OMs with "
        << current_basecount << " bases.\n";
        current_basecount = p.first;
    }
    all_OMs.push_back(p.second);
    base_counts.push_back(p.first);
}
std::cout << "[" << current_basecount << "] Finished reading OMs with " 
<< current_basecount << " bases.\n\n";
std::cout << "Now we compute the face vectors of the upper cones.\n\n";

std::vector<std::array<size_t, binomial_coefficient(N,R)>> upper_cone_face_vectors(all_OMs.size());
EulerCharAnalyzer ec_analyzer;
ec_analyzer.these_are_ecs_of = "upper cones";
EulerCharAnalyzer strange_ec_analyzer;
strange_ec_analyzer.these_are_ecs_of = "upper cones";
for (size_t idx = all_OMs.size() - 1; idx >= 0; --idx) {
    auto larger_OMs = bigger_OMs(
        all_OMs, 
        base_counts, 
        all_OMs[idx],
        all_OMs[idx].countbases()
    );

    auto fvector = face_vector<binomial_coefficient(N, R)>(
        upper_cone_face_vectors,
        larger_OMs
    );
    // Print message:
    if (base_counts[idx] < current_basecount) {
        std::cout << "[" << current_basecount << "] Finished computing upper cones for OMs with "
        << current_basecount << " bases. ";
        ec_analyzer.end_batch();
        std::cout << "Some of these might also belong to a batch of OMs for whom we could not "
        "correctly predict the Euler characteristic. ";
        strange_ec_analyzer.end_batch();
        current_basecount = base_counts[idx];
    }
    // Continue ec counting:
    auto ec = euler_characteristic<binomial_coefficient(N, R)>(fvector);
    ec_analyzer.add_entry(ec);
    if (!predictor(all_OMs[idx], ec)) {
        strange_ec_analyzer.add_entry(ec);
    }
    // Save results:
    upper_cone_face_vectors[idx] = fvector;
}
std::cout << "[" << current_basecount << "] Finished computing upper cones for OMs with "
<< current_basecount << " bases. ";
ec_analyzer.end_batch();
std::cout << "Some of these might also belong to a batch of OMs for whom we could not "
"correctly predict the Euler characteristic. ";
strange_ec_analyzer.end_batch();
return 0;
}

// Using a database of all oriented matroids of a given rank and
// number of elements, compute the Euler characteristic of the
// upper cone of each of them. 
template<int R, int N>
inline int compute_euler_char_of_upper_cones_using_database() {
static_assert((R == 3 && N == 6) || (R == 3 && N == 7),
"This program must be compiled with parameters (3,6) or (3,7)!");

std::cout << "This program reads in the database of rank " << R << " oriented "
"matroids with " << N << " elements, and for each oriented matroid computes the "
"Euler characteristic of its upper cone, i.e. of the poset of oriented matroids "
"which are strictly greater than it.\n\n";

auto input = ReadOMDataFromFiles<Chirotope<R,N>>(
    &database_names::OM_set<R,N>,
    6
);

std::vector<Chirotope<R,N>> all_OMs;
std::vector<int> base_counts;
int current_basecount = 1;
for (auto p: input) {
    if (p.first > current_basecount) {
        std::cout << "[" << current_basecount << "] Finished reading OMs with "
        << current_basecount << " bases.\n";
        current_basecount = p.first;
    }
    all_OMs.push_back(p.second);
    base_counts.push_back(p.first);
}
std::cout << "[" << current_basecount << "] Finished reading OMs with " 
<< current_basecount << " bases.\n\n";
std::cout << "Now we compute the face vectors of the upper cones.\n\n";

std::vector<std::array<size_t, binomial_coefficient(N,R)>> upper_cone_face_vectors(all_OMs.size());
EulerCharAnalyzer ec_analyzer;
ec_analyzer.these_are_ecs_of = "upper cones";
for (long long idx = all_OMs.size() - 1; idx >= 0; --idx) {
    auto larger_OMs = bigger_OMs(
        all_OMs, 
        base_counts, 
        all_OMs[idx],
        all_OMs[idx].countbases()
    );

    auto fvector = face_vector<binomial_coefficient(N, R)>(
        upper_cone_face_vectors,
        larger_OMs
    );
    // Print message:
    if (base_counts[idx] < current_basecount) {
        std::cout << "[" << current_basecount << "] Finished computing upper cones for OMs with "
        << current_basecount << " bases. ";
        ec_analyzer.end_batch();
        current_basecount = base_counts[idx];
    }
    // Continue ec counting:
    auto ec = euler_characteristic<binomial_coefficient(N, R)>(fvector);
    ec_analyzer.add_entry(ec);
    // Save results:
    upper_cone_face_vectors[idx] = fvector;
}
std::cout << "[" << current_basecount << "] Finished computing upper cones for OMs with "
<< current_basecount << " bases. ";
ec_analyzer.end_batch();
return 0;
}

}