#pragma once

#include <float.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include "OMtools.hpp"
#include "abstractly_solvable.hpp"
#include "lowercones.hpp"
#include "researchlib.hpp"
#include "program_template.hpp"
#include "program_utility.hpp"
#include "verboseness.hpp"

namespace programs {

// Call an oriented matroid M 'good' with respect to 
// an element e if all weak insertion problems of the form (?,e,M) 
// are abstractly solvable and are never isolated. This program checks 
// if all (simple) oriented matroids with certain parameters are good
// with respect to some element.
template<int R, int N>
int always_weakly_abstractly_reducible() {
const auto iso_representatives = OMexamples::read_all_Finschi_representatives<R,N>();
std::cout << "Call an oriented matroid M 'good' with respect to "
"an element e if all weak insertion problems of the form (?,e,M) "
"are abstractly solvable and are never isolated. This program checks "
"if all (simple) oriented matroids with certain parameters are good "
"with respect to some element.\n\n";
auto matroids = research::matroids_with_few_loops<R, N>(
    N - research::minimum_N_for_not_abstractly_solvable<R>
);
auto matroids_to_filter_deletion_with = research::matroids_with_few_loops<R, N>(
    N + 1 - research::minimum_N_for_not_abstractly_solvable<R>
);
std::cout << "We iterate over all " << iso_representatives.size() 
<< " isomorphism classes of simple oriented matroids "
"of rank " << R << " and number of elements " << N << ".\n\n";
int idx_of_current_target = 0;
std::vector<int> good_wrt_element;
for (Chirotope<R,N> chi: iso_representatives) {
    // Start announcement
    std::cout << "[TARGET " << idx_of_current_target << "/"
    << iso_representatives.size() << "] " << chi << "\n";
    // Quick setup
    auto cocircuits = OM_operations::cocircuits(chi);
    int loopcount_chi = chi.loopcount();
    good_wrt_element.push_back(-1);
    std::array<sign_vector_operations::Multiply_P0<Chirotope<R, N>>,N> deleters;
    std::vector<bool> not_loop_in_chi;
    for (int e = 0; e < N; ++e) {
        deleters[e] = OM_operations::delete_element<R, N>(e);
        not_loop_in_chi.push_back(!chi.is_loop(e));
    }
    // Slow setup
    std::cout << "Computing lower cone of target... ";
    auto lower_cone_filtered = research::lower_cone_with_few_loops<R,N>(
        chi, 
        matroids,
        N - research::minimum_N_for_not_abstractly_solvable<R>,
        verboseness::result
    );
    // Check each element if chi can be weakly reduced by it
    for (int e = 0; e < N; ++e) {
        if (chi.is_loop(e)) {
            std::cout << "- element " << e << " is a loop\n";
            continue;
        }
        // Test isolatedness
        bool nonisolated_in_all_relevant_deletions = true;
        not_loop_in_chi[e] = false; // see reason 4 lines below
        for (int nr_to_delete = 0; 
        nr_to_delete <= N - research::minimum_N_for_isolation<R> - loopcount_chi
            && nonisolated_in_all_relevant_deletions; 
        ++nr_to_delete) {
            for (auto to_delete: research::tuple_iterator_on_subset(nr_to_delete,N,not_loop_in_chi)) {
                Chirotope<R, N> deletion = chi;
                for (int element_to_delete: to_delete) {
                    deletion = deleters[element_to_delete](deletion);
                }
                auto cocircuits_of_deletion = OM_operations::cocircuits(deletion);
                if (research::is_isolated(
                    deletion, 
                    cocircuits_of_deletion,
                    e
                )) {
                    std::cout << "- element " << e << " is isolated in M\\{";
                    utility::print_comma_separated_iterable_of_ints(to_delete);
                    std::cout << "} = " << deletion << "\n";
                    nonisolated_in_all_relevant_deletions = false;
                    break;
                }
            }
        }
        not_loop_in_chi[e] = true;
        if (!nonisolated_in_all_relevant_deletions) continue;
        std::cout << "- element " << e << " is NOT isolated!\n";
        // Set up abstractly solvability check
        if (research::is_always_abstractly_solvabe(
            chi,
            e,
            lower_cone_filtered,
            matroids_to_filter_deletion_with
        )) {
            good_wrt_element.back() = e;
            break;
        }
    }
    if (good_wrt_element.back() == -1) {
        std::cout << "[FAILURE] The chirotope " << chi << " is not good "
        "with respect to any element.\n\n";
    }
    // Iterate
    idx_of_current_target++;
}

bool was_anything_not_reducible = false;
std::cout << "Result compilation:\n";
for (int idx = 0; idx < iso_representatives.size(); ++idx) {
    if (good_wrt_element[idx] != -1) {
        std::cout << iso_representatives[idx] << " is weakly reducible by "
        << good_wrt_element[idx] << "\n";
    } else {
        std::cout << iso_representatives[idx] << " is not weakly reducible by anything\n";
        was_anything_not_reducible = true;
    }
}
if (was_anything_not_reducible) {
    std::cout << "\nAt least one chirotope was not weakly reducible.\n\n";
} else {
    std::cout << "\nAll chirotopes were weakly reducible.\n\n";
}
std::cout << "Terminating";
return 0;
} 

}