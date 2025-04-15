#pragma once

#include <concepts>
#include <iostream>
#include <vector>
#include "OMtools.hpp"
#include "abstractly_solvable.hpp"
#include "isolation.hpp"
#include "moreRtuples.hpp"
#include "program_utility.hpp"
#include "research_file_template.hpp"

namespace research {

// COMMENT: The implementations are all the same, but with a different call
// to is_always_abstractly_reducible

template<int R, int N, 
std::invocable<const Chirotope<R, N>&, int, enum verboseness> Function>
int _weakly_reducible_by(
    const Chirotope<R, N>& chi,
    Function check_abstractly_solvability,
    enum verboseness verbose
) {
    // Start announcement
    if (verbose >= verboseness::info) {
        std::cout << "Looking for an element by which the target chirotope  " << chi 
        << " is weakly reducible, i.e. all weak insertion problems are non-isolated "
        "and abstractly solvable...\n";
    }
    // Quick setup
    auto cocircuits = OM_operations::cocircuits(chi);
    int loopcount_chi = chi.loopcount();
    std::array<sign_vector_operations::Multiply_P0<Chirotope<R, N>>,N> deleters;
    std::vector<bool> not_loop_in_chi;
    for (int e = 0; e < N; ++e) {
        deleters[e] = OM_operations::delete_element<R, N>(e);
        not_loop_in_chi.push_back(!chi.is_loop(e));
    }
    // Check each element if chi can be weakly reduced by it
    for (int e = 0; e < N; ++e) {
        if (chi.is_loop(e)) {
            if (verbose >= verboseness::info) {
                std::cout << "- element " << e << " is a loop\n";
            }
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
                    if (verbose >= verboseness::info) {
                        std::cout << "- element " << e << " is isolated in M\\{";
                        programs::utility::print_comma_separated_iterable_of_ints(to_delete);
                        std::cout << "} = " << deletion << "\n";
                    }
                    nonisolated_in_all_relevant_deletions = false;
                    break;
                }
            }
        }
        not_loop_in_chi[e] = true;
        if (!nonisolated_in_all_relevant_deletions) continue;
        if (verbose >= verboseness::info) {
            std::cout << "- element " << e << " is NOT isolated!\n";
        }
        // Set up abstractly solvability check
        if (check_abstractly_solvability(
            chi,
            e,
            std::min(verbose, verboseness::result)
        )) {
            if (verbose >= verboseness::result) {
                std::cout << "          and non-isolated, therefore M is weakly reducible by "
                << e << ".\n";
            }
            return e;
        }
    }
    if (verbose >= verboseness::result) {
        std::cout << "[FAILURE] The chirotope " << chi << " is not weakly reducible "
        "by any element.\n\n";
    }
    return -1;
}

}