#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include "OMtools.hpp"
#include "lowercones.hpp"
#include "isolation.hpp"
#include "research_file_template.hpp"
#include "abstractly_solvable.hpp"

namespace research {

template<int R, int N>
bool is_always_abstractly_solvabe(
    const Chirotope<R, N>& chi,
    int element,
    const std::vector<Chirotope<R, N>>& lower_cone_of_chi,
    const std::vector<std::vector<Matroid<R, N>>>& matroids_to_filter_deletion_with,
    enum verboseness verbose
) {
    std::array<bool, N> is_loop_in_chi;
    for (int f = 0; f < N; ++f) {
        is_loop_in_chi[f] = chi.is_loop(f);
    }
    auto delete_e = OM_operations::delete_element<R,N>(element);
    Chirotope<R,N> deletion = delete_e(chi);
    if (verbose >= verboseness::info) {
        std::cout << "Checking that all (loopfree) weak insertion problems"
        " of the form (M', " << element << ", M) with M = " << chi << " are abstractly solvable.\n";
        std::cout << "- computing lower cone of M\\" << element << " = " << deletion << "...\n";
    }
    if (verbose >= verboseness::result) {
        std::cout << "For M\\" << element << ", ";
    }
    auto lc_of_deletion = research::lower_cone_with_few_loops(
        deletion,
        matroids_to_filter_deletion_with,
        N + 1 - minimum_N_for_not_abstractly_solvable<R>, 
        std::min(verbose, verboseness::result)
    );
    if (verbose >= verboseness::info) {
        std::cout << "- sorting them...\n";
    }
    auto compare = [](const Chirotope<R,N>& l, const Chirotope<R,N>& r){
        for (auto i = 0; i < Chirotope<R,N>::NR_INT32; ++i) {
            if ((l.plus[i] | l.minus[i]) == (r.plus[i] | r.minus[i]))
                continue;
            return (l.plus[i] | l.minus[i]) < (r.plus[i] | r.minus[i]);
        }
        return false;
    };
    std::sort(
        lc_of_deletion.begin(),
        lc_of_deletion.end(),
        compare
    );
    if (verbose >= verboseness::info) {
        std::cout << "- deleting " << element << " from lower cone of M, "
        "seeing what is hit in the lower cone of M\\" << element << "...\n";
    }
    std::vector<bool> was_hit(lc_of_deletion.size(),false);
    for (Chirotope<R,N> wmi: lower_cone_of_chi) {
        const Chirotope<R,N> wmi_minus_e = delete_e(wmi);
        if (wmi_minus_e.is_zero()) continue; // element was a coloop!
        // Binary search
        int lower_bound = 0;
        int upper_bound = lc_of_deletion.size() - 1;
        int midpoint = (lower_bound + upper_bound) / 2;
        while (midpoint != upper_bound) {
            if (compare(lc_of_deletion[midpoint], wmi_minus_e)) {
                lower_bound = midpoint + 1;
            } else {
                upper_bound = midpoint;
            }
            midpoint = (lower_bound + upper_bound) / 2;
        }
        if (lc_of_deletion[midpoint] == wmi_minus_e) {
            was_hit[midpoint] = true;
        }
    } 
    int not_hit_idx = -1;
    for (int idx = 0; idx < was_hit.size(); ++idx) {
        if (!was_hit[idx]) {
            not_hit_idx = idx;
            break;
        }
    }
    if (verbose >= verboseness::result) {
        if (not_hit_idx == -1) {
            std::cout << "[SUCCESS] All weak insertion problems of the form (M', "
            << element << ", M) with M = " << chi << " are abstractly solvable\n";
        } else {
            std::cout << "The chirotope " << lc_of_deletion[not_hit_idx]
            << " has no single element extension which is\nsmaller than  " << chi << ".\n";
        }
    }
    return not_hit_idx == -1;
}

}