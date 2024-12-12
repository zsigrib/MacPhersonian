#pragma once
#include <iostream>
#include <vector>
#include <algorithm>
#include "OMtools.hpp"
#include "researchlib.hpp"
#include "program_template.hpp"

namespace programs {

// This program proves the conjecture for `R` = 3 and 
// `N` = 7, by checking that all (isomorphism classes) 
// of simple oriented matroids have an element which is 
// non-isolated, and for which all weak insertion problems
// (where the smaller OM has no unnecessary loops) are
// abstractly solvable.
int prove_r3n7() {
const auto r3n7_representatives = OMexamples::read_all_Finschi_representatives<3,7>();
std::cout << "Proving the conjecture for (3,7).\n"
"We iterate over all " << r3n7_representatives.size() 
<< " isomorphism classes of simple oriented matroids "
"of rank 3 and number of elements 7.\n\n";
int idx_of_current_target = 0;
std::vector<int> weakly_reducible_by_element;
for (Chirotope<3,7> chi: r3n7_representatives) {
    // Start announcement
    std::cout << "[TARGET " << idx_of_current_target << "/"
    << r3n7_representatives.size() << "] " << chi << "\n";
    // Quick setup
    auto cocircuits = OM_operations::cocircuits(chi);
    // Slow setup
    std::cout << "Computing lower cone of target...\n";
    auto lower_cone_filtered = research::loopfree_lower_cone<3,7>(chi, verboseness::result);
    // Check each element if chi can be weakly reduced by it
    weakly_reducible_by_element.push_back(-1);
    for (int e = 0; e < 7; ++e) {
        if (research::is_isolated(chi,cocircuits,e)) {
            std::cout << "- element " << e << " is isolated\n";
            continue;
        }
        std::cout << "- element " << e << " is NOT isolated!\n";
        // Set up abstractly solvability check
        auto delete_e = OM_operations::delete_element<3,7>(e);
        Chirotope<3,7> deletion =  delete_e(chi);
        std::cout << "    - computing lower cone of M/" << e << " = " << deletion << "...\n";
        auto lc_of_deletion_unfiltered = research::lower_cone<3,7>(
            deletion, verboseness::silent
        );
        std::cout << "    - keeping only loopfree...";
        std::vector<Chirotope<3,7>> lc_of_deletion_filtered;
        int total_size = 0;
        for (auto wmis_with_fix_basecount: lc_of_deletion_unfiltered) {
            for (Chirotope<3,7> wmi: wmis_with_fix_basecount) {
                bool loopfree_besides_e = true;
                for (int f = 0; f < 7; ++f) {
                    if (f != e && wmi.is_loop(f)) {
                        loopfree_besides_e = false;
                    }
                }
                if (loopfree_besides_e) {
                    lc_of_deletion_filtered.push_back(wmi);
                }
            }
            total_size += wmis_with_fix_basecount.size();
        }
        std::cout << " kept " << lc_of_deletion_filtered.size() << "/"
        << total_size << "\n";
        std::cout << "    - sorting them...\n";
        auto compare = [](const Chirotope<3,7>& l, const Chirotope<3,7>& r){
            for (auto i = 0; i < Chirotope<3,7>::NR_INT32; ++i) {
                if ((l.plus[i] | l.minus[i]) == (r.plus[i] | r.minus[i]))
                    continue;
                return (l.plus[i] | l.minus[i]) < (r.plus[i] | r.minus[i]);
            }
            return false;
        };
        std::sort(
            lc_of_deletion_filtered.begin(),
            lc_of_deletion_filtered.end(),
            compare
        );
        // Delete e from lower cone
        std::vector<bool> was_hit(lc_of_deletion_filtered.size(),false);
        std::cout << "    - deleting " << e << " from lower cone of M, "
        "seeing what is hit in the lower cone of M/" << e << "...\n";
        for (Chirotope<3,7> wmi: lower_cone_filtered) {
            if (wmi == Chirotope<3,7>("++++++++++++0000++++++000+++0000000")) {
                auto xyz = 7;
            }
            const Chirotope<3,7> wmi_minus_e = delete_e(wmi);
            if (wmi_minus_e.is_zero()) continue; // e was a coloop!
            // Binary search
            int lower_bound = 0;
            int upper_bound = lc_of_deletion_filtered.size() - 1;
            int midpoint = (lower_bound + upper_bound) / 2;
            while (midpoint != upper_bound) {
                if (compare(lc_of_deletion_filtered[midpoint], wmi_minus_e)) {
                    lower_bound = midpoint + 1;
                } else {
                    upper_bound = midpoint;
                }
                midpoint = (lower_bound + upper_bound) / 2;
            }
            was_hit[midpoint] = true;
        } 
        int not_hit_idx = -1;
        for (int idx = 0; idx < was_hit.size(); ++idx) {
            if (!was_hit[idx]) {
                not_hit_idx = idx;
                break;
            }
        }
        if (not_hit_idx == -1) {
            weakly_reducible_by_element[
                weakly_reducible_by_element.size() - 1
            ] = e;
            std::cout << "[SUCCESS] " << chi << " is weakly reducible by " << e << "\n\n";
            break;
        } else {
            std::cout << "    [:(] the chirotope " << lc_of_deletion_filtered[not_hit_idx]
            << " has no single element extension which is smaller than " << chi << "\n";
        }
    }
    if (weakly_reducible_by_element.back() == -1) {
        std::cout << "[FAILURE] the chirotope " << chi << " is not weakly "
        "reducible by any element\n\n";
    }
    // Iterate
    idx_of_current_target++;
}

bool was_anything_not_reducible = false;
std::cout << "Result compilation:\n";
for (int idx = 0; idx < r3n7_representatives.size(); ++idx) {
    if (weakly_reducible_by_element[idx] != -1) {
        std::cout << r3n7_representatives[idx] << " is weakly reducible by "
        << weakly_reducible_by_element[idx] << "\n";
    } else {
        std::cout << r3n7_representatives[idx] << " is not weakly reducible by anything\n";
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