#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include "OMtools.hpp"
#include "researchlib.hpp"
#include "program_template.hpp"

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
    
    good_wrt_element.push_back(research::weakly_reducible_by(
        chi,
        matroids,
        matroids_to_filter_deletion_with,
        verboseness::info
    ));
    std::cout << "\n";
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