#pragma once
#include <iostream>
#include "OMtools.hpp"
#include "program_template.hpp"

namespace programs {

// Tests if the lower cone filtering and generating functions
// produce the same (number of) results (per basecount). This
// means that it can only run for parameters for which a database
// of all oriented matroids is available.
template<int R, int N>
int test_lower_cone_generation(const Chirotope<R, N>& top) {
    auto LC_generated = generate_lower_cone(top);
    std::cout << "\n===============\n\n";
    auto LC_filtered = filter_lower_cone(top);
    std::cout << "Comparing results...\n";
    for (auto b = 0; b < binomial_coefficient(N, R); b++) {
        if (LC_generated[b].size() != LC_filtered[b].size()) {
            std::cout << "(;_;) The two algorithms didn't agree on the "
            "number of weak map images with " << b+1 << " bases:\n";
            std::cout << "--- generate_lower_cone() said " << LC_generated[b].size() << "\n";
            std::cout << "--- filter_lower_cone() said   " << LC_filtered[b].size() << "\n";
            return 1;
        }
    }
    std::cout << "(OuO) The two algorithms agreed on the count of weak map images"
    " for any given base count.\n\n";
    return 0;
}

}