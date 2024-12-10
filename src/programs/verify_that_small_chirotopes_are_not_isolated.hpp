#pragma once
#include "OMtools.hpp"
#include "OMoperations.hpp"
#include "isolation.hpp"
#include "program_utility.hpp"
#include "program_template.hpp"

namespace programs {

// Checks all oriented matroids of a given rank and number 
// of elements, and complains if any of their elements
// are isolated.
template<int R, int N>
int verify_that_small_chirotopes_are_not_isolated() {
//static_assert((R == 3) && (N <= 6), "This function is only available for R = 3 and N <= 6.");
auto input = ReadOMDataFromFiles<Chirotope<R,N>>(
    database_names::OM_set<R,N>,
    6
);
std::vector<std::vector<Chirotope<R,N>>> all_OMs_by_bases(
    binomial_coefficient(N,R),
    std::vector<Chirotope<R,N>>()
);
int count = 0;
int last_basecount = 0;
for (auto p: input) {
    // PARSE
    all_OMs_by_bases[p.first - 1].push_back(p.second);
    for (auto e = 0; e < N; ++e) {
        auto cocircuits = OM_operations::cocircuits(p.second);
        if (!p.second.is_loop(e) &&
        research::is_isolated<R,N>(p.second,cocircuits,e)) {
            std::cout << "The element " << e << " was found "
            "to be isolated in the oriented matroid\n\n";
            utility::print_Rtuples<R,N>();
            std::cout << "\n" << p.second << "\n\nWith cocircuits\n\n";
            for (auto cocircuit: cocircuits)
                std::cout << cocircuit << "\n";
            std::cout << "\nTerminating.";
            return 1;
        }
    }
    // PRINT
    if (last_basecount != p.first) {
        std::cout << "Finished with OMs with " << last_basecount
        << " bases. Running total: " << count << "\n";
        last_basecount = p.first;
    }
    count++;
}
std::cout << "(OuO) Success!! No chirotopes in the database of "
"parameters (" << R << ", " << N << ") have any isolated elements."
"\n\nTerminating.";
return 0;

}

}