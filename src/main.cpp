#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include "OMtools.hpp"
#include "researchlib.hpp"
#include "programs.hpp"

int main() 
{
    std::vector<Chirotope<4,8>> targets;
    for (auto c : OMexamples::NON_REALIZABLES) {
        targets.push_back(c);
    }
    return programs::compute_f_vector_of_independent_elements(
        targets
    );
}