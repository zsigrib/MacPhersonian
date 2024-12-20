#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include "OMtools.hpp"

namespace OMexamples {

// Read all uniform Finschi representatives of a given rank 
// and number of elements into a vector.
template<int R, int N>
std::vector<Chirotope<R, N>> read_uniform_Finschi_representatives() {
    std::ifstream file(database_names::uniform_Finschi<R,N>);
    ignore_n_lines(file, 3);
    std::vector<Chirotope<R,N>> representatives{};
    std::string s;
    while(file) {
        file >> s;
        while(file && s != "=") {
            file >> s;
        }
        if (!file) break;
        file >> s;
        sign_vector<binomial_coefficient(N,R)> encoded(s);
        representatives.push_back(
            OM_operations::decode_Finschi_representative<R,N>(encoded)
        );
    }
    return representatives;
}

// Read all uniform Finschi representatives of a given rank 
// and number of elements into a vector.
template<int R, int N>
std::vector<Chirotope<R, N>> read_all_Finschi_representatives() {
    std::ifstream file(database_names::all_Finschi<R,N>);
    ignore_n_lines(file, 3);
    std::vector<Chirotope<R,N>> representatives{};
    std::string s;
    while(file) {
        file >> s;
        while(file && s != "=") {
            file >> s;
        }
        if (!file) break;
        file >> s;
        sign_vector<binomial_coefficient(N,R)> encoded(s);
        representatives.push_back(
            OM_operations::decode_Finschi_representative<R,N>(encoded)
        );
    }
    return representatives;
}

}