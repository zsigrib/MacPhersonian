#pragma once

#include <iostream>
#include <iomanip>
#include "OMs.hpp"
#include "program_template.hpp"

namespace programs {

// Utility functions, mainly for printing debug info.
namespace utility {

// Print `R`-tuples as a table with `R` rows and
// `binomial_coefficient(N, R)` columns, each column
// being an `R`-tuple.
template<int R, int N>
void print_Rtuples() {
    for (int t = 0; t < R; ++t) {
        for (int idx = 0; idx < Chirotope<R,N>::RTUPLES::NR; ++idx) {
            std::cout << int(Chirotope<R,N>::RTUPLES::LIST::array[idx][t]);
        }
        std::cout << "\n";
    }
}

// Print an iterable producing integers as a horitontal
// right-adjusted single row table, with column width of
// `spacing` many characters.
template<typename Iterable>
void print_iterable_of_ints(const Iterable& iter, int spacing=3) {
    for (auto elem: iter) {
        std::cout << std::setw(spacing) << int(elem);
    }
}

// Print an iterable producing integers as a horizontal
// right-adjusted single row table, with column width of
// `spacing` many characters - excluding the comma.
template<typename Iterable>
void print_comma_separated_iterable_of_ints(const Iterable& iter, int spacing=3) {
    bool first = true;
    for (auto elem: iter) {
        if (!first) {
            std::cout << ",";
        }
        first = false;
        std::cout << std::setw(spacing) << int(elem);
    }
}


}

}