#pragma once
#include <vector>
#include "OMtools.hpp"
#include "research_file_template.hpp"
#include "isolation.hpp"

namespace research {

template<int R, int N>
bool constrains (const Chirotope<R, N>& chi, const sign_vector<N>& cocircuit, char element) {
    std::vector<int> indices_of_zeros = cocircuit.indices_of_zeros();
    for (auto itr = indices_of_zeros.begin(); itr != indices_of_zeros.end(); ++itr) {
        if (*itr == element) {
            indices_of_zeros.erase(itr);
            break;
        }
    }
    if (indices_of_zeros.size() == R - 1) {
        return chi.is_independent(indices_of_zeros);
    } else if (indices_of_zeros.size() < R - 1) {
        return false;
    } else {
        return chi.rank(indices_of_zeros, R - 1)  == R - 1;
    }
}

template<int R, int N, bool loopfree>
bool exists_covariant_nonloop(const Chirotope<R, N>& chi, 
const std::array<const sign_vector<N>*,R>& cocircuits, char element) {
    bit_vector<N> has_no_plus = ~bit_vector<N>();
    bit_vector<N> has_no_minus = ~bit_vector<N>();
    for (auto cocirc_ptr: cocircuits) {
        for (auto i = 0; i < sign_vector<N>::NR_INT32; ++i) {
            has_no_plus[i] &= ~(cocirc_ptr->plus[i]);
            has_no_minus[i] &= ~(cocirc_ptr->minus[i]);
        }
    }
    for (auto covariant_element: has_no_plus.indices_of_ones()) {
        if constexpr (loopfree) {
            if (covariant_element != element && 
            !chi.is_loop(covariant_element)) 
                return true;
        } else {
            if (covariant_element != element)
                return true;
        }
    }
    for (auto covariant_element: has_no_minus.indices_of_ones()) {
        if constexpr (loopfree) {
            if (covariant_element != element && 
            !chi.is_loop(covariant_element)) 
                return true;
        } else {
            if (covariant_element != element)
                return true;
        }
    }
    return false;
}

template<int R, int N, bool loopfree>
bool is_isolated(const Chirotope<R, N>& chi, 
std::vector<sign_vector<N>>& cocircuits, char element) {
    std::vector<int> constraining_cocircuit_indices{};
    std::vector<char> sign_of_element_in_constraining_cocircuit{};
    for (int j = 0; j < cocircuits.size(); ++j) {
        if (constrains(chi, cocircuits[j], element)) {
            constraining_cocircuit_indices.push_back(j);
            sign_of_element_in_constraining_cocircuit.push_back(
                cocircuits[j].get_char(element)
            );
            if (cocircuits[j].minus.get_bit(element)) {
                cocircuits[j].invert();
            }
        }
    }
    for (std::array<char, R> Rtuple: Rtuples::iterator<char,R>(
        constraining_cocircuit_indices.size()
    )) { // Iterate over `R`-element subsets of constraining_cocircuit_indices
        std::array<bool,R> was_already_flipped{};
        std::array<sign_vector<N> const *, R> cocircuit_ptrs;
        for (int t = 0; t < R; ++t) {
            cocircuit_ptrs[t] = &cocircuits[constraining_cocircuit_indices[Rtuple[t]]];
        }
        while (true) { // Iterate over each choice of orienting the given cocircuits
            if (!exists_covariant_nonloop<R,N,loopfree>(chi,cocircuit_ptrs,element)) {
                return true;
            }
            int t = R - 1;
            while(
                t >= 0 && 
                (sign_of_element_in_constraining_cocircuit[Rtuple[t]] != '0'
                || was_already_flipped[t])
            ) {
                was_already_flipped[t] = false;
                --t;
            }
            if (t < 0) break;
            cocircuits[constraining_cocircuit_indices[Rtuple[t]]].invert();
            was_already_flipped[t] = true;
        }
    }
    return false;
}

}