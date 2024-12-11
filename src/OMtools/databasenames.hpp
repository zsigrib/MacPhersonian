#pragma once

#include <string>
#include <format>

namespace database_names {

// `OM_set<R, N>(n_bases, idx)` returns the name of the `idx`th file
// which contains rank `R` chirotopes on `N` elements and with
// `n_bases` many bases.
template<int R, int N>
std::string OM_set(int, int);

template<>
std::string OM_set<3,6>(int n_bases, int idx) {
    if (idx != 0) return "";
    return std::format("../../../resources/oriented_matroid_sets/r3n6/OMs_rank3_6elements_{}bases.txt", n_bases);
}

template<>
std::string OM_set<3,7>(int n_bases, int idx) {
    if (n_bases < 27) {
        if (idx != 0) return "";
        return std::format("../../../resources/oriented_matroid_sets/r3n7/OMs_rank3_7elements_{}bases.txt", n_bases);
    }
    return std::format("../../../resources/oriented_matroid_sets/r3n7/OMs_rank3_7elements_{0}bases_part{1}.txt", n_bases, idx + 1);
}

template<int R, int N>
std::string matroid_set(int n_bases, int idx) {
    if (idx != 0) return "";
    return std::format("../../../resources/matroid_sets/r{0}n{1}/{2}_bases_all.txt", R, N, n_bases);
} 

template<int R, int N>
std::string uniform_Finschi = std::format("../../../resources/Finschi_representatives/uniform/r{}n{}.txt", R, N);

template<int R, int N>
std::string all_Finschi = std::format("../../../resources/Finschi_representatives/all/r{}n{}.txt", R, N);

}