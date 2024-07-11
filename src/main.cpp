#include <iostream>
#include <vector>
#include <array>
#include <format>
#include "OMtools.hpp"


#define R0 3
#define N0 9

int program__compute_euler_char_of_JRG_mod_3()
{
if ((R0 != 3) || (N0 != 9))
    throw std::logic_error("Incorrect compilation! This function must be compiled with parameters (3,9)!");
std::vector<std::vector<Chirotope<R0,N0>>> fixed_OMs_by_bases
(binomial_coefficient(N0, R0), std::vector<Chirotope<R0,N0>>());
auto input = ReadChirotopesFromFile<R0,N0>(
    "../../../resources/fixed_oriented_matroid_sets/fixed_OMs_rank3_9elements_group_Z3.txt",
    6
);
auto JRG = OMexamples::JRG();
for (auto p: input) {
    int base_count = p.countbases();
    fixed_OMs_by_bases[base_count - 1].push_back(p);
}
std::vector<Chirotope<R0,N0>> all_fixed_OMs;
std::vector<int> base_counts;
std::vector<std::vector<size_t>> smaller_OM_indices;
std::vector<std::array<size_t, binomial_coefficient(N0,R0)>> lower_cone_face_vectors;
for (int b = 0; b < binomial_coefficient(N0, R0); b++) {
    std::cout << "Scanning OMs with " << b+1 << " bases...\n";
    for (auto c : fixed_OMs_by_bases[b]) {
        int base_count = c.countbases();
        std::vector<size_t> weak_images = smaller_OMs(
            all_fixed_OMs,
            base_counts,
            c,
            base_count
        );
        auto f_vector = face_vector<binomial_coefficient(N0,R0)>(
            lower_cone_face_vectors, weak_images
        );
        // Save OM to list of all OMS:
        all_fixed_OMs.push_back(c);
        base_counts.push_back(base_count);
        smaller_OM_indices.push_back(weak_images);
        lower_cone_face_vectors.push_back(f_vector);
        if (c == JRG) {
            std::cout << "Euler characteristic of JRG: " << euler_characteristic<binomial_coefficient(N0,R0)>(f_vector) 
            << ", mod 3:" <<  euler_characteristic<binomial_coefficient(N0,R0)>(f_vector) % 3 << "\n";
        }
    }
}
auto f = face_vector<binomial_coefficient(N0,R0)>(lower_cone_face_vectors);
std::cout << "Total f-vector: (" << f[0];
for (auto d = 1; d < binomial_coefficient(N0,R0); d++) {
    std::cout << ", " << f[d];
}
std::cout << ")";
return 0;

}

int main() 
{
    return program__compute_euler_char_of_JRG_mod_3();
}