#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include "OMtools.hpp"
#include "researchlib.hpp"
#include "program_template.hpp"
#include "program_utility.hpp"

namespace programs {

// Helper class.
// Takes in batches of lists of euler characteristics
// and prints to std::cout gradual reports on their
// distributions.
class EulerCharAnalyzer {
private:
    size_t ecs_processed = 0;
    size_t count_of_ec_0 = 0;
    size_t count_of_ec_1 = 0;
    size_t count_of_ec_2 = 0;
    std::vector<long long> ec_batch;
    size_t ec_0_in_this_batch = 0;
    size_t ec_1_in_this_batch = 0;
    size_t ec_2_in_this_batch = 0;
public:
    void add_entry(long long ec) {
        ec_batch.push_back(ec);
        switch (ec) {
            case 0:
                ++count_of_ec_0;
                ++ec_0_in_this_batch;
                break;
            case 1:
                ++count_of_ec_1;
                ++ec_1_in_this_batch;
                break;
            case 2:
                ++count_of_ec_2;
                ++ec_2_in_this_batch;
                break;
        }
    }
    void end_batch() {
        if (ec_batch.size() == 0) {
            std::cout << "This batch of OMs was empty.\n"; return;
        }
        std::cout << "The lower cones of this batch of " << ec_batch.size()
        << " OMs had the following distribution of Euler characteristics:\n";
        std::sort (ec_batch.begin(), ec_batch.end());
        size_t idx_of_start_of_this_euler_char = 0;
        for (size_t i = 1; i < ec_batch.size(); ++i) {
            if (ec_batch[i] == ec_batch[i-1]) continue;
            report_on_1_ec(idx_of_start_of_this_euler_char, i-1);
            idx_of_start_of_this_euler_char = i;
        }
        report_on_1_ec(idx_of_start_of_this_euler_char, ec_batch.size() - 1);
        ecs_processed += ec_batch.size();
        std::cout << "The number of lower cones of given Euler characteristic in this batch and in total so far are:\n";
        utility::ArrayFormatter array_formatter;
        array_formatter.add_cell(std::string("Euler characteristic:"), utility::ArrayFormatter::LEFT);
        array_formatter.add_cell(0);
        array_formatter.add_cell(2);
        array_formatter.add_cell(std::string("0 or 2"));
        array_formatter.add_cell(std::string("any"));
        array_formatter.add_cell(1);
        array_formatter.add_cell(std::string("other"));
        array_formatter.newline();
        array_formatter.add_cell(std::string("Count in this batch:"), utility::ArrayFormatter::LEFT);
        array_formatter.add_cell(ec_0_in_this_batch);
        array_formatter.add_cell(ec_2_in_this_batch);
        array_formatter.add_cell(ec_0_in_this_batch + ec_2_in_this_batch);
        array_formatter.add_cell(ec_batch.size());
        array_formatter.add_cell(ec_1_in_this_batch);
        array_formatter.add_cell(ec_batch.size() - ec_0_in_this_batch - ec_1_in_this_batch - ec_2_in_this_batch);
        array_formatter.newline();
        array_formatter.add_cell(std::string("Running count:"), utility::ArrayFormatter::LEFT);
        array_formatter.add_cell(count_of_ec_0);
        array_formatter.add_cell(count_of_ec_2);
        array_formatter.add_cell(count_of_ec_0 + count_of_ec_2);
        array_formatter.add_cell(ecs_processed);
        array_formatter.add_cell(count_of_ec_1);
        array_formatter.add_cell(ecs_processed - count_of_ec_0 - count_of_ec_1 - count_of_ec_2);
        array_formatter.print(false, " | ", "", "    ");
        ec_batch.clear();
        ec_0_in_this_batch = 0;
        ec_1_in_this_batch = 0;
        ec_2_in_this_batch = 0;
    }
private:
    // Helper
    void report_on_1_ec(size_t index_of_first, size_t index_of_last) {
        switch (ec_batch[index_of_first]) {
            case 0:
                std::cout << "   *0*   "; break;
            case 1:
                std::cout << "   <1>   "; break;
            case 2:
                std::cout << "   *2*   "; break;
            default:
                std::cout << "         ";
        }
        std::cout << ec_batch[index_of_first] << " appeared with multiplicity "
        << index_of_last - index_of_first + 1 << ",\n";
    }
};



// Using a precomputed database of all oriented matroids
// of a given rank and number of elements, compute the euler
// characteristic of the strict lower cone of each one, and
// display the results grouped by basecount and euler
// characteristic.
template<int R, int N>
int euler_chars_of_all_lowercones_by_bases_using_database() 
{
static_assert((R == 3 && N == 6) || (R == 3 && N == 7),
"This program must be compiled with parameters (3,6) or (3,7)!");
auto input = ReadOMDataFromFiles<Chirotope<R,N>>(
    &database_names::OM_set<R,N>,
    6
);

std::vector<Chirotope<R,N>> all_OMs;
std::vector<int> base_counts;
std::vector<std::vector<size_t>> smaller_OM_indices;
std::vector<std::array<size_t, binomial_coefficient(N,R)>> lower_cone_face_vectors;
size_t id = 0;
EulerCharAnalyzer ec_analyzer;
int current_basecount = 1;
for (auto p: input) {
    // Parse new OM:
    auto weak_images = smaller_OMs(
        all_OMs,
        base_counts,
        p.second,
        p.first
    );
    auto fvector = face_vector<binomial_coefficient(N,R)>(
        lower_cone_face_vectors,
        weak_images
    );
    // Print message:
    if (p.first > current_basecount) {
        std::cout << "[" << current_basecount << "] Finished parsing OMs with " 
        << current_basecount << " bases. ";
        ec_analyzer.end_batch();
        current_basecount = p.first;
    }
    // Continue ec counting:
    auto ec = euler_characteristic<binomial_coefficient(N,R)>(fvector);
    ec_analyzer.add_entry(ec);
    // Save results:
    all_OMs.push_back(p.second);
    base_counts.push_back(p.first);
    smaller_OM_indices.push_back(weak_images);
    lower_cone_face_vectors.push_back(fvector);
    // Increment
    id++;
}
std::cout << "[" << current_basecount << "] Finished parsing OMs with " 
<< current_basecount << " bases. ";
ec_analyzer.end_batch();
return 0;
}





}