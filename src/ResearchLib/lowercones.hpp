#pragma once
#include <iostream>
#include <vector>
#include "OMtools.hpp"
#include "research_file_template.hpp"

namespace research {

// Returns the list of all matroids, grouped by basecount.
// The index is shifted by 1 from the actual basecount.
//
// To do so, a database of all oriented matroids is read.
template<int R, int N>
std::vector<std::vector<Matroid<R, N>>> read_matroids(
    enum verboseness verbose = verboseness::checkpoints
) {

if (verbose >= verboseness::info) {
    std::cout << "Reading the set of all rank " << R
    << " matroids on " << N << " elements...";
}
std::vector<std::vector<Matroid<R, N>>> matroids_by_bases(
    binomial_coefficient(N, R),
    std::vector<Matroid<R, N>>()
);
auto input = ReadOMDataFromFiles<Matroid<R,N>>(
    database_names::matroid_set<R, N>, 0
);
int last_basecount = 0;
for (auto p : input) {
    // PRINT
    if (p.first != last_basecount) {
        if (verbose >= verboseness::checkpoints) {
            std::cout << "Finished loading matroids with "
            << last_basecount << " bases. There were " 
            << matroids_by_bases[last_basecount - 1].size()
            << " of them.\n";
        }
        last_basecount = p.first;
    }
    // PARSE
    matroids_by_bases[p.first - 1].push_back(p.second);
}
if (verbose >= verboseness::result) {
    std::cout << "Finished parsing matroids with "
    << last_basecount << " bases. There were " 
    << matroids_by_bases[last_basecount - 1].size()
    << " of them.\n";
}
if (verbose >= verboseness::info) {
    std::cout << "Finished reading in all matroids.\n";
}
return matroids_by_bases;

}



// Return the set of all weak map images of `top`, grouped by basecount.
// The index is shifted from the basecount by 1.
//
// This uses a database of all matroids.
template<int R, int N>
std::vector<std::vector<Chirotope<R,N>>> generate_lower_cone(
    const Chirotope<R, N>& top, 
    enum verboseness verbose = verboseness::checkpoints
) {

if (!top.is_chirotope())
    throw std::invalid_argument("ERROR: top is not a chirotope.");

std::vector<std::vector<Chirotope<R,N>>> all_wmis_by_basecount(
    binomial_coefficient(N, R), std::vector<Chirotope<R,N>>()
);

auto input = ReadOMDataFromFiles<Matroid<R,N>>(
    database_names::matroid_set<R,N>, 0
);

if (verbose >= verboseness::info) {
    std::cout << "Generating lower cone of " << top << "...\n";
}
int last_basecount = 0;
size_t total = 0;
size_t matroids_with_fixed_basecount = 0;
size_t count_of_wmi = 0;
size_t count_of_wmi_with_fixed_basecount = 0;
int top_basecount = top.countbases();
for (auto p : input) {
    if (top_basecount == p.first) {
        if (verbose >= verboseness::info) {
            std::cout << "We have exhausted all basecounts smaller than"
            " top's basecount (" << top_basecount << ").\n";
        }
        all_wmis_by_basecount[top_basecount - 1].push_back(top);
        break;
    }
    // PRINT
    if (p.first != last_basecount) {
        if (verbose >= verboseness::checkpoints) {
             std::cout << "Finished parsing matroids with "
            << last_basecount << " bases.\n";
            std::cout << "--- There were " << count_of_wmi_with_fixed_basecount 
            << "/" << matroids_with_fixed_basecount << " weak map images for this basecount\n";
            std::cout << "--- There are " << count_of_wmi << "/"
            << total << " weak map images in total so far.\n";
        }
        matroids_with_fixed_basecount = 0;
        count_of_wmi_with_fixed_basecount = 0;
        last_basecount = p.first;
    }
    // PARSE
    if (top.weak_maps_to(p.second)) {
        auto restricted = OM_operations::restrict_to_bases<R, N>(p.second) * top;
        if (restricted.is_chirotope()) {
            count_of_wmi++;
            count_of_wmi_with_fixed_basecount++;
            all_wmis_by_basecount[p.first - 1].push_back(restricted);
        }
    }
    // INCREMENT
    matroids_with_fixed_basecount++;
    total++;
}
if (verbose >= verboseness::result) {
    std::cout << "Finished parsing matroids with "
    << last_basecount << " bases.\n";
    std::cout << "--- There were " << count_of_wmi_with_fixed_basecount 
    << "/" << matroids_with_fixed_basecount << " weak map images for this basecount.\n";
    std::cout << "--- There are " << count_of_wmi << "/"
    << total << " weak map images in total so far.\n";
}
if (verbose >= verboseness::info) {
    std::cout << "Generation of the lower cone of " << top << " is complete.\n";
}
return all_wmis_by_basecount;

}



// Calculate all weak map images of `top` using a precomputed
// set of matroids. Both `matroids` and the output are grouped
// by base count, and indices are shifted by 1 from the actual
// base count.
//
// This uses a database of all matroids.
template<int R, int N>
std::vector<std::vector<Chirotope<R, N>>> generate_lower_cone(
    const Chirotope<R, N>& top, 
    const std::vector<std::vector<Matroid<R, N>>>& matroids,
    enum verboseness verbose = verboseness::checkpoints
) {

if (!top.is_chirotope())
    throw std::invalid_argument("ERROR: top is not a chirotope.");

if (verbose >= verboseness::info) {
    std::cout << "Generating lower cone of " << top << 
    ", given the set of appropriate matroids...\n";
}
std::vector<std::vector<Chirotope<R, N>>> all_wmis_by_bases(
    binomial_coefficient(N,R),
    std::vector<Chirotope<R, N>>()
);
size_t count_of_wmi = 0;
size_t total = 0;
int top_basecount = top.countbases();
for (auto b = 0; b < top_basecount - 1; b++) {
    size_t count_of_wmi_with_fixed_basecount = 0;
    for (auto matroid : matroids[b]) {
        if (top.weak_maps_to(matroid)) {
            auto restricted = OM_operations::restrict_to_bases<R,N>(matroid) * top;
            if (restricted.is_chirotope()) {
                count_of_wmi++;
                count_of_wmi_with_fixed_basecount++;
                all_wmis_by_bases[b].push_back(restricted);
            }
        }
    }
    total += matroids[b].size();
    if (verbose >= verboseness::checkpoints) {
        std::cout << "Finished parsing matroids with " << b+1
        << " bases.\n";
        std::cout << "--- There were " << count_of_wmi_with_fixed_basecount
        << "/" << matroids[b].size() << " weak map images for this basecount.\n";
        std::cout << "--- There are " << count_of_wmi << "/"
        << total << " weak map images in total so far.\n";
    }
}
all_wmis_by_bases[top_basecount - 1].push_back(top);
if (verbose >= verboseness::info) {
    std::cout << "Generation of the lower cone of "
    << top << " is complete; we have checked all relevant basecounts.\n";
}
return all_wmis_by_bases;

}



// Reads all OMs of the given rank and number of elements, and selects
// only those which fall inside the given lower cone.
//
// This uses a database of all oriented matroids.
template<int R, int N>
std::vector<std::vector<Chirotope<R, N>>> filter_lower_cone(
    const Chirotope<R, N>& top,
    enum verboseness verbose = verboseness::checkpoints
) {

if (!top.is_chirotope())
    throw std::invalid_argument("ERROR: top is not a chirotope.");

std::vector<std::vector<Chirotope<R,N>>> all_wmis_by_basecount(
    binomial_coefficient(N, R), std::vector<Chirotope<R,N>>()
);

auto input = ReadOMDataFromFiles<Chirotope<R, N>>(
    database_names::OM_set<R, N>, 6
);

if (verbose >= verboseness::info) {
    std::cout << "Reading all OMs with R = " << R << " and N = " << N
    << ", and keeping those which are in the lower cone of " << top
    << "...\n";
}
int last_basecount = 0;
size_t total = 0;
size_t OMs_with_fixed_basecount = 0;
size_t count_of_wmi = 0;
size_t count_of_wmi_with_fixed_basecount = 0;
int top_basecount = top.countbases();
for (auto p : input) {
    if (top_basecount == p.first) {
        if (verbose >= verboseness::info) {
            std::cout << "Finished with exhausting all basecounts smaller than"
            " top's basecount (" << top_basecount << ").\n";
        }
        all_wmis_by_basecount[top_basecount - 1].push_back(top);
        break;
    }
    // PRINT
    if (last_basecount != p.first) {
        if (verbose >= verboseness::checkpoints) {
            std::cout << "Finished parsing OMs with " << last_basecount
            << " bases.\n";
            std::cout << "--- There were " << count_of_wmi_with_fixed_basecount
            << "/" << OMs_with_fixed_basecount << " weak map images for this basecount.\n";
            std::cout << "--- There are " << count_of_wmi << "/"
            << total << " weak map images in total so far.\n";
        }
        last_basecount = p.first;
        OMs_with_fixed_basecount = 0;
        count_of_wmi_with_fixed_basecount = 0;
    }
    // PARSE
    if (top.OM_weak_maps_to(p.second)) {
        count_of_wmi++;
        count_of_wmi_with_fixed_basecount++;
        all_wmis_by_basecount[p.first - 1].push_back(p.second);
    }
    // INCREMENT
    total++;
    OMs_with_fixed_basecount++;
}
if (verbose >= verboseness::result) {
    std::cout << "Finished parsing OMs with " << last_basecount
    << " bases.\n";
    std::cout << "--- There were " << count_of_wmi_with_fixed_basecount
    << "/" << OMs_with_fixed_basecount << " weak map images for this basecount.\n";
    std::cout << "--- There are " << count_of_wmi << "/"
    << total << " weak map images in total so far.\n";
} 
if (verbose >= verboseness::info) {
    std::cout << "Filtering of the lower cone of " << top << " is complete.\n";
}
return all_wmis_by_basecount;

}

}