#pragma once

#include <iostream>
#include <vector>
#include "OMtools.hpp"
#include "research_file_template.hpp"
#include "lowercones.hpp"

namespace research {

template<int R, int N>
std::vector<std::vector<Matroid<R, N>>> read_matroids(
    enum verboseness verbose
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

template<int R, int N>
std::vector<std::vector<Chirotope<R,N>>> generate_lower_cone(
    const Chirotope<R, N>& top, 
    enum verboseness verbose
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

template<int R, int N>
std::vector<std::vector<Chirotope<R, N>>> generate_lower_cone(
    const Chirotope<R, N>& top, 
    const std::vector<std::vector<Matroid<R, N>>>& matroids,
    enum verboseness verbose
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

template<int R, int N>
std::vector<std::vector<Chirotope<R, N>>> filter_lower_cone(
    const Chirotope<R, N>& top,
    enum verboseness verbose
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
            auto corrected = top.weak_maps_to(p.second) ? p.second : p.second.inverse();
            all_wmis_by_basecount[p.first - 1].push_back(corrected);
        }
        // INCREMENT
        total++;
        OMs_with_fixed_basecount++;
    }
    if (verbose >= verboseness::checkpoints) {
        std::cout << "Finished parsing OMs with " << last_basecount
        << " bases.\n";
        std::cout << "--- There were " << count_of_wmi_with_fixed_basecount
        << "/" << OMs_with_fixed_basecount << " weak map images for this basecount.\n";
        std::cout << "--- There are " << count_of_wmi << "/"
        << total << " weak map images in total so far.\n";
    } 
    if (verbose >= verboseness::result) {
        std::cout << "Filtering of the lower cone of " << top << " is complete.\n";
    }
    return all_wmis_by_basecount;
}

template<int R, int N>
std::vector<std::vector<Chirotope<R, N>>> lower_cone(
    const Chirotope<R, N>& top,
    enum verboseness verbose
) {
    std::ifstream first_file(database_names::OM_set<R, N>(1,0));
    if (first_file.good()) {
        return filter_lower_cone(top, verbose);
    } else {
        return generate_lower_cone(top, verbose);
    }
}

template<int R, int N>
std::vector<std::vector<Chirotope<R, N>>> lower_cone(
    const Chirotope<R, N>& top,
    const std::vector<std::vector<Matroid<R, N>>>& matroids,
    enum verboseness verbose
) {
    std::ifstream first_file(database_names::OM_set<R,N>(1,0));
    if (first_file.good()) {
        return filter_lower_cone(top, verbose);
    } else {
        return generate_lower_cone(top, matroids, verbose);
    }
}

template<int R, int N>
std::vector<Chirotope<R, N>> loopfree_lower_cone(
    const Chirotope<R, N>& top,
    enum verboseness verbose
) {
    std::vector<Chirotope<R, N>> loopfrees;
    int kept_total = 0;
    int basecount = 1;
    int wmis_total = 0;
    for (auto wmis_with_fixed_basecount: lower_cone(top, verbose)) {
        int kept_with_basecount = 0;
        for (Chirotope<R, N> wmi: wmis_with_fixed_basecount) {
            if (wmi.is_loopfree()) {
                loopfrees.push_back(wmi);
                kept_with_basecount++;
                kept_total++;
            }
        }
        if (verbose >= verboseness::checkpoints) {
            std::cout << kept_with_basecount << "/" << wmis_with_fixed_basecount.size()
            << " weak map images with " << basecount << " bases were loopfree.\n";
        }
        basecount++;
        wmis_total += wmis_with_fixed_basecount.size();
    }
    if (verbose >= verboseness::result) {
        std::cout << kept_total << "/" << wmis_total << " weak map images were loopfree.\n";
    }
    return loopfrees;
}

template<int R, int N>
std::vector<Chirotope<R, N>> loopfree_lower_cone(
    const Chirotope<R, N>& top,
    const std::vector<std::vector<Matroid<R, N>>>& matroids,
    enum verboseness verbose
) {
    std::vector<Chirotope<R, N>> loopfrees;
    int kept_total = 0;
    int basecount = 1;
    int wmis_total = 0;
    for (auto wmis_with_fixed_basecount: lower_cone(top, matroids, verbose)) {
        int kept_with_basecount = 0;
        for (Chirotope<R, N> wmi: wmis_with_fixed_basecount) {
            if (wmi.is_loopfree()) {
                loopfrees.push_back(wmi);
                kept_with_basecount++;
                kept_total++;
            }
        }
        if (verbose >= verboseness::checkpoints) {
            std::cout << kept_with_basecount << "/" << wmis_with_fixed_basecount.size()
            << " weak map images with " << basecount << " bases were loopfree.\n";
        }
        basecount++;
        wmis_total += wmis_with_fixed_basecount.size();
    }
    if (verbose >= verboseness::result) {
        std::cout << kept_total << "/" << wmis_total << " weak map images were loopfree.\n";
    }
    return loopfrees;
}

}