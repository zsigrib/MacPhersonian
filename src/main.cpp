#include <iostream>
#include <vector>
#include <array>
#include <format>
#include <random>
#include <algorithm>
#include "OMtools.hpp"


#define R0 3
#define N0 9


enum verboseness {silent, result, info, checkpoints};

// PROGRAM

int program__compute_fvector_of_MacP()
{
if (!((R0 == 3 && N0 == 6) || (R0 == 3 && N0 == 7)))
    throw std::logic_error("Incorrect compilation! This function must be compiled with parameters (3,6) or (3,7)!");
auto input = ReadOMDataFromFiles<Chirotope<R0,N0>>(
    &database_names::OM_set<R0,N0>,
    6
);

std::vector<Chirotope<R0,N0>> all_OMs;
std::vector<int> base_counts;
std::vector<std::vector<size_t>> smaller_OM_indices;
std::vector<std::array<size_t, binomial_coefficient(N0,R0)>> lower_cone_face_vectors;
size_t id = 0;
int current_basecount = 1;
size_t OMs_with_fixed_basecount = 0;
size_t OMs_with_fixed_basecount_and_good_ec = 0;
size_t total_good_ec = 0;
for (auto p: input) {
    // Parse new OM:
    auto weak_images = smaller_OMs(
        all_OMs,
        base_counts,
        p.second,
        p.first
    );
    auto fvector = face_vector<binomial_coefficient(N0,R0)>(
        lower_cone_face_vectors,
        weak_images
    );
    // Print message:
    if (p.first > current_basecount) {
        std::cout << "Finished parsing OMs with " << current_basecount << " bases. There were "
        << OMs_with_fixed_basecount_and_good_ec << "/" << OMs_with_fixed_basecount 
        << " many of them with non-1 Euler-characteristic.\n";
        total_good_ec += OMs_with_fixed_basecount_and_good_ec;
        OMs_with_fixed_basecount = 0;
        OMs_with_fixed_basecount_and_good_ec = 0;
        current_basecount = p.first;
    }
    // Continue ec counting:
    auto ec = euler_characteristic<binomial_coefficient(N0,R0)>(fvector);
    if (ec != 1) OMs_with_fixed_basecount_and_good_ec++;
    // Save results:
    all_OMs.push_back(p.second);
    base_counts.push_back(p.first);
    smaller_OM_indices.push_back(weak_images);
    lower_cone_face_vectors.push_back(fvector);
    // Increment
    OMs_with_fixed_basecount++;
    id++;
}
std::cout << "Finished parsing OMs with " << current_basecount 
<< " bases. There were " << OMs_with_fixed_basecount_and_good_ec << "/" 
<< OMs_with_fixed_basecount << " many of them with non-1 Euler-characteristic.\n";
total_good_ec += OMs_with_fixed_basecount_and_good_ec;

std::cout << "There were " << total_good_ec << "/" << all_OMs.size()
<< " many non-1 Euler characteristic lower cones.\n";

auto f = face_vector<binomial_coefficient(N0,R0)>(lower_cone_face_vectors);
std::cout << "Total f-vector: (" << f[0];
for (auto d = 1; d < binomial_coefficient(N0,R0); d++) {
    std::cout << ", " << f[d];
}
std::cout << ")";
return 0;
}



// PROGRAM

int program__compute_euler_char_of_JRG_mod_3()
{
if ((R0 != 3) || (N0 != 9))
    throw std::logic_error("Incorrect compilation! This function must be compiled with parameters (3,9)!");
std::vector<std::vector<Chirotope<R0,N0>>> fixed_OMs_by_bases
(binomial_coefficient(N0, R0), std::vector<Chirotope<R0,N0>>());
auto input = ReadOMDataFromFile<Chirotope<R0,N0>>(
    "../../../resources/fixed_oriented_matroid_sets/fixed_OMs_rank3_9elements_group_Z3.txt",
    6
);
int count = 0;
for (auto p: input) {
    int base_count = p.countbases();
    fixed_OMs_by_bases[base_count - 1].push_back(p);
    count++;
}
std::vector<Chirotope<R0,N0>> all_fixed_OMs;
std::vector<int> base_counts;
std::vector<std::vector<size_t>> smaller_OM_indices;
std::vector<std::array<size_t, binomial_coefficient(N0,R0)>> lower_cone_face_vectors;
for (int b = 0; b < binomial_coefficient(N0, R0); b++) {
    std::cout << "Scanning OMs with " << b+1 << " bases...\n";
    size_t non_contr = 0;
    for (auto c : fixed_OMs_by_bases[b]) {
        std::vector<size_t> weak_images = smaller_OMs(
            all_fixed_OMs,
            base_counts,
            c,
            b+1
        );
        auto f_vector = face_vector<binomial_coefficient(N0,R0)>(
            lower_cone_face_vectors, weak_images
        );
        auto ec = euler_characteristic<binomial_coefficient(N0,R0)>(f_vector);
        if (ec % 3 != 1 && ec % 3 != -2) non_contr++;
        // Save OM to list of all OMS:
        all_fixed_OMs.push_back(c);
        base_counts.push_back(b+1);
        smaller_OM_indices.push_back(weak_images);
        lower_cone_face_vectors.push_back(f_vector);
        /*if (c == OMexamples::RIN9) {
            std::cout << "Euler characteristic of RIN9: " << euler_characteristic<binomial_coefficient(N0,R0)>(f_vector) 
            << ", mod 3:" <<  euler_characteristic<binomial_coefficient(N0,R0)>(f_vector) % 3 << "\n";
        }*/
    }
    std::cout << "# of non 1 (mod 3) Euler-characteristic lower cones: "
    << non_contr << "/" << fixed_OMs_by_bases[b].size() << " (# of bases: "
    << b+1 <<").\n";
}
auto f = face_vector<binomial_coefficient(N0,R0)>(lower_cone_face_vectors);
std::cout << "Total f-vector: (" << f[0];
for (auto d = 1; d < binomial_coefficient(N0,R0); d++) {
    std::cout << ", " << f[d];
}
std::cout << ")";
return 0;

}



// PROGRAM

int program__verify_ischirotope_port() 
{

if ((R0 != 3) || (N0 != 6))
    throw std::logic_error("Incorrect compilation! This function must be compiled with parameters (3,6)!");
auto input = ReadOMDataFromFiles<Chirotope<R0,N0>>(
    database_names::OM_set<R0,N0>,
    6
);
std::vector<std::vector<Chirotope<R0,N0>>> all_OMs_by_bases(
    binomial_coefficient(N0,R0),
    std::vector<Chirotope<R0,N0>>()
);
int count = 0;
int last_basecount = 0;
for (auto p: input) {
    // PARSE
    all_OMs_by_bases[p.first - 1].push_back(p.second);
    if (!p.second.is_chirotope()) {
        std::cout << "(;_;) Found a chirotope in the database for which"
        " is_chirotope() evaluates to false: \n" << p.second << "\n";
        for (int r = 0; r < R0; r++) {
            for (int i = 0; i < binomial_coefficient(N0, R0); i++) {
                std::cout << Chirotope<R0,N0>::RTUPLES_LIST::array[i][r];
            }
            std::cout << "\n";
        }
        return 1;
    }
    // PRINT
    if (last_basecount != p.first) {
        std::cout << "Finished with OMs with " << last_basecount
        << " bases. Running total: " << count << "\n";
        last_basecount = p.first;
    }
    count++;
}
std::cout << "(OuO) Success!! All chirotopes in the database are "
"recognised as chirotopes by the program.\n";
while (true) {
    Chirotope<R0,N0> chi;
    std::cout << "                   ";
    for (int x = 0; x < binomial_coefficient(N0,R0); x++) {
        std::cout << "_";
    }
    std::cout << "\n";
    std::cout << "Enter a chirotope: ";
    try {
        std::cin >> chi;
    }
    catch(const std::exception& e)
    {
        std::cout << "That was not a valid chirotope.\n";
        break;
    }
    bool is_in_database = false;
    for (auto c : all_OMs_by_bases[chi.countbases() - 1]) {
        if (c.is_same_OM_as(chi)) {
            is_in_database = true;
            break;
        }
    }
    std::cout << "--- is_chirotope() == " << chi.is_chirotope() << "\n";
    std::cout << "--- is_in_database == " << is_in_database << "\n";
}

size_t TESTS = 100000;
std::cout << "Now testing " << TESTS << " random functions if they are chirotopes...\n";
size_t was_chirotope = 0;
std::random_device dev;
std::mt19937 rng(dev());
const char CHARS[3] {'0','-','+'};
std::uniform_int_distribution<std::mt19937::result_type> dist3(0,2);
for (auto i = 0; i < TESTS; i++) {
    Chirotope<R0,N0> c;
    for (auto idx = 0; idx < binomial_coefficient(N0,R0); idx++) {
        c.set(idx, CHARS[dist3(rng)]);
    }
    bool is_in_database = false;
    for (auto c2 : all_OMs_by_bases[c.countbases() - 1]) {
        if (c2.is_same_OM_as(c)) {
            is_in_database = true;
            break;
        }
    }
    if (c.is_chirotope() != is_in_database) {
        std::cout << "(;_;) Found a chirotope in the database for which"
        " is_chirotope() != is_in_database: \n" << c << "\n";
        for (int r = 0; r < R0; r++) {
            for (int i = 0; i < binomial_coefficient(N0, R0); i++) {
                std::cout << Chirotope<R0,N0>::RTUPLES_LIST::array[i][r];
            }
            std::cout << "\n";
        }
        std::cout << "--- is_chirotope() == " << c.is_chirotope() << "\n";
        std::cout << "--- is_in_database == " << is_in_database << "\n";
        return 1;
    }
    if (is_in_database) was_chirotope++;
}
std::cout << "All tests passed successfully. There were "
<< was_chirotope << "/" << TESTS << " functions which were chirotopes.\n";
return 0;

}


// Returns the list of all matroids, grouped by basecount.
// The index is shifted by 1 from the actual basecount.
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
        auto restricted = top.restriction_to_matroid(p.second);
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
            auto restricted = top.restriction_to_matroid(matroid);
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
            std::cout << "We have exhausted all basecounts smaller than"
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



// PROGRAM

template<int R, int N>
int program__test_lower_cone_generation(const Chirotope<R, N>& top) {
    auto LC_generated = generate_lower_cone(top);
    std::cout << "\n===============\n\n";
    auto LC_filtered = filter_lower_cone(top);
    std::cout << "Comparing results...\n";
    for (auto b = 0; b < binomial_coefficient(N, R); b++) {
        if (LC_generated[b].size() != LC_filtered[b].size()) {
            std::cout << "(;_;) The two algorithms didn't agree on the "
            "number of weak map images with " << b+1 << " bases:\n";
            std::cout << "--- generate_lower_cone() said " << LC_generated[b].size() << "\n";
            std::cout << "--- filter_lower_cone() said   " << LC_filtered[b].size() << "\n";
            return 1;
        }
    }
    std::cout << "(OuO) The two algorithms agreed on the count of weak map images"
    " for any given base count.\n\n";
    return 0;
}



// PROGRAM

template<int R, int N>
int program__compute_fvector_of_lowercone(
    const Chirotope<R, N>& top,
    const std::vector<std::vector<Chirotope<R, N>>>& targets,
    const std::vector<std::vector<Chirotope<R,N>>>& wmis_by_bases_={},
    bool wmis_precalculated = true,
    enum verboseness verbose = verboseness::checkpoints
) {

for (auto target_list: targets) {
    for (auto target : target_list) {
        if (!top.OM_weak_maps_to(target)) {
            std::cerr << "The list of targets included        " << target
            << ",\nbut this is not a weak map image of " << top;
            return 1;
        }
    }
}

if (!wmis_precalculated && verbose >= verboseness::info) {
    std::cout << "The set of weak map images is not precalculated. "
    "Calculating it...\n\n";
}
const auto& wmis_by_bases(wmis_precalculated? wmis_by_bases_ : generate_lower_cone(top));

std::vector<Chirotope<R, N>> all_wmis;
std::vector<int> base_counts;
//std::vector<std::vector<size_t>> smaller_OM_indices;
std::vector<std::array<size_t, binomial_coefficient(N0,R0)>> lower_cone_face_vectors;
int top_basecount = top.countbases();
for (int b = 0; b < top_basecount; b++) {
    if (verbose >= verboseness::checkpoints) {
        std::cout << "Computing face vectors of OMs with " << b+1 << " bases...\n";
        std::cout << ">.< " << verbose << "\n";
    }
    size_t non_contr = 0;
    for (auto c : wmis_by_bases[b]) {
        // PARSE
        std::vector<size_t> weak_images = smaller_OMs(
            all_wmis,
            base_counts,
            c,
            b+1
        );
        auto f_vector = face_vector<binomial_coefficient(N0,R0)>(
            lower_cone_face_vectors, weak_images
        );
        auto ec = euler_characteristic<binomial_coefficient(N0,R0)>(f_vector);
        if (ec != 1) non_contr++;
        // PRINT
        for (auto target : targets[b]) {
            if (target.is_same_OM_as(c) && verbose >= verboseness::result) {
                std::cout << "Found target OM " << target << "! "
                "Euler-characteristic: " << ec << ", face vector: ("
                << f_vector[0];
                for (auto d = 1; d < binomial_coefficient(N0,R0); d++) {
                    std::cout << ", " << f_vector[d];
                }
                std::cout << ")\n";
            }
        }
        // SAVE
        all_wmis.push_back(c);
        base_counts.push_back(b+1);
        //smaller_OM_indices.push_back(weak_images);
        lower_cone_face_vectors.push_back(f_vector);
    }
    if (verbose >= verboseness::checkpoints) {
        std::cout << "# of non 1 Euler-characteristic lower cones: "
        << non_contr << "/" << wmis_by_bases[b].size() << " (# of bases: "
        << b+1 << ").\n";
    }
}
if (verbose >= verboseness::info) {
    std::cout << "All face vectors computed successfully. Terminating.\n";
}
return 0;

}



// PROGRAM

template<int R, int N>
int program__compute_f_vector_of_single_element(
    const Chirotope<R,N>& chi,
    enum verboseness verbose = verboseness::checkpoints
) {

std::vector<std::vector<Chirotope<R,N>>> targets(
    binomial_coefficient(N,R), 
    std::vector<Chirotope<R,N>>()
);
targets[chi.countbases() - 1].push_back(chi);
return program__compute_fvector_of_lowercone<R,N>(
    chi,
    targets,
    std::vector<std::vector<Chirotope<R,N>>>{},
    false,
    verbose
);

}



// PROGRAM

template<int R, int N>
int program__compute_f_vector_of_independent_elements(
    const std::vector<Chirotope<R, N>>& targets,
    enum verboseness verbose = verboseness::checkpoints
) {

auto matroids = read_matroids<R, N>(verbose);

int target_idx = 0;
for (auto chi: targets) {
    std::cout << "\n[TARGET " << target_idx+1 << "/" 
    << targets.size() << "]\n";
    std::vector<std::vector<Chirotope<R,N>>> targets_(
        binomial_coefficient(N,R), 
        std::vector<Chirotope<R,N>>()
    );
    targets_[chi.countbases() - 1].push_back(chi);
    auto lc = generate_lower_cone(
        chi, 
        matroids, 
        std::min(verboseness::result, verbose)
    );
    std::cout << "\n";
    auto ret = program__compute_fvector_of_lowercone(
        chi,
        targets_,
        lc,
        true,
        std::min(verboseness::result, verbose)
    );
    if (ret) return ret;
    target_idx++;
}
std::cout << "\nFinished with all targets.";
return 0;

}

int main() 
{
    std::vector<Chirotope<4,8>> targets;
    for (auto c : OMexamples::NON_REALIZABLES) {
        targets.push_back(c);
    }
    return program__compute_f_vector_of_independent_elements(
        targets
    );
    //return program__compute_f_vector_of_single_element(Chirotope<2,5>("+++0000000"));
}