#include <iostream>
#include <vector>
#include <array>
#include <format>
#include <random>
#include <algorithm>
#include "OMtools.hpp"


#define R0 3
#define N0 6



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
auto JRG = OMexamples::JRG();
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
        auto ec = euler_characteristic<binomial_coefficient(N0,R0)>(f_vector);
        if (ec % 3 != 1 && ec % 3 != -2) non_contr++;
        // Save OM to list of all OMS:
        all_fixed_OMs.push_back(c);
        base_counts.push_back(base_count);
        smaller_OM_indices.push_back(weak_images);
        lower_cone_face_vectors.push_back(f_vector);
        /*if (c == JRG) {
            std::cout << "Euler characteristic of JRG: " << euler_characteristic<binomial_coefficient(N0,R0)>(f_vector) 
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



template<int R, int N>
std::vector<std::vector<Chirotope<R,N>>> generate_lower_cone(const Chirotope<R, N>& top)
{

std::vector<std::vector<Chirotope<R,N>>> all_wmis_by_basecount(
    binomial_coefficient(N, R), std::vector<Chirotope<R,N>>()
);

auto input = ReadOMDataFromFiles<Matroid<R,N>>(
    database_names::matroid_set<R,N>, 0
);

std::cout << "Generating lower cone of " << top << "...\n\n";
int last_basecount = 0;
int total = 0;
int matroids_with_fixed_basecount = 0;
int count_of_wmi = 0;
int count_of_wmi_with_fixed_basecount = 0;
for (auto p : input) {
    // PRINT
    if (p.first != last_basecount) {
        std::cout << "Finished parsing matroids with "
        << last_basecount << " bases.\n";
        std::cout << "--- There were " << count_of_wmi_with_fixed_basecount 
        << "/" << matroids_with_fixed_basecount << " weak map images for this basecount\n";
        std::cout << "--- There are " << count_of_wmi << "/"
        << total << " weak map images in total so far.\n";
        matroids_with_fixed_basecount = 0;
        count_of_wmi_with_fixed_basecount = 0;
        last_basecount = p.first;
    }
    // PARSE
    auto restriced = top.restrict_to_matroid(p.second);
    if (top.weak_maps_to(p.second) && restriced.is_chirotope()) {
        count_of_wmi++;
        count_of_wmi_with_fixed_basecount++;
        all_wmis_by_basecount[p.first - 1].push_back(restriced);
    }
    // INCREMENT
    matroids_with_fixed_basecount++;
    total++;
}
std::cout << "Finished parsing matroids with "
<< last_basecount << " bases.\n";
std::cout << "--- There were " << count_of_wmi_with_fixed_basecount 
<< "/" << matroids_with_fixed_basecount << " weak map images for this basecount.\n";
std::cout << "--- There are " << count_of_wmi << "/"
<< total << " weak map images in total so far.\n";
std::cout << "Generation of the lower cone of " << top << " is complete.\n";
return all_wmis_by_basecount;

}



// Reads all OMs of the given rank and number of elements, and selects
// only those which fall inside the given lower cone.
template<int R, int N>
std::vector<std::vector<Chirotope<R, N>>> filter_lower_cone(const Chirotope<R, N>& top) 
{

std::vector<std::vector<Chirotope<R,N>>> all_wmis_by_basecount(
    binomial_coefficient(N, R), std::vector<Chirotope<R,N>>()
);

auto input = ReadOMDataFromFiles<Chirotope<R, N>>(
    database_names::OM_set<R, N>, 6
);

std::cout << "Reading all OMs with R = " << R << " and N = " << N
<< ", and keeping those which are in the lower cone of " << top
<< "...\n\n";
int last_basecount = 0;
int total = 0;
int OMs_with_fixed_basecount = 0;
int count_of_wmi = 0;
int count_of_wmi_with_fixed_basecount = 0;
for (auto p : input) {
    // PRINT
    if (last_basecount != p.first) {
        std::cout << "Finished parsing OMs with " << last_basecount
        << " bases.\n";
        std::cout << "--- There were " << count_of_wmi_with_fixed_basecount
        << "/" << OMs_with_fixed_basecount << " weak map images for this basecount.\n";
        std::cout << "--- There are " << count_of_wmi << "/"
        << total << " weak map images in total so far.\n";
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
std::cout << "Finished parsing OMs with " << last_basecount
<< " bases.\n";
std::cout << "--- There were " << count_of_wmi_with_fixed_basecount
<< "/" << OMs_with_fixed_basecount << " weak map images for this basecount.\n";
std::cout << "--- There are " << count_of_wmi << "/"
<< total << " weak map images in total so far.\n";
std::cout << "Filtering of the lower cone of " << top << " is complete.\n";
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

int main() 
{
    Chirotope<3,6> chi("++--+---+-+-------++");
    return program__test_lower_cone_generation(chi);
}