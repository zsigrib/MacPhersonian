#pragma once

#include <iostream>
#include <vector>
#include <random>
#include "OMtools.hpp"
#include "program_template.hpp"

namespace programs {

// This is a test to check if the `Chirotope<3,6>::is_chirotope()`
// function is working properly. It parses the database of all
// oriented matroid, and verifies that `.is_chirotope()` is true for
// all of them. Next, the user is requested to enter chirotopes
// of parameters (3,6), and the program tests if they are in the 
// database (+-) and if `.is_chirotope()` evaluates to true. Lastly,
// once the user entered an invalid chirotope, 100000 tests are
// performed on randomly generated sign-vectors to see that
// being in the database and `.is_chirotope()` agree on all of them.
inline int verify_ischirotope_port() 
{
constexpr int R0 = 3;
constexpr int N0 = 6;
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
                std::cout << Chirotope<R0,N0>::RTUPLES::LIST::array[i][r];
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
                std::cout << Chirotope<R0,N0>::RTUPLES::LIST::array[i][r];
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

}