#pragma once
#include <float.h>
#include <stdexcept>
#include <vector>
#include <list>
#include "OMtools.hpp"
#include "ordercomplexes.hpp"
#include "research_file_template.hpp"
#include "signvectors.hpp"

namespace research {

// Helper class for collapsibility tests of order
// complexes. Keeps track of the recursive calls
// of the collapsibility tests, and based on this
// info enables optimizations and determines how 
// much effort should be put into a single 
// collapsibility test.
struct _gcProtocol {
    enum LastCone {INIT, NONE, LOWER, UPPER};
    LastCone last_cone = INIT;
    void start_collapsibiliy_test() {
        switch (last_cone) {
            case INIT:
                last_cone = NONE;
                break;
            case NONE:
                last_cone = LOWER;
                break;
            case LOWER:
                last_cone = UPPER;
                break;
            case UPPER:
                last_cone = LOWER;
                break;
        }
    }
    void end_collapsibility_test(bool with_result) {
        switch (last_cone) {
            case INIT:
                throw std::runtime_error("Ending a collapsibility test"
                " which hasn't started yet!");
                break;
            case NONE:
                break; // If a 'NONE' type test finishes, then the program is ending
            case LOWER:
                last_cone = UPPER;
                break;
            case UPPER:
                last_cone = LOWER;
                break;
        }
    }
};

// Helper function for is_greedy_collapsible with the same functionality.
template<typename Comparable, bool (*less_than)(const Comparable&, const Comparable&)>
bool _is_greedy_collapsible(std::list<const Comparable*>& objects, _gcProtocol protocol=_gcProtocol()) {
    protocol.start_collapsibiliy_test();
    if (objects.empty()) {
        protocol.end_collapsibility_test(false);
        return false;
    }
    switch (protocol.last_cone) {
        case _gcProtocol::INIT:
            throw std::runtime_error("_gcProtocol.last_cone == INIT, but"
            " the protocol should have been initialized by now!");
            break;
        case _gcProtocol::NONE:
            while (++objects.cbegin() != objects.cend()) {
                // Look for element for which either the upper or lower cone is
                // collapsible.
                bool found_object_with_contractible_link = false;
                for (auto iter = objects.cbegin(); iter != objects.cend(); ++iter) {
                    // Lower cone
                    std::list<const Comparable*> lower_cone; 
                    for (auto smaller_iter = objects.begin(); smaller_iter != iter; ++smaller_iter) {
                        if (less_than(**smaller_iter, **iter)) {
                            lower_cone.push_back(*smaller_iter);
                        }
                    }
                    if (_is_greedy_collapsible<Comparable, less_than>(lower_cone)) {
                        objects.erase(iter);
                        found_object_with_contractible_link = true;
                        break;
                    }
                    // Upper cone
                    std::list<const Comparable*> upper_cone; 
                    for (auto bigger_iter = --objects.cend(); bigger_iter != iter; --bigger_iter) {
                        if (less_than(**iter, **bigger_iter)) {
                            upper_cone.push_front(*bigger_iter);
                        }
                    }
                    if (_is_greedy_collapsible<Comparable, less_than>(upper_cone)) {
                        objects.erase(iter);
                        found_object_with_contractible_link = true;
                        break;
                    }
                }
                if (!found_object_with_contractible_link) {
                    protocol.end_collapsibility_test(false);
                    return false;
                }
            }
            break;
        case research::_gcProtocol::LOWER:
            while (++objects.cbegin() != objects.cend()) {
                // Look for element for which either the upper or lower cone is
                // collapsible.
                bool found_object_with_contractible_link = false;
                for (auto iter = --objects.cend();true; --iter) {
                    // Upper cone
                    std::list<const Comparable*> upper_cone; 
                    for (auto bigger_iter = --objects.cend(); bigger_iter != iter; --bigger_iter) {
                        if (less_than(**iter, **bigger_iter)) {
                            upper_cone.push_front(*bigger_iter);
                        }
                    }
                    if (_is_greedy_collapsible<Comparable, less_than>(upper_cone)) {
                        objects.erase(iter);
                        found_object_with_contractible_link = true;
                        break;
                    }
                    if (iter == objects.cbegin()) {
                        break;
                    }
                }
                if (!found_object_with_contractible_link) {
                    protocol.end_collapsibility_test(false);
                    return false;
                }
            }
            break;
        case _gcProtocol::UPPER:
            while (++objects.cbegin() != objects.cend()) {
                // Look for element for which either the upper or lower cone is
                // collapsible.
                bool found_object_with_contractible_link = false;
                for (auto iter = objects.cbegin(); iter != objects.cend(); ++iter) {
                    // Lower cone
                    std::list<const Comparable*> lower_cone; 
                    for (auto smaller_iter = objects.begin(); smaller_iter != iter; ++smaller_iter) {
                        if (less_than(**smaller_iter, **iter)) {
                            lower_cone.push_back(*smaller_iter);
                        }
                    }
                    if (_is_greedy_collapsible<Comparable, less_than>(lower_cone)) {
                        objects.erase(iter);
                        found_object_with_contractible_link = true;
                        break;
                    }
                }
                if (!found_object_with_contractible_link) {
                    protocol.end_collapsibility_test(false);
                    return false;
                }
            }
            break;
    }
    
    protocol.end_collapsibility_test(true);
    return true;
}

// If this returns true, then the order complex of the given set
// of oriented matroids is collapsible, and therefore contractible.
//
// Input: a sequence of distinct objects topologically ordered to 
// be non-decreasing.
//
// Algorithm: find the first object for which either
// the lower cone or the upper cone tests as collapsible using this
// function, and remove it, and call this function recursively on the
// remaining list.
template<typename Comparable, bool (*less_than)(const Comparable&, const Comparable&)>
bool is_greedy_collapsible(const std::vector<Comparable>& objects) {
    std::list<const Comparable*> objects_list;
    for (const Comparable& object: objects) {
        objects_list.push_back(&object);
    }
    return _is_greedy_collapsible<Comparable, less_than>(objects_list);
}

template<int R, int N>
bool is_greedy_collapsible(const std::vector<Chirotope<R, N>>& objects) {
    return is_greedy_collapsible<Chirotope<R, N>, &is_OM_weak_map_of>(objects);
}

template<int L>
bool is_greedy_collapsible(const std::vector<bit_vector<L>>& objects) {
    return is_greedy_collapsible<bit_vector<L>, &less_than>(objects);
}

}