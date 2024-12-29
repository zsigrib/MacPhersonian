#pragma once

#include <vector>
#include "research_file_template.hpp"

namespace research {

// Iterates over the `R`-element subsets of `0..N-1`
// (also called `R`-tuples) as monotonely increasing
// `R`-long sequences in `0..N-1` in lexicographic order, 
// and returns them as sorted arrays, assuming all elements
// of the tuple fall inside the subset given by `filter`
//
// The type `LabelType` should be able to store non-negative
// integers up to `N`.
template<typename LabelType=char>
struct tuple_iterator_on_subset {
	// The current `R`-tuple that the iterator is pointing to.
	std::vector<LabelType> current;
	// The number of elements in a tuple.
    const LabelType R;
    // The number of labels, i.e. an `R`-tuple can contain elements
	// in `0..N-1`.
	LabelType N;
    // The subset of `0..N-1` whose `R`-element subsets we list.
    const std::vector<bool> filter;
	// Construct an iterator, pointing to the first `R`-tuple
	// `(0,1,2,...,R-1)`, inside the set of `R`-tuples in `0..n-1`.
	constexpr tuple_iterator_on_subset(LabelType r, LabelType n, const std::vector<bool>& f): 
    R(r), N(n), filter(f) {
		current = std::vector<LabelType>(R,0);
        LabelType s = 0;
        LabelType t = 0;
        while (t < R) {
            while (s < N && !filter[s]) ++s;
            if (s == N) {
                current[0] = N;
                return;
            }
            current[t] = s; 
            ++s;
            ++t;
        }
		if (R == 0) N = 0;
	}
	// Construct an iterator, pointing to an invalid `R`-tuple,
	// consisting solely of `n`s.
	constexpr tuple_iterator_on_subset(LabelType r, LabelType n, const std::vector<bool>& f, int): R(r), N(n), filter(f) {
		current = std::vector<LabelType>(R,0);
		if (R > 0) current[0] = N;
		else N = 1;
	}
	// Construct an iterator, pointing to the specified `R`-tuple.
	constexpr tuple_iterator_on_subset(LabelType r, LabelType n, const std::vector<bool>& f, const std::vector<LabelType>& target): 
    R(r), N(n), filter(f) {
		current = target;
	}
	// Construct an iterator, pointing to the first `R`-tuple
	// `(0,1,2,...,R-1)`, inside the set of `R`-tuples in `0..n-1`.
	constexpr tuple_iterator_on_subset begin() const {
		return *this;
	}
	// Construct an iterator, pointing to an invalid `R`-tuple, which
	// signifies the end of the iteration.
	constexpr tuple_iterator_on_subset end() const {
		return tuple_iterator_on_subset(R,N,filter,0);
	}
	// Return the `R`-tuple this iterator is currently pointing to.
	constexpr std::vector<LabelType> operator*() const {
		return current;
	}
	constexpr tuple_iterator_on_subset& operator++() {
		LabelType s = N - 1;
        LabelType r = R - 1;
        while(s >= 0 && r >= 0) {
            if (current[r] < s && filter[s]) break;
            else if (current[r] == s) --r;
            --s;
        }
        if (r < 0) {
			if (R > 0) current[0] = N;
			else N = 1;
			return *this;
		}
        s = current[r] + 1;
        while (r < R) {
            while (!filter[s]) ++s;
            current[r] = s;
            ++r;
            ++s;
        }
		return *this;
	}
	constexpr bool operator==(const tuple_iterator_on_subset& other) {
		return (R == other.R) && (N == other.N) && (filter == other.filter) && (
			(current == other.current) ||
			(current[0] == N && other.current[0] == N)
		);
	}
	constexpr bool operator!=(const tuple_iterator_on_subset& other) {
		return (R != other.R) || (N != other.N) || (filter != other.filter) || (
			(current != other.current) &&
			(current[0] != N || other.current[0] != N)
		);
	}
};

}