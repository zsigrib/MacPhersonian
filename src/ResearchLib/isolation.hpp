#pragma once
#include "OMtools.hpp"
#include "research_file_template.hpp"

namespace research {

// Returns whether the given cocircuit constrains the given
// element in the given oriented matroid, i.e. if kernel 
// (without possibly the element) is of rank `R-1`.
template<int R, int N>
bool constrains(const Chirotope<R, N>&, const sign_vector<N>&, char);

// Returns whether there is a non-loop element, distinct form
// the given one, for which none of the given cocircuits (given 
// through pointers) of the given chirotope take `+` or none of
// them take `-`.
//
// If the template parameter `loopfree` is true, it is assumed that
// the given oriented matroid is loopfree.
template<int R, int N, bool loopfree=false>
bool exists_covariant_nonloop(const Chirotope<R, N>&, const std::array<const sign_vector<N>*,R>&, char);

// Returns whether the given element is isolated by any set of
// the given cocircuits (allowing sign reversal in said cocircuits)
// of the given chirotope.
// 
// In other words, the result is `false` if there are `R` cocircuits
// `C[0],...C[1]` in the given vector, together with `R` many signs 
// `s[0],...,s[R-1]` such that `C[i]*s[i]` is not `-` on the given
// element, each of the `C[i]` constrain it (in the sense of 
// `constrains(...)`), and there is no non-loop `w` with the property
// that the `C[i]` is either never `+` or never `-` on `w`.
// Note that this means that loops can easily be isolated.
//
// Changes the sign of certain cocircuits in `cocircuits`.
// If the template parameter `loopfree` is true, it is assumed that
// the given oriented matroid is loopfree.
template<int R, int N, bool loopfree=false>
bool is_isolated(const Chirotope<R, N>&, std::vector<sign_vector<N>>&, char);

}

#include "isolation.cpp"