#pragma once

#include <iostream>
#include <string>
#include <array>
#include <vector>
#include "OMs.hpp"

// =============
// Matroid<R, N>
// =============

// Here "idx >> 5" is just "idx / 32" and "idx & 31" is "idx % 32"

template<int R, int N>
constexpr bool Matroid<R, N>::is_loop(int element) const {
	for (auto i = 0; i < BASE::NR_INT32; i++) {
		if (BASE::bits[i] & RTUPLES::LIST::contained_mask32[element][i])
			return false;
	}
	return true;
}

template<int R, int N>
constexpr bool Matroid<R, N>::is_coloop(int element) const {
	for (auto i = 0; i < BASE::NR_INT32; i++) {
		if (BASE::bits[i] & ~RTUPLES::LIST::contained_mask32[element][i])
			return false;
	}
	return true;
}

template<int R, int N>
constexpr bool Matroid<R, N>::is_loopfree() const {
	for (int e = 0; e < N; ++e) {
		if (is_loop(e)) return false;
	}
	return true;
}

template<int R, int N>
constexpr int Matroid<R, N>::loopcount(int max_nr) const {
	int count = 0;
	for (int e = 0; e < N; ++e) {
		if (is_loop(e)) {
			++count;
			if (count >= max_nr)
				return max_nr;
		}
	}
	return count;
}

template<int R, int N>
template<typename Iterable>
constexpr bool Matroid<R, N>::is_independent(const Iterable& elements) const {
	for (auto i = 0; i < BASE::NR_INT32; ++i) {
		uint32_t bases_with_all_elements = BASE::bits[i];
		for (char element: elements) {
			bases_with_all_elements &= RTUPLES::LIST::contained_mask32[element][i];
		}
		if (bases_with_all_elements) {
			return true;
		}
	}
	return false;
}

template<int R, int N>
template<typename Iterable>
constexpr int Matroid<R, N>::rank(const Iterable& elements, int max_rank) const {
	int rank = 0;
	std::array<uint32_t, BASE::NR_INT32> bases_with_selected_elements;
	std::array<uint32_t, BASE::NR_INT32> bases_with_selected_elements_and_one;
	for (auto i = 0; i < BASE::NR_INT32; ++i) {
		bases_with_selected_elements[i] = BASE::bits[i];
	}
	for (char element: elements) {
		bool is_element_independent_of_selected_elements = false;
		for (auto i = 0; i < BASE::NR_INT32; ++i) {
			bases_with_selected_elements_and_one[i] = 
				bases_with_selected_elements[i]
				& RTUPLES::LIST::contained_mask32[element][i];
			is_element_independent_of_selected_elements =
				is_element_independent_of_selected_elements || 
				bool(bases_with_selected_elements_and_one[i]);
		}
		if (is_element_independent_of_selected_elements) {
			rank += 1;
			if (rank == max_rank) {
				return rank;
			}
			for (auto i = 0; i < BASE::NR_INT32; ++i) {
				bases_with_selected_elements[i] =
					bases_with_selected_elements_and_one[i];
			}
		}
	}
	return rank;
}

template<int R, int N>
template<typename Iterable>
constexpr std::vector<char> Matroid<R, N>::maximal_independent_subset(const Iterable& elements, int max_rank) const {
	int rank = 0;
	std::array<uint32_t, BASE::NR_INT32> bases_with_selected_elements;
	std::array<uint32_t, BASE::NR_INT32> bases_with_selected_elements_and_one;
	std::vector<char> maximal_independent_subset_so_far{};
	for (auto i = 0; i < BASE::NR_INT32; ++i) {
		bases_with_selected_elements[i] = BASE::bits[i];
	}
	for (char element: elements) {
		bool is_element_independent_of_selected_elements = false;
		for (auto i = 0; i < BASE::NR_INT32; ++i) {
			bases_with_selected_elements_and_one[i] = 
				bases_with_selected_elements[i]
				& RTUPLES::LIST::contained_mask32[element][i];
			is_element_independent_of_selected_elements =
				is_element_independent_of_selected_elements || 
				bool(bases_with_selected_elements_and_one[i]);
		}
		if (is_element_independent_of_selected_elements) {
			rank += 1;
			maximal_independent_subset_so_far.push_back(element);
			if (rank == max_rank) {
				return maximal_independent_subset_so_far;
			}
			for (auto i = 0; i < BASE::NR_INT32; ++i) {
				bases_with_selected_elements[i] =
					bases_with_selected_elements_and_one[i];
			}
		}
	}
	return maximal_independent_subset_so_far;
}

template<int R, int N>
template<int target_rank, typename Iterable>
constexpr std::array<char, target_rank> Matroid<R, N>::maximal_independent_subset_of_rank(const Iterable& elements) const {
	int rank = 0;
	std::array<uint32_t, BASE::NR_INT32> bases_with_selected_elements;
	std::array<uint32_t, BASE::NR_INT32> bases_with_selected_elements_and_one;
	std::array<char, target_rank> maximal_independent_subset_so_far{};
	for (auto i = 0; i < BASE::NR_INT32; ++i) {
		bases_with_selected_elements[i] = BASE::bits[i];
	}
	for (char element: elements) {
		bool is_element_independent_of_selected_elements = false;
		for (auto i = 0; i < BASE::NR_INT32; ++i) {
			bases_with_selected_elements_and_one[i] = 
				bases_with_selected_elements[i]
				& RTUPLES::LIST::contained_mask32[element][i];
			is_element_independent_of_selected_elements =
				is_element_independent_of_selected_elements || 
				bool(bases_with_selected_elements_and_one[i]);
		}
		if (is_element_independent_of_selected_elements) {
			maximal_independent_subset_so_far[rank] = element;
			rank += 1;
			if (rank == target_rank) {
				return maximal_independent_subset_so_far;
			}
			for (auto i = 0; i < BASE::NR_INT32; ++i) {
				bases_with_selected_elements[i] =
					bases_with_selected_elements_and_one[i];
			}
		}
	}
	return maximal_independent_subset_so_far;
}

template<int R, int N>
std::ostream& operator<<(std::ostream& os, const Matroid<R, N>& matroid) {
    return os << static_cast<const Matroid<R, N>::BASE&>(matroid);
}

template<int R, int N>
std::ofstream& operator<<(std::ofstream& of, const Matroid<R, N>& matroid) {
    return of << static_cast<const Matroid<R, N>::BASE&>(matroid);
}

template<int R, int N>
std::istream& operator>>(std::istream& is, Matroid<R, N>& matroid) {
    std::string str;
    is >> str;
    matroid.read(str);
    return is;
}

template<int R, int N>
std::ifstream& operator>>(std::ifstream& ifs, Matroid<R, N>& matroid) {
    std::string str;
    ifs >> str;
    matroid.read(str);
    return ifs;
}

// ===============
// Chirotope<R, N>
// ===============

template<int R, int N>
constexpr Matroid<R, N> Chirotope<R, N>::underlying_matroid() const {
    Matroid<R, N> matroid;
    for (auto i = 0; i < BASE::NR_INT32; i++) {
        matroid.bits[i] = BASE::plus[i] | BASE::minus[i];
    }
    return matroid;
}

// Here "idx >> 5" is just "idx / 32" and "idx & 31" is "idx % 32"

template<int R, int N>
constexpr bool Chirotope<R, N>::is_loop(int element) const {
	for (auto i = 0; i < BASE::NR_INT32; i++) {
		if ((BASE::plus[i] | BASE::minus[i]) & RTUPLES::LIST::contained_mask32[element][i])
			return false;
	}
	return true;
}

template<int R, int N>
constexpr bool Chirotope<R, N>::is_coloop(int element) const {
	for (auto i = 0; i < BASE::NR_INT32; i++) {
		if ((BASE::plus[i] | BASE::minus[i]) & ~RTUPLES::LIST::contained_mask32[element][i])
			return false;
	}
	return true;
}

template<int R, int N>
constexpr bool Chirotope<R, N>::is_loopfree() const {
	for (int e = 0; e < N; ++e) {
		if (is_loop(e)) return false;
	}
	return true;
}

template<int R, int N>
constexpr int Chirotope<R, N>::loopcount(int max_nr) const {
	int count = 0;
	for (int e = 0; e < N; ++e) {
		if (is_loop(e)) {
			++count;
			if (count >= max_nr)
				return max_nr;
		}
	}
	return count;
}

template<int R, int N>
template<typename Iterable>
constexpr bool Chirotope<R, N>::is_independent(const Iterable& elements) const {
	for (auto i = 0; i < BASE::NR_INT32; ++i) {
		uint32_t bases_with_all_elements = BASE::plus[i] | BASE::minus[i];
		for (char element: elements) {
			bases_with_all_elements &= RTUPLES::LIST::contained_mask32[element][i];
		}
		if (bases_with_all_elements) {
			return true;
		}
	}
	return false;
}

template<int R, int N>
template<typename Iterable>
constexpr int Chirotope<R, N>::rank(const Iterable& elements, int max_rank) const {
	int rank = 0;
	std::array<uint32_t, BASE::NR_INT32> bases_with_selected_elements;
	std::array<uint32_t, BASE::NR_INT32> bases_with_selected_elements_and_one;
	for (auto i = 0; i < BASE::NR_INT32; ++i) {
		bases_with_selected_elements[i] = BASE::plus[i] | BASE::minus[i];
	}
	for (char element: elements) {
		bool is_element_independent_of_selected_elements = false;
		for (auto i = 0; i < BASE::NR_INT32; ++i) {
			bases_with_selected_elements_and_one[i] = 
				bases_with_selected_elements[i]
				& RTUPLES::LIST::contained_mask32[element][i];
			is_element_independent_of_selected_elements =
				is_element_independent_of_selected_elements || 
				bool(bases_with_selected_elements_and_one[i]);
		}
		if (is_element_independent_of_selected_elements) {
			rank += 1;
			if (rank == max_rank) {
				return rank;
			}
			for (auto i = 0; i < BASE::NR_INT32; ++i) {
				bases_with_selected_elements[i] =
					bases_with_selected_elements_and_one[i];
			}
		}
	}
	return rank;
}

template<int R, int N>
template<typename Iterable>
constexpr std::vector<char> Chirotope<R, N>::maximal_independent_subset(const Iterable& elements, int max_rank) const {
	int rank = 0;
	std::array<uint32_t, BASE::NR_INT32> bases_with_selected_elements;
	std::array<uint32_t, BASE::NR_INT32> bases_with_selected_elements_and_one;
	std::vector<char> maximal_independent_subset_so_far{};
	for (auto i = 0; i < BASE::NR_INT32; ++i) {
		bases_with_selected_elements[i] = BASE::plus[i] | BASE::minus[i];
	}
	for (char element: elements) {
		bool is_element_independent_of_selected_elements = false;
		for (auto i = 0; i < BASE::NR_INT32; ++i) {
			bases_with_selected_elements_and_one[i] = 
				bases_with_selected_elements[i]
				& RTUPLES::LIST::contained_mask32[element][i];
			is_element_independent_of_selected_elements =
				is_element_independent_of_selected_elements || 
				bool(bases_with_selected_elements_and_one[i]);
		}
		if (is_element_independent_of_selected_elements) {
			rank += 1;
			maximal_independent_subset_so_far.push_back(element);
			if (rank == max_rank) {
				return maximal_independent_subset_so_far;
			}
			for (auto i = 0; i < BASE::NR_INT32; ++i) {
				bases_with_selected_elements[i] =
					bases_with_selected_elements_and_one[i];
			}
		}
	}
	return maximal_independent_subset_so_far;
}

template<int R, int N>
template<int target_rank, typename Iterable>
constexpr std::array<char, target_rank> Chirotope<R, N>::maximal_independent_subset_of_rank(const Iterable& elements) const {
	int rank = 0;
	std::array<uint32_t, BASE::NR_INT32> bases_with_selected_elements;
	std::array<uint32_t, BASE::NR_INT32> bases_with_selected_elements_and_one;
	std::array<char, target_rank> maximal_independent_subset_so_far{};
	for (auto i = 0; i < BASE::NR_INT32; ++i) {
		bases_with_selected_elements[i] = BASE::plus[i] | BASE::minus[i];
	}
	for (char element: elements) {
		bool is_element_independent_of_selected_elements = false;
		for (auto i = 0; i < BASE::NR_INT32; ++i) {
			bases_with_selected_elements_and_one[i] = 
				bases_with_selected_elements[i]
				& RTUPLES::LIST::contained_mask32[element][i];
			is_element_independent_of_selected_elements =
				is_element_independent_of_selected_elements || 
				bool(bases_with_selected_elements_and_one[i]);
		}
		if (is_element_independent_of_selected_elements) {
			maximal_independent_subset_so_far[rank] = element;
			rank += 1;
			if (rank == target_rank) {
				return maximal_independent_subset_so_far;
			}
			for (auto i = 0; i < BASE::NR_INT32; ++i) {
				bases_with_selected_elements[i] =
					bases_with_selected_elements_and_one[i];
			}
		}
	}
	return maximal_independent_subset_so_far;
}

template<int R, int N>
constexpr bool Chirotope<R, N>::weak_maps_to(const Matroid<R, N>& matroid) const {
    for (auto i = 0; i < BASE::NR_INT32; i++) {
        if (~(BASE::plus[i] | BASE::minus[i]) & matroid.bits[i]) 
			return false;
    }
    return true;
}

template<int R, int N>
constexpr bool Chirotope<R, N>::OM_weak_maps_to(const Chirotope& chi) const {
    bool is_wm = true;
    for (auto i = 0; i < BASE::NR_INT32; i++) {
        if (~BASE::plus[i] & chi.minus[i] | ~BASE::minus[i] & chi.plus[i]) {
            is_wm = false;
            break;
        }
    }
    if (is_wm) return true;
    return weak_maps_to(chi);
}

template<int R, int N>
constexpr bool Chirotope<R, N>::is_same_OM_as(const Chirotope& chi) const {
    bool is_same = true;
    for (auto i = 0; i < BASE::NR_INT32; i++) {
        if (BASE::plus[i] != chi.minus[i] || BASE::minus[i] != chi.plus[i]) {
            is_same = false;
            break;
        }
    }
    if (is_same) return true;
    return (*this) == chi;
}

template<int R, int N>
std::ostream& operator<<(std::ostream& os, const Chirotope<R, N>& chi) {
    return os << static_cast<const Chirotope<R, N>::BASE&>(chi);
}

template<int R, int N>
std::ofstream& operator<<(std::ofstream& of, const Chirotope<R, N>& chi) {
    return of << static_cast<const Chirotope<R, N>::BASE&>(chi);
}

template<int R, int N>
std::istream& operator>>(std::istream& is, Chirotope<R, N>& chi) {
    std::string str;
    is >> str;
    chi.read(str);
    return is;
}

template<int R, int N>
std::ifstream& operator>>(std::ifstream& ifs, Chirotope<R, N>& chi) {
    std::string str;
    ifs >> str;
    chi.read(str);
    return ifs;
}

// =======================
// LEGACY CODE FROM NEVENA
// =======================
// Ported to better interact with the C++ codebase.

template<int R, int N>
constexpr bool Chirotope<R, N>::is_chirotope() const {
	char i,j,k,l,p,q;
	long long int h;

	for (k=0;k<BASE::NR_INT32;k++)
		if (BASE::plus[k] & BASE::minus[k]) 		//if the same basis is both positive and negative
			return false;
	

	for (k=0;k<BASE::NR_INT32;k++)
		if (BASE::plus[k]!=0 || BASE::minus[k]!=0) 		//(B0)
			break;

	if (k==BASE::NR_INT32)
		return false; 						

	std::array<char,R> x;						//(B2') Lemma 3.5.4
	std::array<char,R> y;					

	char sign,limit_i,limit_j;
	int pi,pj,mi,mj;

	for (k=0;k<BASE::NR_INT32;k++)
	{
		if (k==BASE::NR_INT32)
			limit_i=RTUPLES::NR &31;
		else limit_i=32;
		for (i=0;i<limit_i;i++)
		{
			
			h=(uint32_t)1 <<i;
			pi=BASE::plus[k] & h;	 		//pi!=0 iff \chi(bases[i])=1
			mi = BASE::minus[k] & h;		 	//mi!=0 iff \chi(bases[i])=-1
			if ((pi == 0) && (mi == 0)) 	//\chi(bases[i])=0, we do not have to worry about this basis
				continue;
			
			for (l=k;l<BASE::NR_INT32;l++)
			{
				if (l==k)
					j=i+1;
				else j=0;
				if (l==BASE::NR_INT32-1)
					limit_j=RTUPLES::NR&31;
				else limit_j=32;
				for (;j<limit_j;j++)
				{	
					h=(uint32_t)1<<j;
					pj = BASE::plus[l] & h;
					mj = BASE::minus[l] & h;
					if (pj == 0 && mj == 0)
						continue;
					
					if ((pi && mj) || (mi && pj)) 		 //\chi(x_1,x_2,x_3)* \chi(y_1,y_2,y_3)=-1
						sign = -1;
					else if ((pi && pj) || (mi && mj)) 	 //\chi(x_1,x_2,x_3)* \chi(y_1,y_2,y_3)=1
						sign = 1;

					for (p=0;p<R;p++)			//we have to check B2' for all permutations of x1,x2,x3, but it suffices to have all entries of x in the first position
					{
						for (q=0;q<R;q++)
						{
							x[q]=RTUPLES::LIST::array[i+(k<<5)][q];
							y[q]=RTUPLES::LIST::array[j+(l<<5)][q];
						}	
						x[0]=RTUPLES::LIST::array[i+(k<<5)][p];
						x[p]=RTUPLES::LIST::array[i+(k<<5)][0];
						
						if (p==1)
							sign=-sign;

						if (!b2prime(sign,x,y)) 		//checks B2'
						{
							return false;
						}
						
					}					    
			
				}
			
			}
		}
	}
	return true;
}

template<int R, int N>
constexpr bool Chirotope<R, N>::b2prime(char sign, const std::array<char,R>& X, const std::array<char,R>& Y) const {
	int s1,s2,in1,in2,i,j,q,sx,sy;
	std::array<char,R> x;
	std::array<char,R> y;

	// Bálint: what's the point of this?
	for (i=0;i<R;i++)
	{
		x[i]=X[i];
		y[i]=Y[i];
	}

	for (j=0;j<R;j++)
	{
		for (i=0;i<R;i++)
		{
			x[i]=X[i];
			y[i]=Y[i];
		}
		x[0]=Y[j];
		y[j]=X[0];

		
		//checks \chi(y1,x2,x3)*\chi(x1,y2,y3)
		auto s1_and_idx1 = RTUPLES::sign_and_index_of_unordered(x);
		auto s2_and_idx2 = RTUPLES::sign_and_index_of_unordered(y);
		s1 = s1_and_idx1.first;
		s2 = s2_and_idx2.first;

		if (s1!=0 && s2!=0)			//s1==0 means that two of y1,x2,x3 are the same, its chirotope value is 0 and we want \chi(y1,x2,x3)*\chi(x1,y2,y3)<0
		{
			
			if (axB2(sign,s1,s2,s1_and_idx1.second,s2_and_idx2.second))
			{
				return true;
			}
			
		}
	}	
	
	return false;
}

template<int R, int N>
constexpr bool Chirotope<R, N>::axB2(char sign, char s1, char s2, int in1, int in2) const {
	Chirotope<R, N> M1, M2;
	char i,j,p1,p2;
	long long int i1,i2;
	char res;
	if (s1==1)					 //if sign(y1,x2,x3)=1, then chi(y1,x2,x3) is exactly what we give by plus,minus; 
		M1 = *this;
		/*for (i=0;i<NR_INT32;i++)
		{
			M1.plus[i]=plus[i];
			M1.minus[i]=M.minus[i];
		}*/		
	else						//if sign(y1,x2,x3)=-1, then we have to switch the signs of the bases (chirotopes are alternating)
		M1 = inverse();
		/*for (i=0;i<nr_ints;i++)
		{
			M1.plus[i]=M.minus[i];
			M1.minus[i]=M.plus[i];
		}*/

	if (s2==1)					//the same as for s1
		M2 = *this;
		/*for (i=0;i<nr_ints;i++)
		{
			M2.plus[i]=M.plus[i];
			M2.minus[i]=M.minus[i];
		}*/		
	else						//if sign(y1,x2,x3)=-1, then we have to switch the signs of the bases (chirotopes are alternating)
		M2 = inverse();
		/*for (i=0;i<nr_ints;i++)
		{
			M2.plus[i]=M.minus[i];
			M2.minus[i]=M.plus[i];
		}*/

	

	p1=in1>>5;
	p2=in2>>5;
	in1=in1&31;
	in2=in2&31;

	i1 = 1l<<in1;
	i2 = 1l<<in2;

	
	res=-1;
	if (sign == -1)
		res = (((M1.plus[p1] & i1)&&(M2.minus[p2] & i2))||((M1.minus[p1] & i1)&&(M2.plus[p2] & i2)));		 //\chi(y1,x2,x3)*\chi(x1,y2,y3)
	else if (sign == 1)
		res = (((M1.plus[p1] & i1)&&(M2.plus[p2] & i2))||((M1.minus[p1] & i1)&&(M2.minus[p2] & i2))); 		//\chi(y1,x2,x3)*\chi(x1,y2,y3)	
	
	return res;
}