#pragma once
#include <array>
#include <utility>
#include "mymath.hpp"

namespace Rtuples {

// Iterates over the `R`-element subsets of `0..N-1`
// (also called `R`-tuples) as monotonely increasing
// `R`-long sequences in `0..N-1` in lexicographic order, 
// and returns them as sorted arrays.
//
// The type `LabelType` should be able to store non-negative
// integers up to `N`.
template<typename LabelType, LabelType R>
struct iterator {
	// The current `R`-tuple that the iterator is pointing to.
	std::array<LabelType, R> current;
	// The number of labels, i.e. an `R`-tuple can contain elements
	// in `0..N-1`.
	const LabelType N;
	// Construct an iterator, pointing to the first `R`-tuple
	// `(0,1,2,...,R-1)`, inside the set of `R`-tuples in `0..n-1`.
	constexpr iterator(LabelType n): N(n) {
		if (N < R) {
			current[0] = N;
		} else {
			for (LabelType i = 0; i < R; ++i) {
				current[i] = i;
			}
		}
	}
	// Construct an iterator, pointing to an invalid `R`-tuple,
	// consisting solely of `n`s.
	constexpr iterator(LabelType n, int): N(n) {
		current[0] = N;
	}
	// Construct an iterator, pointing to the specified `R`-tuple.
	constexpr iterator(LabelType n, const std::array<LabelType, R>& target): N(n) {
		current = target;
	}
	// Construct an iterator, pointing to the first `R`-tuple
	// `(0,1,2,...,R-1)`, inside the set of `R`-tuples in `0..n-1`.
	constexpr iterator begin() const {
		return *this;
	}
	// Construct an iterator, pointing to an invalid `R`-tuple, which
	// signifies the end of the iteration.
	constexpr iterator end() const {
		return iterator(N,0);
	}
	// Return the `R`-tuple this iterator is currently pointing to.
	constexpr std::array<LabelType,R> operator*() const {
		return current;
	}
	constexpr iterator& operator++() {
		LabelType r = R - 1;
		while(r > 0 && current[r] >= N-R+r)
			r--;
		if (r == 0 && current[0] >= N-R) {
			current[0] = N;
			return *this;
		}
		LabelType current_r = current[r];
		for (LabelType i = r; i < R; i++)
			current[i] = current_r + (i - r) + 1;
		return *this;
	}
	constexpr bool operator==(const iterator& other) {
		return (N == other.N) && (
			(current[0] == N && other.current[0] == N) || 
			(current == other.current)
		);
	}
	constexpr bool operator!=(const iterator& other) {
		return (N != other.N) || (
			(current[0] != N || other.current[0] != N) && 
			(current != other.current)
		);
	}
};

// Iterates over the `R`-element subsets of `0..N-1`
// (also called `R`-tuples) as monotonely increasing
// `R`-long sequences in `0..N-1` in "reverse lexicographic 
// order", and returns them as sorted arrays. In other
// words if one were to reverse each resulting `R`-tuple
// to instead be monotonely decreasing, then those
// `R`-tuples would be lexicographically sorted.
//
// The type `LabelType` should be able to store non-negative
// integers up to `N`.
template<typename LabelType, LabelType R>
struct iterator_reverse_lexicographic {
	// The current `R`-tuple that the iterator is pointing to.
	std::array<LabelType, R> current;
	// The number of labels, i.e. an `R`-tuple can contain elements
	// in `0..N-1`.
	const LabelType N;
	// Construct an iterator, pointing to the first `R`-tuple
	// `(0,1,2,...,R-1)`, inside the set of `R`-tuples in `0..n-1`.
	constexpr iterator_reverse_lexicographic(LabelType n): N(n) {
		if (N < R) {
			current[0] = N;
		} else {
			for (LabelType i = 0; i < R; ++i) {
				current[i] = i;
			}
		}
	}
	// Construct an iterator, pointing to an invalid `R`-tuple,
	// consisting solely of `n`s.
	constexpr iterator_reverse_lexicographic(LabelType n, int): N(n) {
		current[0] = N;
	}
	// Construct an iterator, pointing to the specified `R`-tuple.
	constexpr iterator_reverse_lexicographic(LabelType n, const std::array<LabelType, R>& target): N(n) {
		current = target;
	}
	// Construct an iterator, pointing to the first `R`-tuple
	// `(0,1,2,...,R-1)`, inside the set of `R`-tuples in `0..n-1`.
	constexpr iterator_reverse_lexicographic begin() const {
		return *this;
	}
	// Construct an iterator, pointing to an invalid `R`-tuple, which
	// signifies the end of the iteration.
	constexpr iterator_reverse_lexicographic end() const {
		return iterator_reverse_lexicographic(N,0);
	}
	// Return the `R`-tuple this iterator is currently pointing to.
	constexpr std::array<LabelType,R> operator*() const {
		return current;
	}
	constexpr iterator_reverse_lexicographic& operator++() {
		LabelType r = 0;
		while(r < R - 1 && current[r] == current[r+1] - 1)
			r++;
		if (r == R - 1 && current[r] >= N-1) {
			current[0] = N;
			return *this;
		}
		++current[r];
		for (LabelType i = 0; i < r; ++i) {
			current[i] = i;
		}
		return *this;
	}
	constexpr bool operator==(const iterator_reverse_lexicographic& other) {
		return (N == other.N) && (
			(current[0] == N && other.current[0] == N) || 
			(current == other.current)
		);
	}
	constexpr bool operator!=(const iterator_reverse_lexicographic& other) {
		return (N != other.N) || (
			(current[0] != N || other.current[0] != N) && 
			(current != other.current)
		);
	}
};

}

// NchooseK<N,k>::array lists all k-element subsets of {0..N-1}.
// Auxulary constants are provided to make querying this list a
// smoother experience. All constants are computed compile-time.
//
// The template parameter `LabelType` should be able to store
// non-negative integers up to and including `N`. `R<=N` is assumed.
// `IndexType` should be able to store non-negative integers up to
// `binomial_coefficient(N,R)`.
template<typename LabelType, LabelType N, LabelType k, typename IndexType = LabelType,
    std::enable_if_t<std::is_integral_v<LabelType>, bool> = true>
struct NchooseK {
    // `(N choose k)=N!/(k!*(N-k)!)`, the number of distinct `k`-element 
    // subsets of `0..N-1`.
    constexpr static const auto NR_RTUPLES = binomial_coefficient<IndexType>(N,k);
    // Number of (unsigned) 32-bit integers needed to have 1 bit for each 
    // `k` element subset of `0..N-1`.
    constexpr static const auto NR_INT32 = division_rounded_up<IndexType>(NR_RTUPLES, 32);
    // Number of bits used in the last 32-bit unsigned integer used to 
    // store the characteristic vector.
    constexpr static const auto NR_REMAINING_BITS = NR_RTUPLES % 32;
    // `array[i][t]` is the `t`th element of the `i`th `k`-element 
    // subset of `0..N-1`. This is computed at compile time.
    constexpr static const auto array{[]() constexpr{
        std::array<std::array<LabelType, k>, NR_RTUPLES> a{}; 
        IndexType idx = 0;
		for (auto Rtuple: Rtuples::iterator<LabelType,k>(N)) {
			a[idx] = Rtuple;
			++idx;
		}
		return a;
    }()};
    // `contained[e][i]` is true if and only if the `i`th `k`-element
    // subset of `0..N-1` contains the element `e`.
    constexpr static const auto contained{[]() constexpr{
        std::array<std::array<bool, NR_RTUPLES>, N> c{}; 
        for (LabelType e = 0; e < k; e++) {
            for (IndexType i = 0; i < NR_RTUPLES; i++) {
                c[array[i][e]][i] = true;
            }
        }
        return c;
    }()};
    // The `r`th lowest bit of `contained_mask32[e][idx]` is 1
    // if and only if the `idx * 32 + r`th `k`-element subset
    // of `0..N-1` contains the element `e`.
    constexpr static const auto contained_mask32{[]() constexpr{
        std::array<std::array<uint32_t, NR_INT32>, N> m{};
        for (LabelType e = 0; e < N; e++) {
            for (IndexType idx = 0; idx < NR_INT32 - 1; idx++) {
                for (auto r = 0; r < 32; r++) {
                    m[e][idx] |= ((uint32_t)contained[e][idx * 32 + r]
                        & (uint32_t)1) << r;
                }
            }
            for (auto r = 0; r < NR_REMAINING_BITS; r++) {
                m[e][NR_INT32 - 1] |= 
                    ((uint32_t)contained[e][(NR_INT32 - 1) * 32 + r]
                    & (uint32_t)1) << r;
            }
        }
        return m;
    }()};
    // If one writes the elements of `0..N-1` not present in the `i`th
    // `k`-tuple right before this `k`-tuple (in increasing order), is the 
    // sign of the resulting permutation of `0..N-1` equal to `-`?
    // This is answered by `is_dual_inverted[i]`.
    constexpr static const auto is_dual_inverted{[]() constexpr{
        std::array<bool, NR_RTUPLES> di{};
        for (IndexType i = 0; i < NR_RTUPLES; i++) {
            IndexType sum = (k * (2*N-k-1)) / 2; // binomial_coefficient(k,2)+k*(N-k)
            for (LabelType t = 0; t < k; t++) {
                sum += (array[i][t]) % 2;
            }
            di[i] = sum % 2;
        }
        return di;
    }()};
    // The `r`th lowest bit of `is_dual_inverted[idx]` is 1 if and only
    // if `is_dual_inverted[idx * 32 + r]` is true.
    constexpr static const auto is_dual_inverted_mask32{[]() constexpr{
        std::array<uint32_t, NR_INT32> di{};
        for (IndexType idx = 0; idx < NR_INT32 - 1; idx++) {
            for (auto r = 0; r < 32; r++) {
                di[idx] |= ((uint32_t)is_dual_inverted[idx * 32 + r]
                    & (uint32_t)1) << r;
            }
        }
        for (auto r = 0; r < NR_REMAINING_BITS; r++) {
            di[NR_INT32 - 1] |= ((uint32_t)is_dual_inverted[(NR_INT32 - 1) * 32 + r]
                &(uint32_t)1) << r;
        }
        return di;
    }()};

	constexpr static const auto char_vector_of_Rtuple{[]() constexpr {
		constexpr auto N_NR_INT32 = division_rounded_up<LabelType>(N, 32);
		std::array<std::array<uint32_t, N_NR_INT32>, NR_RTUPLES> cv{};
		for (IndexType idx = 0; idx < NR_RTUPLES; ++idx) {
			for (LabelType t = 0; t < k; ++t) {
				LabelType entry = array[idx][t];
				cv[idx][entry >> 5] |= (uint32_t)1 << (entry & 31);
			}
		}
		return cv;
	}()};
};

// ================
// MANAGING RTUPLES
// ================

// NchooseK provided a way to assign to each index an (orderd) 
// Rtuple. This namespace provides multiple ways to bring R-tuples 
// into canonical (ordered) form, find their indices, and determine
// the sign of their permutations.
//
// Use the type `RTUPLES` to configure all template parameters once
// and for all.
namespace Rtuples {

// Nevena's legacy code for managing R-tuples. Sometimes it is
// slower, sometimes it is faster than the other solutions.
namespace legacy {
//NOLINTBEGIN

// Sorts the given array, and returns the sign of the permutation
// used for sorting. If the array does not consist of distinct 
// elements, then leaves the array in a valid but unspecified state,
// and returns 0.
//
// The template parameter `LabelType` should be able to store
// non-negative integers up to and including `N`. `R<=N` is assumed.
// `SignType` should be able to store `-1`, `0`, or `1`.
template<typename LabelType, LabelType R, LabelType N, typename SignType=int>
constexpr SignType sort_array(std::array<LabelType, R>& a) {
    LabelType p,q;
	if constexpr (R==3)			//quick code for three integers
		if(a[1]<a[0])
			if(a[2]<a[1])
			{
				p=a[2];
				a[2]=a[0];
				a[0]=p;
				return -1;
			}
			else	if (a[2]==a[1])
					return 0;
				else	if (a[0]<a[2])
					{
						p=a[1];
						a[1]=a[0];
						a[0]=p;
						return -1;
					}
					else	if (a[0]==a[2])
							return 0;
						else    {
								p=a[0];
								a[0]=a[1];
								a[1]=a[2];
								a[2]=p;
								return 1;
							}		
		else    if (a[1]==a[0])
				return 0;
			else    if (a[2]<a[1])
					if(a[0]<a[2])
					{
						p=a[1];
						a[1]=a[2];
						a[2]=p;
						return -1;
					}
					else    if (a[0]==a[2])
							return 0;
						else    {
								p=a[0];
								a[0]=a[2];
								a[2]=a[1];
								a[1]=p;
								return 1;
							}
				else    if (a[2]==a[1])
						return 0;
					else    return 1;
	

	if constexpr (R==4)						//quick code for four integers
		if (a[1]<a[0])
			if (a[2]<a[1])
				if (a[3]<a[2])
				{
					p=a[0];
					a[0]=a[3];
					a[3]=p;
					p=a[1];
					a[1]=a[2];
					a[2]=p;
					return 1;
				}
				else	if (a[3]==a[2])
						return 0;
					else 	if (a[3]<a[1])
						{
							p=a[0];
							a[0]=a[2];
							a[2]=a[1];
							a[1]=a[3];
							a[3]=p;
							return -1;
						}
						else 	if (a[1]==a[3])
								return 0;
							else	if (a[3]<a[0])
								{
									p=a[0];
									a[0]=a[2];
									a[2]=a[3];
									a[3]=p;
									return 1;
								}
								else if (a[3]==a[0])
									return 0;
								else 	
								{
									p=a[0];
									a[0]=a[2];
									a[2]=p;
									return -1;
								}
			else	if (a[2]==a[1])
					return 0;
				else	if (a[3]<a[1])
						if (a[2]<a[0])
						{
							p=a[0];
							a[0]=a[3];
							a[3]=p;
							return -1;
						}
						else	if (a[2]==a[0])
								return 0;
							else
							{
								p=a[0];
								a[0]=a[3];
								a[3]=a[2];
								a[2]=p;
								return 1;
							}
					else 	if (a[3]==a[1])
							return 0;
						else	if (a[3]<a[2])
								if (a[2]<a[0])
								{
									p=a[0];
									a[0]=a[1];
									a[1]=a[3];
									a[3]=p;
									return 1;
								}
								else 	if (a[2]==a[0])
										return 0;
									else 	if (a[3]<a[0])
										{
											p=a[0];
											a[0]=a[1];
											a[1]=a[3];
											a[3]=a[2];
											a[2]=p;
											return -1;
										}
										else	if (a[3]==a[0])
												return 0;
											else
											{
												p=a[0];
												a[0]=a[1];
												a[1]=p;
												p=a[2];
												a[2]=a[3];
												a[3]=p;
												return 1;
											}							
						
							else	if (a[3]==a[2])
									return 0;
								else 	if (a[0]<a[2])
									{
										p=a[0];
										a[0]=a[1];
										a[1]=p;
										return -1;
									}
									else 	if (a[2]==a[0])
											return 0;
										else	if (a[0]<a[3])
											{
												p=a[0];
												a[0]=a[1];
												a[1]=a[2];
												a[2]=p;
												return 1;
											}
											else	if (a[0]==a[3])
													return 0;
												else
												{
													p=a[0];
													a[0]=a[1];
													a[1]=a[2];
													a[2]=a[3];
													a[3]=p;
													return -1;
												}
		else 	if (a[1]==a[0])
				return 0;
			else 	if (a[3]<a[1])
					if (a[0]<a[3])
						if (a[2]<a[0])
						{
							p=a[0];
							a[0]=a[2];
							a[2]=a[3];
							a[3]=a[1];
							a[1]=p;
							return -1;
						}
						else	if (a[2]==a[0])
								return 0;
							else	if (a[2]<a[3])
								{
									p=a[1];
									a[1]=a[2];
									a[2]=a[3];
									a[3]=p;
									return 1;
								}
								else	if (a[2]==a[3])
										return 0;
									else	if (a[2]<a[1])
										{
											p=a[1];
											a[1]=a[3];
											a[3]=p;
											return -1;
										}
										else	if (a[2]==a[1])
												return 0;
											else
											{
												p=a[1];
												a[1]=a[3];
												a[3]=a[2];
												a[2]=p;
												return 1;
											}			
					else	if (a[0]==a[3])
							return 0;
						else 	if (a[2]<a[3])
							{
								p=a[0];
								a[0]=a[2];
								a[2]=p;
								p=a[1];
								a[1]=a[3];
								a[3]=p;
								return 1;
							}
							else	if (a[2]==a[3])
									return 0;
								else	if (a[2]<a[0])
									{
										p=a[0];
										a[0]=a[3];
										a[3]=a[1];
										a[1]=a[2];
										a[2]=p;
										return -1;
									}
									else	if (a[2]==a[0])	
											return 0;
										else	if (a[2]<a[1])
											{
												p=a[0];
												a[0]=a[3];
												a[3]=a[1];
												a[1]=p;
												return 1;
											}
											else 	if (a[2]==a[1])
													return 0;
												else
												{
													p=a[0];
													a[0]=a[3];
													a[3]=a[2];
													a[2]=a[1];
													a[1]=p;
													return -1;
												}
				else	if (a[3]==a[1])
						return 0;
					else	if (a[2]<a[0])
						{
							p=a[0];
							a[0]=a[2];
							a[2]=a[1];
							a[1]=p;
							return 1;
						}
						else 	if (a[2]==a[0])
								return 0;
							else	if (a[2]<a[1])
								{
									p=a[1];
									a[1]=a[2];
									a[2]=p;
									return -1;
								}
								else 	if (a[2]==a[1])
										return 0;
									else	if (a[2]<a[3])
											return 1;
										else	if (a[2]==a[3])
												return 0;
											else
											{
												p=a[2];
												a[2]=a[3];
												a[3]=p;
												return -1;
											}
	

	LabelType sorted,min,ind_min,i,j;				//slow code, but works for all R
    SignType sign;

	sorted=0;
	sign=1;
	while (sorted<R)
	{
		min=N;
		for (i=sorted;i<R;i++)
			if (a[i]<min)
			{
				min=a[i];
				ind_min=i;
			}
		if (ind_min!=sorted)
		{
			j=a[sorted];
			a[sorted]=min;
			a[ind_min]=j;
			sign=-sign;
		}
		sorted++;
	}

	for (i=1;i<R;i++)
		if (a[i-1]==a[i])
			return 0;	
	return sign;
}

// Determines the index of an `R`-tuple in the list of 
// `R`-tuples with entries in `0..N-1` given in 
// `NchooseK<LabelType,N,R>::array`. The input `R`-tuple
// is assumed to be sorted.
//
// The template parameter `LabelType` should be able to store
// non-negative integers up to and including `N`. `R<=N` is assumed.
// `IndexType` should be able to store non-negative integers up to
// `binomial_coefficient(N,R)`.
template<typename LabelType, LabelType R, LabelType N, typename IndexType=LabelType>
constexpr IndexType index_of(const std::array<LabelType, R>& a) {
	if constexpr (R==3 && N<=9)		//quick code for common used parameters
	{
		if constexpr (N==9)
		{	
			if (a[0]==0)
				if (a[1]==1)
					return a[2]-IndexType(2);
				else if (a[1]==2)
					return a[2]+IndexType(4);
				else if (a[1]==3)
					return a[2]+IndexType(9);
				else if (a[1]==4)
					return a[2]+IndexType(13);
				else if (a[1]==5)
					return a[2]+IndexType(16);			
				else	return IndexType(12)+a[1]+a[2];
			if (a[0]==1)
				if (a[1]==2)
					return IndexType(25)+a[2];
				else if (a[1]==3)
					return IndexType(30)+a[2];
				else if (a[1]==4)
					return IndexType(34)+a[2];
				else if (a[1]==5)
					return IndexType(37)+a[2];
				else    return IndexType(33)+a[1]+a[2];
			if (a[0]==2)
				if (a[1]==3)
					return IndexType(45)+a[2];
				else if (a[1]==4)
					return IndexType(49)+a[2];
				else if (a[1]==5)
					return IndexType(52)+a[2];
				else    return IndexType(48)+a[1]+a[2];
			if (a[0]==3)
				if (a[1]==4)
					return IndexType(59)+a[2];
				else if (a[1]==5)
					return IndexType(62)+a[2];
				else    return IndexType(58)+a[1]+a[2];
			if (a[0]==4)
				if (a[1]==5)
					return IndexType(68)+a[2];
				else return IndexType(64)+a[1]+a[2];
			return IndexType(62)+a[0]+a[1]+a[2];
		}
		if constexpr (N==8)
		{	
			if (a[0]==0)
				if (a[1]==1)
					return IndexType(-2)+a[2];
				else if (a[1]==2)
					return IndexType(3)+a[2];
				else if (a[1]==3)
					return IndexType(7)+a[2];
				else if (a[1]==4)
					return IndexType(10)+a[2];			
				else	return IndexType(7)+a[1]+a[2];
			if (a[0]==1)
				if (a[1]==2)
					return IndexType(18)+a[2];
				else if (a[1]==3)
					return IndexType(22)+a[2];
				else if (a[1]==4)
					return IndexType(25)+a[2];
				else    return IndexType(22)+a[1]+a[2];
			if (a[0]==2)
				if (a[1]==3)
					return IndexType(32)+a[2];
				else if (a[1]==4)
					return IndexType(35)+a[2];
				else    return IndexType(32)+a[1]+a[2];
			if (a[0]==3)
				if (a[1]==4)
					return IndexType(41)+a[2];
				else return IndexType(38)+a[1]+a[2]; 
			return IndexType(37)+a[0]+a[1]+a[2];
		}
		if constexpr (N==7)
		{
			if (a[0]==0)
				if (a[1]==1)
					return IndexType(-2)+a[2];
				else if (a[1]==2)
					return IndexType(2)+a[2];
				else if (a[1]==3)
					return IndexType(5)+a[2];			
				else	return IndexType(3)+a[1]+a[2];
			if (a[0]==1)
				if (a[1]==2)
					return IndexType(12)+a[2];
				else if (a[1]==3)
					return IndexType(15)+a[2];
				else    return IndexType(13)+a[1]+a[2];
			if (a[0]==2)
				if (a[1]==3)
					return IndexType(21)+a[2];
				else    return IndexType(19)+a[1]+a[2]; 
			return IndexType(19)+a[0]+a[1]+a[2];
		}
		
		if constexpr (N==6)
		{
			if (a[0]==0)
				if (a[1]==1)
					return IndexType(-2)+a[2];
				else if (a[1]==2)
					return IndexType(1)+a[2];
				else	return IndexType(a[1])+a[2];
			if (a[0]==1)
				if (a[1]==2)
					return IndexType(7)+a[2];
				else    return IndexType(6)+a[1]+a[2]; 
			return IndexType(7)+a[0]+a[1]+a[2];
		
		}
		
		if constexpr (N==5)
		{
			if (a[0]==0)
				if (a[1]==1)
					return IndexType(-2)+a[2];
				else	return IndexType(-2)+a[1]+a[2];
			return IndexType(a[0])+a[1]+a[2];		
		}
		if constexpr (N==4)
		{
			return IndexType(-3)+a[0]+a[1]+a[2];
		}
	}
	if constexpr (R==4 && N<=9)
	{
		if constexpr (N==9)
		{
			if (a[0]==0)
				if (a[1]==1)
					if (a[2]==2)
						return IndexType(-3)+a[3];
					else if (a[2]==3)
						return IndexType(2)+a[3];
					else if (a[2]==4)
						return IndexType(6)+a[3];
					else if (a[2]==5)
						return IndexType(9)+a[3];
					else return IndexType(5)+a[2]+a[3];
				else if (a[1]==2)
					if (a[2]==3)
						return IndexType(17)+a[3];
					else if (a[2]==4)
						return IndexType(21)+a[3];
					else if (a[2]==5)
						return IndexType(24)+a[3];
					else return IndexType(20)+a[2]+a[3];
				else if (a[1]==3)
					if (a[2]==4)
						return IndexType(31)+a[3];
					else if (a[2]==5)
						return IndexType(34)+a[3];
					else return IndexType(30)+a[2]+a[3];
				else if (a[1]==4)
					if (a[2]==5)
						return IndexType(40)+a[3];
					else return IndexType(36)+a[2]+a[3];
				else return IndexType(34)+a[1]+a[2]+a[3];
			if (a[0]==1)
				if (a[1]==2)
					if (a[2]==3)
						return IndexType(52)+a[3];
					else if (a[2]==4)
						return IndexType(56)+a[3];
					else if (a[2]==5)
						return IndexType(59)+a[3];
					else return IndexType(55)+a[2]+a[3];
				else if (a[1]==3)
					if (a[2]==4)
						return IndexType(66)+a[3];
					else if (a[2]==5)
						return IndexType(69)+a[3];
					else return IndexType(65)+a[2]+a[3];
				else if (a[1]==4)
					if (a[2]==5)
						return IndexType(75)+a[3];
					else return IndexType(71)+a[2]+a[3];
				else return IndexType(69)+a[1]+a[2]+a[3];
			if (a[0]==2)
				if (a[1]==3)
					if (a[2]==4)
						return IndexType(86)+a[3];
					else if (a[2]==5)
						return IndexType(89)+a[3];
					else return IndexType(85)+a[2]+a[3];
				else if (a[1]==4)
					if (a[2]==5)
						return IndexType(95)+a[3];
					else return IndexType(91)+a[2]+a[3];
				else return IndexType(89)+a[1]+a[2]+a[3];
			if (a[0]==3)
				if (a[1]==4)
					if (a[2]==5)
						return IndexType(105)+a[3];
					else return IndexType(101)+a[2]+a[3];
				else return IndexType(99)+a[1]+a[2]+a[3];
			return IndexType(99)+a[0]+a[1]+a[2]+a[3];
		}

		if constexpr (N==8)
		{
			if (a[0]==0)
				if (a[1]==1)
					if (a[2]==2)
						return IndexType(-3)+a[3];
					else if (a[2]==3)
						return IndexType(1)+a[3];
					else if (a[2]==4)
						return IndexType(4)+a[3];
					else return IndexType(1)+a[2]+a[3];
				else if (a[1]==2)
					if (a[2]==3)
						return IndexType(11)+a[3];
					else if (a[2]==4)
						return IndexType(14)+a[3];
					else return IndexType(11)+a[2]+a[3];
				else if (a[1]==3)
					if (a[2]==4)
						return IndexType(20)+a[3];
					else return IndexType(17)+a[2]+a[3];
				else return IndexType(16)+a[1]+a[2]+a[3];
			if (a[0]==1)
				if (a[1]==2)
					if (a[2]==3)
						return IndexType(31)+a[3];
					else if (a[2]==4)
						return IndexType(34)+a[3];
					else return IndexType(31)+a[2]+a[3];
				else if (a[1]==3)
					if (a[2]==4)
						return IndexType(40)+a[3];
					else return IndexType(37)+a[2]+a[3];
				else return IndexType(36)+a[1]+a[2]+a[3];
			if (a[0]==2)
				if (a[1]==3)
					if (a[2]==4)
						return IndexType(50)+a[3];
					else return IndexType(47)+a[2]+a[3];
				else return IndexType(46)+a[1]+a[2]+a[3];
			return IndexType(47)+a[0]+a[1]+a[2]+a[3];
		}

		if constexpr (N==7)
		{
			if (a[0]==0)
				if (a[1]==1)
					if (a[2]==2)
						return IndexType(-3)+a[3];
					else if (a[2]==3)
						return IndexType(a[3]);
					else return IndexType(-2)+a[2]+a[3];
				else if (a[1]==2)
					if (a[2]==3)
						return IndexType(6)+a[3];
					else return IndexType(4)+a[2]+a[3];
				else return IndexType(4)+a[1]+a[2]+a[3];
			if (a[0]==1)
				if (a[1]==2)
					if (a[2]==3)
						return IndexType(16)+a[3];
					else return IndexType(14)+a[2]+a[3];
				else return IndexType(14)+a[1]+a[2]+a[3];
			return IndexType(16)+a[0]+a[1]+a[2]+a[3];
		}
			
		if constexpr (N==6)
		{
			if (a[0]==0)
				if (a[1]==1)
					if (a[2]==2)
						return IndexType(-3)+a[3];
					else return IndexType(-4)+a[2]+a[3];
				else return IndexType(-3)+a[1]+a[2]+a[3];
			return IndexType(a[0])+a[1]+a[2]+a[3];
		}			
		if constexpr (N==5)
			return IndexType(-6)+a[0]+a[1]+a[2]+a[3];
	}
	
	
	IndexType l,u,m;						//code that works for every parameter, slow...
	LabelType i;
	l=0;
	u=binomial_coefficient<IndexType>(N,R)-1;
	m=(l+u)>>1;
	// TODO: this line has to be rewritten once NchooseK also gets
	// two integral type parameters
	const auto& bases = NchooseK<LabelType, N, R>::array;

	for (i=0;i<R;i++)
	{
		while (a[i]!=bases[m][i])
			if (a[i]>bases[m][i])
			{
				l=m+1;
				if (l==u)
					return l;
				m=(u+l)>>1;
				
			}
			else if (a[i]<bases[m][i])
			{
				u=m-1;
				if (l==u)
					return l;
				m=(u+l)>>1;
			}
		l=m;
		u=m;
		while (l>=0 && a[i]==bases[l][i])
			l--;
		l++;
		while (u< binomial_coefficient<IndexType>(N,R) && a[i]==bases[u][i])
			u++;
		u--;
		if (l==u)
			return l;
		m=(u+l)>>1;
	}
	throw std::invalid_argument("We shouldn't have reached this line in legacy::index_of.");
}

//NOLINTEND
}

// Given any `R`-tuple of possibly non-distinct and unordered elements,
// `unordered_to_index_of_ordered_TABLE::lookup(...)` can be used
// to look up the index of the sorted version of this `R`-tuple
// inside the lexicographic ordering of `R`-element subsets of
// `0..N-1`, given in `NchooseK::array`. This is achieved by
// precomputing all these indices for all possible inputs at
// compile time, and storing them in a flattened (1D) array
// (`flattened_table`). `...::flattened_index(...)` can be used
// to find the corresponding index to an `R`-tuple inside this
// flattened array.
//
// Note that `flattened_index(...)` is identical to
// `unordered_sign_of_sorting_TABLE::flattened_index(...)`.
//
// The template parameter `LabelType` should be able to store
// non-negative integers up to (=and including) `N`. `R<=N` is assumed.
// `IndexType` should be able to store non-negative integers up to
// `binomial_coefficient(N,R)`. `FlatIndexType` should be able to
// store non-negative integers up to the `R`th power of `N`.
template<typename LabelType, LabelType R, LabelType N, typename IndexType=LabelType, typename FlatIndexType=IndexType>
struct unordered_to_index_of_ordered_TABLE {
	// `flattened_table[flattened_index(Rtuple)]` is the index of
	// the ordered version of `Rtuple` inside `NchooseK::array`.
	// If `Rtuple` does not consist of distinct elements, then this
	// value is `0`.
	constexpr static const auto flattened_table{[]() constexpr {
		std::array<IndexType,power<FlatIndexType>(N,R)> table{};
        std::array<LabelType, R> Rtuple;
        for (FlatIndexType i = 0; i < power<FlatIndexType>(N,R); ++i) {
            FlatIndexType remainder = i;
            for (LabelType t = 0; t < R; ++t) {
                Rtuple[R-1-t] = remainder % N;
                remainder /= N;
            }
            int sign = legacy::sort_array<LabelType,R,N>(Rtuple);
            if (sign != 0) {
                table[i] = legacy::index_of<LabelType,R,N,IndexType>(Rtuple);
            }
        }
        return table;
	}()};
	// Given an `R`-tuple of not necessarily distinct elements in `0..N-1`
	// in not necessarily increasing order, find the index of the corresponding
	// entry in `flattened_table`.
	constexpr static FlatIndexType flattened_index(const std::array<LabelType,R>& Rtuple) {
		FlatIndexType idx = 0;
        for (LabelType t = 0; t < R; ++t) {
            idx = idx * N + Rtuple[t];
        }
		return idx;
	}
	// `lookup(Rtuple)` is the index of the ordered version of `Rtuple` 
	// inside `NchooseK::array`. If `Rtuple` does not consist of distinct 
	// elements, then this value is `0`.
	constexpr static IndexType lookup(const std::array<LabelType,R>& Rtuple) {
        return flattened_table[flattened_index(Rtuple)];
	}
};

// Given any `R`-tuple of possibly non-distinct and unordered elements,
// `unordered_to_sign_of_sorting_TABLE::lookup(...)` can be used
// to look up the sign of the permutation that must be used in order
// to sort an `R`-tuple consisting of entries in `0..N-1`. This is 
// achieved by precomputing all these signs for all possible inputs at
// compile time, and storing them in a flattened (1D) array
// (`flattened_table`). `...::flattened_index(...)` can be used
// to find the corresponding sign to an `R`-tuple inside this
// flattened array.
//
// Note that `flattened_index(...)` is identical to
// `unordered_to_index_of_ordered_TABLE::flattened_index(...)`.
//
// The template parameter `LabelType` should be able to store
// non-negative integers up to (=and including) `N`. `R<=N` is assumed.
// `SignType` should be able to store the values `0`, `1` and `-1`. 
// `FlatIndexType` should be able to store non-negative integers up to 
// the `R`th power of `N`.
template<typename LabelType, LabelType R, LabelType N, typename FlatIndexType=LabelType, typename SignType=int>
struct unordered_to_sign_of_sorting_TABLE {
	// `flattened_table[flattened_index(Rtuple)]` is the sign of
	// the permutation that must be used in order to sort `Rtuple`.
	// If `Rtuple` does not consist of distinct elements, then this
	// value is `0`.
	constexpr static const auto flattened_table{[]() constexpr{
        std::array<SignType,power<FlatIndexType>(N,R)> table{};
        std::array<LabelType, R> Rtuple{};
        for (FlatIndexType i = 0; i < power<FlatIndexType>(N,R); ++i) {
            FlatIndexType remainder = i;
            for (LabelType t = 0; t < R; ++t) {
                Rtuple[R-1-t] = remainder % N;
                remainder /= N;
            }
            table[i] = legacy::sort_array<LabelType,R,N,SignType>(Rtuple);
        }
        return table;
    }()};
	// Given an `R`-tuple of not necessarily distinct elements in `0..N-1`
	// in not necessarily increasing order, find the index of the corresponding
	// entry in `flattened_table`.
	constexpr static FlatIndexType flattened_index(const std::array<LabelType,R>& Rtuple) {
		FlatIndexType idx = 0;
        for (LabelType t = 0; t < R; ++t) {
            idx = idx * N + Rtuple[t];
        }
		return idx;
	}
	// `lookup(Rtuple)` is the sign of the permutation that must be used in
	// order to sort `Rtuple`. If `Rtuple` does not consist of distinct 
	// elements, then this value is `0`.
	constexpr static SignType lookup(const std::array<LabelType,R>& Rtuple) {
        return flattened_table[flattened_index(Rtuple)];
	}
};

// Given any `R`-tuple of ordered and distinct elements in `0..N-1`,
// `ordered_to_index_SMALL_TABLE::lookup(...)` can be used
// to look up the index of this `R`-tuple inside the lexicographic 
// ordering of `R`-element subsets of `0..N-1`, given in 
// `NchooseK::array`. This is achieved by precomputing a helper
// table of size `N*R*N`, and summing `R` specific entries of this
// table specified (somehow) by the `R`-tuple. 
//
// The template parameter `LabelType` should be able to store
// non-negative integers up to (=and including) `N`. `R<=N` is assumed.
// `IndexType` should be able to store non-negative integers up to
// `binomial_coefficient(N,R)`. `FlatIndexType` should be able to
// store non-negative integers up to `N*R*N`.
template<typename LabelType, LabelType R, LabelType N, typename IndexType=LabelType, typename FlatIndexType=IndexType>
struct ordered_to_index_SMALL_TABLE {
	// This is a flattened version of an `N`×`R`×`N` table,
	// where the entry at `(n,r,x)` tells how many ordered
	// `r`-tuples are there with distinct elements in `0..n-1`
	// whose first element is strictly smaller than `x`.
	constexpr static const auto flattened_helper_table{[]() constexpr{
        std::array<IndexType, FlatIndexType(N)*R*N> table{};
        for (LabelType n = 0; n < N; ++n) {
            for (LabelType r = 0; r < R; ++r) {
                for (LabelType x = 0; x <= n; ++x) {
                    IndexType total = 0;
                    for (LabelType i = 0; i < x; ++i) {
                        total += binomial_coefficient(n-i,r);
                    }
                    table[FlatIndexType(n)*R*N+FlatIndexType(r)*N+x] = total;
                }
            }
        }
        return table;
    }()};
	// Computes the index of `Rtuple` inside the lexicographic ordering
	// of `R`-tuples with entries in `0..N-1`, as they are listed in
	// `NchooseK::array`. The time complexity is O(`R`), and `R` lookups
	// into `flattened_helper_table` are taken.
	constexpr static IndexType lookup(const std::array<LabelType, R>& Rtuple) {
        const std::array<LabelType, R>& a = Rtuple;
		IndexType sum = 0;
        LabelType last_a = -1;
        for (LabelType j = 0; j < R; ++j) {
            sum += flattened_helper_table[
				FlatIndexType(N-2-last_a)*R*N
				+FlatIndexType(R-1-j)*N
				+a[j]-last_a-1];
            last_a = a[j];
        }
        return sum;
    }
};

// Find the index of a given ordered `R`-tuple of distinct elements
// in `0..N-1` inside the lexicographic ordering of `R`-element
// subsets of `0..N-1`, as they are listed in `NchooseK<LabelType,R,N>::array`.
//
// The template parameter `LabelType` should be able to store
// non-negative integers up to (=and including) `N`. `R<=N` is assumed.
// `IndexType` should be able to store non-negative integers up to
// `binomial_coefficient(N,R)`. `FlatIndexType` should be able to
// store large integers, depending on the specific implementation
// (currently: up to the `R`th power of `N`).
template<typename LabelType, LabelType R, LabelType N, typename IndexType=LabelType, typename FlatIndexType=IndexType>
constexpr IndexType index_of_ordered(const std::array<LabelType, R>& Rtuple) {
	return unordered_to_index_of_ordered_TABLE<LabelType,R,N,IndexType,FlatIndexType>::lookup(Rtuple);
	//return ordered_to_index_SMALL_TABLE<LabelType,R,N,IndexType,FlatIndexType>::lookup(Rtuple);
	//return legacy::index_of<LabelType,R,N,IndexType>(Rtuple);
}
// Find the index of a given not necessarily ordered `R`-tuple of 
// distinct elements in `0..N-1` inside the lexicographic ordering 
// of `R`-element subsets of `0..N-1`, as they are listed in 
// `NchooseK<LabelType,R,N>::array`.
//
// The template parameter `LabelType` should be able to store
// non-negative integers up to (=and including) `N`. `R<=N` is assumed.
// `IndexType` should be able to store non-negative integers up to
// `binomial_coefficient(N,R)`. `FlatIndexType` should be able to
// store large integers, depending on the specific implementation
// (currently: up to the `R`th power of `N`).
template<typename LabelType, LabelType R, LabelType N, typename IndexType=LabelType, typename FlatIndexType=IndexType>
constexpr IndexType index_of_unordered(const std::array<LabelType, R>& Rtuple) {
	return unordered_to_index_of_ordered_TABLE<LabelType,R,N,IndexType,FlatIndexType>::lookup(Rtuple);
	//std::array<LabelType, R> copy = Rtuple;
	//legacy::sort_array(copy);
	//return legacy::index_of<LabelType,R,N,IndexType>(copy);
}
// Find the sign of the permutation that must be used in order to
// sort a given `R`-tuple of distinct elements in `0..N-1`.
// of `R`-element subsets of `0..N-1`, as they are listed in 
// `NchooseK<LabelType,R,N>::array`.
//
// The template parameter `LabelType` should be able to store
// non-negative integers up to (=and including) `N`. `R<=N` is assumed.
// `FlatIndexType` should be able to store large integers, depending on 
// the specific implementation (currently: up to the `R`th power of `N`).
// `SignType` should be able to store the numbers `0`, `1`, and `-1`.
template<typename LabelType, LabelType R, LabelType N, typename FlatIndexType=LabelType, typename SignType=int>
constexpr SignType sign_of_sorting(const std::array<LabelType, R>& Rtuple) {
	return unordered_to_sign_of_sorting_TABLE<LabelType,R,N,FlatIndexType,SignType>::lookup(Rtuple);
	//std::array<LabelType, R> copy = Rtuple;
	//return legacy::sort_array<LabelType,R,N,SignType>(copy);
}
// Return a pair of `sign_of_sorting(Rtuple)` and
// `index_of_unordered(Rtuple)` with the appropriate template
// parameters.
template<typename LabelType, LabelType R, LabelType N, typename IndexType, typename FlatIndexType=IndexType, typename SignType=int>
constexpr std::pair<SignType, IndexType> sign_and_index_of_unordered(const std::array<LabelType, R>& Rtuple) {
	FlatIndexType idx = unordered_to_index_of_ordered_TABLE<LabelType,R,N,IndexType,FlatIndexType>::flattened_index(Rtuple);
	return {
		unordered_to_sign_of_sorting_TABLE<LabelType,R,N,FlatIndexType,SignType>::flattened_table[idx],
		unordered_to_index_of_ordered_TABLE<LabelType,R,N,IndexType,FlatIndexType>::flattened_table[idx]
	};
	//std::array<LabelType, R> copy = Rtuple;
	//return {
	//	legacy::sort_array<LabelType,R,N,SignType>(copy),
	//	legacy::index_of<LabelType,R,N,IndexType>(copy)
	//};
}
// Return a sorted copy of the given `R`-tuple with elements in
// `0..N-1`.
//
// The template parameter `LabelType` should be able to store
// non-negative integers up to (=and including) `N`. `R<=N` is assumed.
// `IndexType` should be able to store non-negative integers up to
// `binomial_coefficient(N,R)`. `FlatIndexType` should be able to
// store non-negative integers up to the `R`th power of `N`.
template<typename LabelType, LabelType R, LabelType N, typename IndexType=LabelType, typename FlatIndexType=LabelType>
constexpr std::array<LabelType, R> sorted(const std::array<LabelType, R>& Rtuple) {
	return NchooseK<LabelType,R,N,IndexType>::array[
		unordered_to_index_of_ordered_TABLE<LabelType,R,N,IndexType,FlatIndexType>::lookup(Rtuple)
	];
	//std::array<LabelType, R> copy = Rtuple;
	//legacy::sort_array<LabelType,R,N,SignType>(copy);
	//return copy; 
}
// Sort the given `R`-tuple with elements in `0..N-1`, and return
// the sign of the permutation used. If there are two equal elements
// in the tuple, then leave the tuple in a valid but unspecified state,
// and return 0.
//
// The template parameter `LabelType` should be able to store
// non-negative integers up to (=and including) `N`. `R<=N` is assumed.
// `IndexType` should be able to store non-negative integers up to
// `binomial_coefficient(N,R)`. `FlatIndexType` should be able to
// store non-negative integers up to the `R`th power of `N`.
// `SignType` should be able to store the numbers `0`, `-1`, and `1`.
template<typename LabelType, LabelType R, LabelType N, typename IndexType=LabelType, typename FlatIndexType=LabelType, typename SignType=int>
constexpr SignType sort(std::array<LabelType, R>& Rtuple) {
	using IndexLookupType = unordered_to_index_of_ordered_TABLE<LabelType,R,N,IndexType,FlatIndexType>;
	using SignLookupType = unordered_to_sign_of_sorting_TABLE<LabelType,R,N,FlatIndexType,SignType>;
	FlatIndexType idx = IndexLookupType::flattened_index(Rtuple);
	Rtuple =  NchooseK<LabelType,N,R,IndexType>::array[
		IndexLookupType::flattened_table[idx]
	];
	return SignLookupType::flattened_table[idx];
	//return legacy::sort_array<LabelType,R,N,SignType>(Rtuple);
}

// This is a helper struct which enables the specification of all
// template parameters once and for all. 
//
// The entry `...::LIST` is a subtype which contains a list 
// (`...::LIST::array`) of all ordered `R`-tuples of distinct
// elements in `0..N-1`, representing `R`-element subsets thereof.
// `...::LIST` contains additional helper objects which assign a
// value to every such `R`-tuple. To manipulate `R`-tuples, like
// sort them or find their index in `...::LIST::array`, use the other
// functions inside this struct.
template<typename LabelType, LabelType R, LabelType N, typename IndexType=LabelType, typename FlatIndexType=IndexType, typename SignType=int>
struct RTUPLES {
	// Contains a list of all `R`-tuples (`...::array`), and auxulary
	// data to help query them.
	using LIST=NchooseK<LabelType,N,R,IndexType>;

	// Number of `R`-element subsets of `0..N-1`.
	static constexpr int NR = LIST::NR_RTUPLES;

	// Wrapper for `Rtuples::index_of_ordered`.
	constexpr static IndexType index_of_ordered(const std::array<LabelType, R>& Rtuple) {
		return Rtuples::index_of_ordered<LabelType,R,N,IndexType,FlatIndexType>(Rtuple);
	}
	// Wrapper for `Rtuples::index_of_unordered`.
	constexpr static IndexType index_of_unordered(const std::array<LabelType, R>& Rtuple) {
		return Rtuples::index_of_unordered<LabelType,R,N,IndexType,FlatIndexType>(Rtuple);
	}
	// Wrapper for `Rtuples::sign_of_sorting`.
	constexpr static SignType sign_of_sorting(const std::array<LabelType, R>& Rtuple) {
		return Rtuples::sign_of_sorting<LabelType,R,N,FlatIndexType,SignType>(Rtuple);
	}
	// Wrapper for `Rtuples::sign_and_index_of_unordered`
	constexpr static std::pair<SignType,IndexType> sign_and_index_of_unordered(const std::array<LabelType, R>& Rtuple) {
		return Rtuples::sign_and_index_of_unordered<LabelType,R,N,IndexType,FlatIndexType,SignType>(Rtuple);
	}
	// Wrapper for `Rtuples::sorted`.
	constexpr static std::array<LabelType, R> sorted(const std::array<LabelType, R>& Rtuple) {
		return Rtuples::sorted<LabelType,R,N,IndexType,FlatIndexType>(Rtuple);
	}
	// Wrapper for `Rtuples::sort`.
	constexpr static SignType sort(std::array<LabelType, R>& Rtuple) {
		return Rtuples::sort<LabelType,R,N,IndexType,FlatIndexType,SignType>(Rtuple);
	}
};
}