#include <iostream>
#include <string>
#include "OMs.hpp"

// =======
// HELPERS 
// =======

constexpr int count_1bits(uint32_t n)
{
    n = (((uint32_t)0xaaaaaaaa & n) >> 1) + ((uint32_t)0x55555555 & n);
    n = (((uint32_t)0xcccccccc & n) >> 2) + ((uint32_t)0x33333333 & n);
    n = (((uint32_t)0xf0f0f0f0 & n) >> 4) + ((uint32_t)0x0f0f0f0f & n);
    n = (((uint32_t)0xff00ff00 & n) >> 8) + ((uint32_t)0x00ff00ff & n);
    n = (((uint32_t)0xffff0000 & n) >> 16) + ((uint32_t)0x0000ffff & n);
    return n;
}

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
constexpr int sort_array(std::array<char, R>& a) {
    char p,q;
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
	


	int sorted,min,ind_min,i,j,sign;				//slow code, but works for all R

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

template<int R, int N>
constexpr int index_of_Rtuple(const std::array<char, R>& a) {
	if (R==3 && N<=9)		//quick code for common used parameters
	{
		if (N==9)
		{	
			if (a[0]==0)
				if (a[1]==1)
					return a[2]-2;
				else if (a[1]==2)
					return a[2]+4;
				else if (a[1]==3)
					return a[2]+9;
				else if (a[1]==4)
					return a[2]+13;
				else if (a[1]==5)
					return a[2]+16;			
				else	return a[1]+a[2]+12;
			if (a[0]==1)
				if (a[1]==2)
					return a[2]+25;
				else if (a[1]==3)
					return a[2]+30;
				else if (a[1]==4)
					return a[2]+34;
				else if (a[1]==5)
					return a[2]+37;
				else    return a[1]+a[2]+33;
			if (a[0]==2)
				if (a[1]==3)
					return a[2]+45;
				else if (a[1]==4)
					return a[2]+49;
				else if (a[1]==5)
					return a[2]+52;
				else    return a[1]+a[2]+48;
			if (a[0]==3)
				if (a[1]==4)
					return a[2]+59;
				else if (a[1]==5)
					return a[2]+62;
				else    return a[1]+a[2]+58;
			if (a[0]==4)
				if (a[1]==5)
					return a[2]+68;
				else return a[1]+a[2]+64;
			return a[0]+a[1]+a[2]+62;
		}
		if (N==8)
		{	
			if (a[0]==0)
				if (a[1]==1)
					return a[2]-2;
				else if (a[1]==2)
					return a[2]+3;
				else if (a[1]==3)
					return a[2]+7;
				else if (a[1]==4)
					return a[2]+10;			
				else	return a[1]+a[2]+7;
			if (a[0]==1)
				if (a[1]==2)
					return a[2]+18;
				else if (a[1]==3)
					return a[2]+22;
				else if (a[1]==4)
					return a[2]+25;
				else    return a[1]+a[2]+22;
			if (a[0]==2)
				if (a[1]==3)
					return a[2]+32;
				else if (a[1]==4)
					return a[2]+35;
				else    return a[1]+a[2]+32;
			if (a[0]==3)
				if (a[1]==4)
					return a[2]+41;
				else return a[1]+a[2]+38; 
			return a[0]+a[1]+a[2]+37;
		}
		if (N==7)
		{
			if (a[0]==0)
				if (a[1]==1)
					return a[2]-2;
				else if (a[1]==2)
					return a[2]+2;
				else if (a[1]==3)
					return a[2]+5;			
				else	return a[1]+a[2]+3;
			if (a[0]==1)
				if (a[1]==2)
					return a[2]+12;
				else if (a[1]==3)
					return a[2]+15;
				else    return a[1]+a[2]+13;
			if (a[0]==2)
				if (a[1]==3)
					return a[2]+21;
				else    return a[1]+a[2]+19; 
			return a[0]+a[1]+a[2]+19;
		}
		
		if (N==6)
		{
			if (a[0]==0)
				if (a[1]==1)
					return a[2]-2;
				else if (a[1]==2)
					return a[2]+1;
				else	return a[1]+a[2];
			if (a[0]==1)
				if (a[1]==2)
					return a[2]+7;
				else    return a[1]+a[2]+6; 
			return a[0]+a[1]+a[2]+7;
		
		}
		
		if (N==5)
		{
			if (a[0]==0)
				if (a[1]==1)
					return a[2]-2;
				else	return a[1]+a[2]-2;
			return a[0]+a[1]+a[2];		
		}
		if (N==4)
		{
			return a[0]+a[1]+a[2]-3;
		}
	}
	if (R==4 && N<=9)
	{
		if (N==9)
		{
			if (a[0]==0)
				if (a[1]==1)
					if (a[2]==2)
						return a[3]-3;
					else if (a[2]==3)
						return a[3]+2;
					else if (a[2]==4)
						return a[3]+6;
					else if (a[2]==5)
						return a[3]+9;
					else return a[2]+a[3]+5;
				else if (a[1]==2)
					if (a[2]==3)
						return a[3]+17;
					else if (a[2]==4)
						return a[3]+21;
					else if (a[2]==5)
						return a[3]+24;
					else return a[2]+a[3]+20;
				else if (a[1]==3)
					if (a[2]==4)
						return a[3]+31;
					else if (a[2]==5)
						return a[3]+34;
					else return a[2]+a[3]+30;
				else if (a[1]==4)
					if (a[2]==5)
						return a[3]+40;
					else return a[2]+a[3]+36;
				else return a[1]+a[2]+a[3]+34;
			if (a[0]==1)
				if (a[1]==2)
					if (a[2]==3)
						return a[3]+52;
					else if (a[2]==4)
						return a[3]+56;
					else if (a[2]==5)
						return a[3]+59;
					else return a[2]+a[3]+55;
				else if (a[1]==3)
					if (a[2]==4)
						return a[3]+66;
					else if (a[2]==5)
						return a[3]+69;
					else return a[2]+a[3]+65;
				else if (a[1]==4)
					if (a[2]==5)
						return a[3]+75;
					else return a[2]+a[3]+71;
				else return a[1]+a[2]+a[3]+69;
			if (a[0]==2)
				if (a[1]==3)
					if (a[2]==4)
						return a[3]+86;
					else if (a[2]==5)
						return a[3]+89;
					else return a[2]+a[3]+85;
				else if (a[1]==4)
					if (a[2]==5)
						return a[3]+95;
					else return a[2]+a[3]+91;
				else return a[1]+a[2]+a[3]+89;
			if (a[0]==3)
				if (a[1]==4)
					if (a[2]==5)
						return a[3]+105;
					else return a[2]+a[3]+101;
				else return a[1]+a[2]+a[3]+99;
			return a[0]+a[1]+a[2]+a[3]+99;
		}

		if (N==8)
		{
			if (a[0]==0)
				if (a[1]==1)
					if (a[2]==2)
						return a[3]-3;
					else if (a[2]==3)
						return a[3]+1;
					else if (a[2]==4)
						return a[3]+4;
					else return a[2]+a[3]+1;
				else if (a[1]==2)
					if (a[2]==3)
						return a[3]+11;
					else if (a[2]==4)
						return a[3]+14;
					else return a[2]+a[3]+11;
				else if (a[1]==3)
					if (a[2]==4)
						return a[3]+20;
					else return a[2]+a[3]+17;
				else return a[1]+a[2]+a[3]+16;
			if (a[0]==1)
				if (a[1]==2)
					if (a[2]==3)
						return a[3]+31;
					else if (a[2]==4)
						return a[3]+34;
					else return a[2]+a[3]+31;
				else if (a[1]==3)
					if (a[2]==4)
						return a[3]+40;
					else return a[2]+a[3]+37;
				else return a[1]+a[2]+a[3]+36;
			if (a[0]==2)
				if (a[1]==3)
					if (a[2]==4)
						return a[3]+50;
					else return a[2]+a[3]+47;
				else return a[1]+a[2]+a[3]+46;
			return a[0]+a[1]+a[2]+a[3]+47;
		}

		if (N==7)
		{
			if (a[0]==0)
				if (a[1]==1)
					if (a[2]==2)
						return a[3]-3;
					else if (a[2]==3)
						return a[3];
					else return a[2]+a[3]-2;
				else if (a[1]==2)
					if (a[2]==3)
						return a[3]+6;
					else return a[2]+a[3]+4;
				else return a[1]+a[2]+a[3]+4;
			if (a[0]==1)
				if (a[1]==2)
					if (a[2]==3)
						return a[3]+16;
					else return a[2]+a[3]+14;
				else return a[1]+a[2]+a[3]+14;
			return a[0]+a[1]+a[2]+a[3]+16;
		}
			
		if (N==6)
		{
			if (a[0]==0)
				if (a[1]==1)
					if (a[2]==2)
						return a[3]-3;
					else return a[2]+a[3]-4;
				else return a[1]+a[2]+a[3]-3;
			return a[0]+a[1]+a[2]+a[3];
		}			
		if (N==5)
			return a[0]+a[1]+a[2]+a[3]-6;
	}
	
	
	int l,u,m,i;						//code that works for every parameter, slow...
	l=0;
	u=Chirotope<R, N>::RTUPLES::NR-1;
	m=(l+u)>>1;
	const auto& bases = NchooseK<int, N, R>::array;

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
		while (u< Chirotope<R, N>::RTUPLES::NR && a[i]==bases[u][i])
			u++;
		u--;
		if (l==u)
			return l;
		m=(u+l)>>1;
	}
	throw std::invalid_argument("We shouldn't have reached this line in index_of_Rtuple.");
}

template<int R, int N>
constexpr bool Chirotope<R, N>::is_chirotope() const {
	char i,j,k,l,p,q;
	long long int h;

	for (k=0;k<BASE::NR_INT32;k++)
		if (BASE::plus[k] & BASE::minus[k]) 		//if the same basis is both positive and negative
			return 0;
	

	for (k=0;k<BASE::NR_INT32;k++)
		if (BASE::plus[k]!=0 || BASE::minus[k]!=0) 		//(B0)
			break;

	if (k==BASE::NR_INT32)
		return 0; 						


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
				for (j;j<limit_j;j++)
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
	if (X[0] == 1 && X[1] == 0 && X[2] == 2
	&& Y[0] == 0 && Y[1] == 2 && Y[2] == 3) {
		//std::cout << "!!!\n";
	}
		

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