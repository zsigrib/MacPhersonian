#include <algorithm>
#include "OMoperations.hpp"

namespace OM_operations {

// =========
// UTILITIES
// =========

template<typename NumberType, NumberType N>
std::array<NumberType, N> monotone_function_skipping_number(NumberType number) {
    std::array<NumberType, N> func{};
    number = (number <= N) ? number : N;
    for (NumberType idx = 0; idx < number; ++idx) {
        func[idx] = idx;
    }
    for (NumberType idx = number; idx < N; ++idx) {
        func[idx] = idx + 1;
    }
    return func;
}

template<typename NumberType, NumberType N, typename Iterable>
std::array<NumberType, N> monotone_function_skipping_numbers(const Iterable& numbers) {
    std::vector<NumberType> numbers_sorted{};
    for (NumberType number: numbers) {
        numbers_sorted.push_back(number);
    }
    std::sort(numbers_sorted.begin(),numbers_sorted.end());
    std::array<NumberType, N> func{};
    NumberType idx = 0;
    NumberType next_to_write = 0;
    for (NumberType to_skip: numbers_sorted) {
        for (NumberType value = next_to_write; (value < to_skip) && (idx < N); ++value) {
            func[idx] = value;
            ++idx;
        }
        next_to_write = to_skip + 1;
    }
    for (; idx < N; ++idx) {
        func[idx] = next_to_write;
        ++next_to_write;
    }
    return func;
}

template<typename SmallLabelType, typename BigLabelType, 
SmallLabelType R, SmallLabelType N0, BigLabelType N1, typename Iterable, 
typename IndexType, typename FlatIndexType>
constexpr std::array<IndexType, binomial_coefficient(N0, R)>
Rtuple_function_induced_by_relabeling(const Iterable& f) {
    using RTUPLES0 = Rtuples::RTUPLES<SmallLabelType,R,N0,IndexType,FlatIndexType>;
    using RTUPLES1 = Rtuples::RTUPLES<BigLabelType,R,N1,IndexType,FlatIndexType>;
    std::array<IndexType, RTUPLES0::NR> induced{};
    std::array<BigLabelType, N0> function{};
    SmallLabelType x = 0;
    for (BigLabelType entry: f) {
        function[x] = entry;
        ++x;
    }
    for (IndexType idx = 0; idx < RTUPLES0::NR; ++idx) {
        std::array<BigLabelType, R> mapped_Rtuple{};
        for (BigLabelType t = 0; t < R; ++t) {
            mapped_Rtuple[t] = function[RTUPLES0::LIST::array[idx][t]];
        }
        induced[idx] = RTUPLES1::index_of_unordered(mapped_Rtuple);
    }
    return induced;
}

template<typename SmallLabelType, typename BigLabelType, 
SmallLabelType R, SmallLabelType N0, BigLabelType N1, typename Iterable, 
typename IndexType, typename FlatIndexType>
constexpr std::array<IndexType, binomial_coefficient(N0, R)>
Rtuple_function_induced_by_monotone_relabeling(const Iterable& f) {
    using RTUPLES0 = Rtuples::RTUPLES<SmallLabelType,R,N0,IndexType,FlatIndexType>;
    using RTUPLES1 = Rtuples::RTUPLES<BigLabelType,R,N1,IndexType,FlatIndexType>;
    std::array<IndexType, RTUPLES0::NR> induced{};
    std::array<BigLabelType, N0> function{};
    SmallLabelType x = 0;
    for (BigLabelType entry: f) {
        function[x] = entry;
        ++x;
    }
    for (IndexType idx = 0; idx < RTUPLES0::NR; ++idx) {
        std::array<BigLabelType, R> mapped_Rtuple{};
        for (BigLabelType t = 0; t < R; ++t) {
            mapped_Rtuple[t] = function[RTUPLES0::LIST::array[idx][t]];
        }
        induced[idx] = RTUPLES1::index_of_ordered(mapped_Rtuple);
    }
    return induced;
}

// =============
// REORIENTATION
// =============

template<int R, int N>
constexpr svo::Multiply<Chirotope<R, N>> reorient_element(int element) {
    return {{Chirotope<R,N>::RTUPLES::LIST::contained_mask32[element]}};
}

template<int R, int N, typename Iterable>
constexpr svo::Multiply<Chirotope<R, N>> reorient_elements(const Iterable& elements) {
    using RTUPLES = Chirotope<R,N>::RTUPLES;
    bit_vector<RTUPLES::NR> to_reorient{};
    for (auto element: elements) {
        for (auto i = 0; i < bit_vector<RTUPLES::NR>::NR_INT32; ++i) {
            to_reorient[i] ^= RTUPLES::LIST::contained_mask32[element][i];
        }
    }
    return {to_reorient};
}

// ==========
// RELABELING
// ==========

template<int R, int N0, int N1, typename Iterable>
constexpr svo::composes_to<svo::PushForward<Chirotope<R, N0>, Chirotope<R, N1>>,svo::Multiply<Chirotope<R,N0>>>
relabel_elements(const Iterable& f) {
    using RTUPLES0 = Chirotope<R, N0>::RTUPLES;
    using RTUPLES1 = Chirotope<R, N1>::RTUPLES;
    std::array<int, RTUPLES0::NR> induced{};
    bit_vector<RTUPLES0::NR> to_reorient{};
    std::array<char, N0> function{};
    char x = 0;
    for (char entry: f) {
        function[x] = entry;
        ++x;
    }
    for (int idx = 0; idx < RTUPLES0::NR; ++idx) {
        std::array<char, R> mapped_Rtuple{};
        for (char t = 0; t < R; ++t) {
            mapped_Rtuple[t] = function[RTUPLES0::LIST::array[idx][t]];
        }
        auto sgn_and_idx = RTUPLES1::sign_and_index_of_unordered(mapped_Rtuple);
        induced[idx] = sgn_and_idx.second;
        to_reorient.set_bit(idx, sgn_and_idx.first == -1);
    }
    return {{{induced}},{to_reorient}};
}

template<int R, int N0, int N1, typename Iterable>
constexpr svo::PushForward<Chirotope<R, N0>, Chirotope<R, N1>>
relabel_monotonely(const Iterable& f) {
    return {Rtuple_function_induced_by_monotone_relabeling<char,char,R,N0,N1,Iterable,int>(f)};
}

template<int R0, int N0, int R1, int N1, typename Iterable>
constexpr svo::composes_to<svo::Multiply<Chirotope<R1,N1>>,svo::PullBack<Chirotope<R0,N0>,Chirotope<R1,N1>>>
general_minor(const Iterable& f) {
    using RTUPLES0 = Chirotope<R0, N0>::RTUPLES;
    using RTUPLES1 = Chirotope<R1, N1>::RTUPLES;
    using RTUPLES01 = Chirotope<R1, N0>::RTUPLES;
    std::array<char, N1> function{};
    std::array<bool, N0> in_image{};
    std::array<char, R0 - R1> not_in_image{};
    char x = 0;
    for (char entry: f) {
        function[x] = entry;
        in_image[entry] = true;
        ++x;
    }
    x = 0;
    for (char t = 0; t < N0; ++t) {
        if (!in_image[t]) {
            not_in_image[x] = t;
            ++x;
        }
    }
    std::array<int, RTUPLES1::NR> induced{};
    bit_vector<RTUPLES1::NR> to_reorient{};
    for (int idx = 0; idx < RTUPLES1::NR; ++idx) {
        std::array<char, R0> mapped_Rtuple_with_new_coloops{};
        std::array<char, R1> mapped_Rtuple_without_new_coloops{};
        for (char t = 0; t < R1; ++t) {
            mapped_Rtuple_with_new_coloops[t] = function[RTUPLES1::LIST::array[idx][t]];
            mapped_Rtuple_without_new_coloops[t] = function[RTUPLES1::LIST::array[idx][t]];
        }
        for (char t = R1; t < R0; ++t) {
            mapped_Rtuple_with_new_coloops[t] = not_in_image[t - R1];
        }
        induced[idx] = RTUPLES0::index_of_unordered(mapped_Rtuple_with_new_coloops);
        to_reorient.set_bit(idx, -1 == RTUPLES01::sign_of_sorting(mapped_Rtuple_without_new_coloops));
    }
    return {{to_reorient},{{induced}}};
}

template<int R0, int N0, int R1, int N1, typename Iterable>
constexpr svo::PullBack<Chirotope<R0, N0>, Chirotope<R1, N1>>
monotone_minor(const Iterable& f) {
    using RTUPLES0 = Chirotope<R0, N0>::RTUPLES;
    using RTUPLES1 = Chirotope<R1, N1>::RTUPLES;
    using RTUPLES01 = Chirotope<R1, N0>::RTUPLES;
    std::array<char, N1> function{};
    std::array<bool, N0> in_image{};
    std::array<char, R0 - R1> not_in_image{};
    char x = 0;
    for (char entry: f) {
        function[x] = entry;
        in_image[entry] = true;
        ++x;
    }
    x = 0;
    for (char t = 0; t < N0; ++t) {
        if (!in_image[t]) {
            not_in_image[x] = t;
            ++x;
        }
    }
    std::array<int, RTUPLES1::NR> induced{};
    for (int idx = 0; idx < RTUPLES1::NR; ++idx) {
        std::array<char, R0> mapped_Rtuple_with_new_coloops{};
        for (char t = 0; t < R1; ++t) {
            mapped_Rtuple_with_new_coloops[t] = function[RTUPLES1::LIST::array[idx][t]];
        }
        for (char t = R1; t < R0; ++t) {
            mapped_Rtuple_with_new_coloops[t] = not_in_image[t - R1];
        }
        induced[idx] = RTUPLES0::index_of_ordered(mapped_Rtuple_with_new_coloops);
    }
    return {{induced}};
}

template<int R, int N>
constexpr svo::PushForward<Chirotope<R,N>,Chirotope<R,N+1>>
insert_loop_at(int element) {
    return relabel_monotonely<R,N,N+1>(monotone_function_skipping_number<char,N>(element));
}

template<int R, int N>
constexpr svo::PushForward<Chirotope<R,N>,Chirotope<R+1,N+1>>
insert_coloop_at(int element) {
    using RTUPLES0 = Chirotope<R, N>::RTUPLES;
    using RTUPLES1 = Chirotope<R+1, N+1>::RTUPLES;
    std::array<int, RTUPLES0::NR> induced{};
    for (int idx = 0; idx < RTUPLES0::NR; ++idx) {
        std::array<char, R+1> mapped_Rtuple_with_new_coloops{};
        char t = 0;
        while (RTUPLES0::LIST::array[idx][t] < element && t < R) {
            mapped_Rtuple_with_new_coloops[t] = RTUPLES0::LIST::array[idx][t];
            ++t;
        }
        mapped_Rtuple_with_new_coloops[t] = element;
        while (t < R) {
            mapped_Rtuple_with_new_coloops[t + 1] = RTUPLES0::LIST::array[idx][t] + 1;
            ++t;
        }
        induced[idx] = RTUPLES1::index_of_ordered(mapped_Rtuple_with_new_coloops);
    }
    return {{induced}};
}

// =====================
// DELETION, CONTRACTION
// =====================

template<int R, int N>
constexpr svo::Multiply_P0<Chirotope<R, N>> delete_element(int element) {
    return {~bit_vector<Chirotope<R,N>::RTUPLES::NR>{
        Chirotope<R,N>::RTUPLES::LIST::contained_mask32[element]
    }};
}

template<int R, int N>
constexpr svo::PullBack<Chirotope<R,N>,Chirotope<R,N-1>> delete_label(int label) {
    return {{
        Rtuple_function_induced_by_monotone_relabeling
        <char,char,R,N-1,N,std::array<char,N-1>,int>(
            monotone_function_skipping_number<char,N-1>(label)
        )
    }};
}

template<int R, int N, typename Iterable>
constexpr svo::Multiply_P0<Chirotope<R, N>> delete_elements(const Iterable& elements) {
    using RTUPLES = Chirotope<R,N>::RTUPLES;
    bit_vector<RTUPLES::NR> to_delete{};
    for (int element: elements) {
        for (auto i = 0; i < bit_vector<RTUPLES::NR>::NR_INT32; ++i) {
            to_delete[i] |= RTUPLES::LIST::contained_mask32[element][i];
        }
    }
    return {~to_delete};
}

template<int R, int N0, int N1, typename Iterable>
constexpr svo::PullBack<Chirotope<R,N0>,Chirotope<R,N1>> delete_labels(const Iterable& labels) {
    return {{
        Rtuple_function_induced_by_monotone_relabeling
        <char,char,R,N1,N0,std::array<char, N1>,int>(
            monotone_function_skipping_numbers<char,N1>(labels)
        )
    }};
}

template<int R, int N>
constexpr svo::Multiply_P0<Chirotope<R, N>> contract_element(int element) {
    using RTUPLES = Chirotope<R, N>::RTUPLES;
    return {bit_vector<RTUPLES::NR>(RTUPLES::LIST::contained_mask32[element])};
}

template<int R, int N>
constexpr svo::PullBack<Chirotope<R,N>,Chirotope<R-1,N-1>>
contract_label(int label) {
    return monotone_minor<R,N,R-1,N-1>(monotone_function_skipping_number<char,N-1>(label));
}

template<int R, int N, typename Iterable>
constexpr svo::Multiply_P0<Chirotope<R ,N>> contract_elements(const Iterable& elements) {
    using RTUPLES = Chirotope<R, N>::RTUPLES;
    auto to_keep = ~bit_vector<RTUPLES::NR>();
    for (auto element: elements) {
        for (auto i = 0; i < bit_vector<RTUPLES::NR>::NR_INT32; ++i) {
            to_keep[i] &= RTUPLES::LIST::contained_mask32[element][i];
        }
    }
    return {to_keep};
}

template<int R0, int N0, int R1, int N1, typename Iterable>
constexpr svo::PullBack<Chirotope<R0 ,N0>, Chirotope<R1, N1>> 
contract_labels(const Iterable& labels) {
    return monotone_minor<R0,N0,R1,N1>(monotone_function_skipping_numbers<char,N1>(labels));
}

template<int R, int N>
constexpr std::vector<sign_vector<N>> cocircuits(const Chirotope<R, N>& chi) {
    using R1TUPLES = Rtuples::RTUPLES<char,R-1,N,int>; // (R-1)-tuples
    constexpr auto N_NR_INT32 = division_rounded_up(N, 32);
    std::vector<sign_vector<N>> cocircuits_so_far{};
    auto candidate = cocircuit_extractors_from_chirotope<R,N>[0] * chi;
    if (!candidate.is_zero()) {
        cocircuits_so_far.push_back(candidate);
    }
    for (int idx = 1; idx < R1TUPLES::NR; ++idx) {
        auto candidate = cocircuit_extractors_from_chirotope<R,N>[idx] * chi;
        bool seen_before = false;
        for (auto i = 0; i < (R1TUPLES::LIST::array[idx][R-2] >> 5); ++i) {
            if (~(candidate.plus[i] | candidate.minus[i] 
            | R1TUPLES::LIST::char_vector_of_Rtuple[idx][i])) {
                seen_before = true;
                break;
            }
        }
        if (N % 32 != 0 && ~(candidate.plus[N_NR_INT32 - 1] 
        | candidate.minus[N_NR_INT32 - 1] 
        | R1TUPLES::LIST::char_vector_of_Rtuple[idx][N_NR_INT32 - 1]
        | (-((uint32_t)1 << ((R1TUPLES::LIST::array[idx][R-2] + 1) & 31))))) {
            seen_before = true;
        }
        if (!seen_before) {
            cocircuits_so_far.push_back(candidate);
        }
    }
    return cocircuits_so_far;
}

}