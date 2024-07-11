#include <type_traits>
#include <iostream>
#include <format>
#include "mymath.hpp"
#include "OMs.hpp"

// ===============
// Chirotope<R, N>
// ===============

template<int R, int N>
std::ostream& operator<<(std::ostream& os, const Chirotope<R, N>& chi) {
    for (auto i = 0; i < chi.NR_INT32 - 1; i++) {
        for (auto c = 0; c < 32; c++) {
            if (chi.plus[i] & ((uint32_t)1 << c))
                os << '+';
            else if (chi.minus[i] & ((uint32_t)1 << c))
                os << '-';
            else
                os << '0';
        }
    }
    for (auto c = 0; c < chi.NR_REMAINING_BITS; c++) {
        if (chi.plus[chi.NR_INT32 - 1] & ((uint32_t)1 << c))
            os << '+';
        else if (chi.minus[chi.NR_INT32 - 1] & ((uint32_t)1 << c))
            os << '-';
        else
            os << '0';
    }
    return os;
}

template<int R, int N>
std::ofstream& operator<<(std::ofstream& of, const Chirotope<R, N>& chi) {
    for (auto i = 0; i < chi.NR_INT32 - 1; i++) {
        for (auto c = 0; c < 32; c++) {
            if (chi.plus[i] & ((uint32_t)1 << c))
                of << '+';
            else if (chi.minus[i] & ((uint32_t)1 << c))
                of << '-';
            else
                of << '0';
        }
    }
    for (auto c = 0; c < chi.NR_REMAINING_BITS; c++) {
        if (chi.plus[chi.NR_INT32 - 1] & ((uint32_t)1 << c))
            of << '+';
        else if (chi.minus[chi.NR_INT32 - 1] & ((uint32_t)1 << c))
            of << '-';
        else
            of << '0';
    }
    return of;
}

// ======================
// ReadChirotopesFromFile
// ======================

std::ifstream& ignore_n_lines(std::ifstream& fs, size_t ignored_lines) {
    std::string ignore_this_str;
    for (auto i = 0; i < ignored_lines; i++) {
        if (!std::getline(fs, ignore_this_str)) {
            throw std::invalid_argument(std::format(
                "The file stream does not contain at least {0} lines, which we expect to ignore.",
                ignored_lines
            ));
        }
    }
    return fs;
}

template<int R, int N>
ReadChirotopesFromFile<R,N>::ReadChirotopesFromFile
(std::string path_, int ignored_lines, bool exhausted_): 
path(path_), file_stream(path_.c_str()), 
file_found(file_stream.good()), ignore_first_n_lines(ignored_lines) {
    if (exhausted_) {
        exhausted = true;
    } else {
        exhausted = !file_found;
        if (file_found) {
            ignore_n_lines(file_stream, ignored_lines);
        }
        current_chirotope = load_next_chirotope();
    }  
}

template<int R, int N>
ReadChirotopesFromFile<R, N>::ReadChirotopesFromFile
(const ReadChirotopesFromFile<R, N>& other):
path(other.path), file_stream(other.path.c_str()),
file_found(file_stream.good()), ignore_first_n_lines(other.ignore_first_n_lines) {
    if (file_found != other.file_found) {
        throw std::invalid_argument("Could not copy construct ReadChirotopesFromFiles, as the target file was created/destroyed since the creation of the original instance.");
    } else if (other.is_exhausted()) {
        exhausted = true;
    } else {
        exhausted = false;
        ignore_n_lines(file_stream, ignore_first_n_lines);
        current_chirotope = load_next_chirotope();
        while (!exhausted && current_chirotope != other.current_chirotope) {
            current_chirotope = load_next_chirotope();
        }
        if (exhausted) {
            throw std::invalid_argument("Could not copy construct ReadChirotopesFromFiles, because could not find current Chirotope of the original instance in the file.");
        }
    }
}

template<int R, int N>
ReadChirotopesFromFile<R, N>& ReadChirotopesFromFile<R, N>::change_path(std::string new_path) {
    path = new_path;
    file_stream.close();
    file_stream.open(path.c_str());
    file_found = file_stream.good();
    exhausted = !file_found;
    if (file_found) {
        ignore_n_lines(file_stream, ignore_first_n_lines);
    }
    current_chirotope = load_next_chirotope();
    return (*this);
}

template<int R, int N>
Chirotope<R, N> ReadChirotopesFromFile<R, N>::load_next_chirotope() {
    if (exhausted) return Chirotope<R, N>();
    std::string next_line;
    if (!std::getline(file_stream, next_line)) {
        exhausted = true;
        return Chirotope<R, N>();
    } else if (next_line.empty()) {
        return load_next_chirotope();
    } else {
        return Chirotope<R,N>().read(next_line);
    }
}

template<int R, int N>
constexpr int ReadChirotopesFromFile<R, N>::ignored_number_of_lines() const {
    return ignore_first_n_lines;
}

template<int R, int N>
bool ReadChirotopesFromFile<R, N>::file_exists() const {
    return file_found;
}

template<int R, int N>
bool ReadChirotopesFromFile<R, N>::is_exhausted() const {
    return exhausted;
}

template<int R, int N>
Chirotope<R, N> ReadChirotopesFromFile<R, N>::operator*() const {
    return current_chirotope;
}

template<int R, int N>
ReadChirotopesFromFile<R, N>& ReadChirotopesFromFile<R, N>::operator++() {
    current_chirotope = load_next_chirotope();
    return (*this);
}

template<int R, int N>
bool ReadChirotopesFromFile<R, N>::operator==(const ReadChirotopesFromFile<R, N>& other) const {
    return (exhausted && other.exhausted
        && (path == other.path) 
    ) || (!exhausted && !other.exhausted
        && (path == other.path)
        && (current_chirotope == other.current_chirotope)
    );
}

template<int R, int N>
bool ReadChirotopesFromFile<R, N>::operator!=(const ReadChirotopesFromFile<R, N>& other) const {
    return (!exhausted || !other.exhausted
        || (path != other.path) 
    ) && (exhausted || other.exhausted
        || (path != other.path)
        || (current_chirotope != other.current_chirotope)
    );
}

template<int R, int N>
ReadChirotopesFromFile<R, N>& ReadChirotopesFromFile<R, N>::begin() {
    return (*this);
}

template<int R, int N>
ReadChirotopesFromFile<R, N> ReadChirotopesFromFile<R, N>::end() const {
    return ReadChirotopesFromFile<R,N>(path, ignore_first_n_lines, true);
}

// =======================
// ReadChirotopesFromFiles
// =======================

template<int R, int N>
ReadChirotopesFromFiles<R, N>::ReadChirotopesFromFiles(
    std::string(*file_path_constructor)(int,int), 
    int ignored_number_of_lines,
    int number_of_bases_, 
    int idx_of_file_
):  number_of_bases(number_of_bases_),
    idx_of_file(idx_of_file_),
    path_constructor(file_path_constructor),
    file_reader(file_path_constructor(number_of_bases_, idx_of_file_), ignored_number_of_lines)
{
    update_file_reader();
}

template<int R, int N>
ReadChirotopesFromFiles<R, N>::ReadChirotopesFromFiles
(const ReadChirotopesFromFiles& other):
file_reader(other.path_constructor(other.number_of_bases, other.idx_of_file), other.ignored_number_of_lines()) {
    path_constructor = other.path_constructor;
    if (other.is_exhausted()) {
        number_of_bases = binomial_coefficient(N, R) + 1;
        idx_of_file = 0;
        update_file_reader();
    } else {
        number_of_bases = other.number_of_bases;
        idx_of_file = other.idx_of_file;
        update_file_reader();
    }
}

template<int R, int N>
ReadChirotopesFromFiles<R, N>::ReadChirotopesFromFiles
(std::string(*file_path_constructor)(int, int), int ignored_number_of_lines): 
file_reader(file_path_constructor(1, 0), ignored_number_of_lines) {
    path_constructor = file_path_constructor;
    number_of_bases = 1;
    idx_of_file = 0;
    update_file_reader();
}

template<int R, int N>
ReadChirotopesFromFile<R, N>& ReadChirotopesFromFiles<R, N>::update_file_reader() {
    file_reader.change_path(path_constructor(number_of_bases, idx_of_file));
    if (file_reader.file_exists() || idx_of_file == 0) return file_reader;
    idx_of_file = 0;
    number_of_bases++;
    file_reader.change_path(path_constructor(number_of_bases, idx_of_file));
    while (file_reader.file_exists() && file_reader.is_exhausted()) {
        number_of_bases++;
        file_reader.change_path(path_constructor(number_of_bases, idx_of_file));
    }
    return file_reader;
}

template<int R, int N>
constexpr bool ReadChirotopesFromFiles<R, N>::is_exhausted() const {
    return !file_reader.file_exists();
}

template<int R, int N>
constexpr int ReadChirotopesFromFiles<R, N>::ignored_number_of_lines() const {
    return file_reader.ignored_number_of_lines();
}

template<int R, int N>
std::pair<int, Chirotope<R, N>> ReadChirotopesFromFiles<R, N>::operator*() const {
    return std::pair<int, Chirotope<R, N>>(number_of_bases, *file_reader);
}

template<int R, int N>
ReadChirotopesFromFiles<R, N>& ReadChirotopesFromFiles<R, N>::operator++() {
    ++file_reader;
    if (file_reader.is_exhausted()) {
        idx_of_file++;
        update_file_reader();
    } 
    return (*this);
}

template<int R, int N>
bool ReadChirotopesFromFiles<R, N>::operator==(const ReadChirotopesFromFiles& other) const {
    return (!is_exhausted() && !other.is_exhausted()
        && (number_of_bases == other.number_of_bases)
        && (idx_of_file == other.idx_of_file)
        && (file_reader == other.file_reader)
        && (path_constructor) == (other.path_constructor)
    ) || (is_exhausted() && other.is_exhausted()
        && path_constructor == other.path_constructor
    );
}

template<int R, int N>
bool ReadChirotopesFromFiles<R, N>::operator!=(const ReadChirotopesFromFiles& other) const {
    return (is_exhausted() || other.is_exhausted()
        || (number_of_bases != other.number_of_bases)
        || (idx_of_file != other.idx_of_file)
        || (file_reader != other.file_reader)
        || (path_constructor) != (other.path_constructor)
    ) && (!is_exhausted() || !other.is_exhausted()
        || path_constructor != other.path_constructor
    );
}

template<int R, int N>
ReadChirotopesFromFiles<R, N>& ReadChirotopesFromFiles<R, N>::begin() {
    return (*this);
}

template<int R, int N>
ReadChirotopesFromFiles<R, N> ReadChirotopesFromFiles<R, N>::end() const {
    return ReadChirotopesFromFiles<R, N>(
        path_constructor,
        file_reader.ignored_number_of_lines(),
        binomial_coefficient(N, R) + 1,
        0
    );
}