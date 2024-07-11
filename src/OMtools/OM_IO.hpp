#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include "OMs.hpp"

// ============================
// ReadChirotopesFromFile<R, N>
// ============================

// Reads all chirotopes in a given file, assuming each chirotope
// occupies a full line, ignoring empty lines and the first few
// lines as specified by the user. Usage:
// ```
// for (Chirotope<R, N> p : ReadChirotopesFromFile<R,N>("path/to/file", ig)) {
//     /* USER CODE*/
// }
// ```
// where the first `ig` lines of the given file will be ignored.
template<int R, int N>
struct ReadChirotopesFromFile {
    using value_type = Chirotope<R, N>;
    
    private:
    std::string path;
    const int ignore_first_n_lines;
    std::ifstream file_stream;
    bool file_found;
    bool exhausted;
    Chirotope<R, N> current_chirotope;

    public:
    ReadChirotopesFromFile(std::string path, int ignored_lines, bool exhausted_=false);
    ReadChirotopesFromFile(const ReadChirotopesFromFile&);

    ReadChirotopesFromFile& change_path(std::string);
    Chirotope<R, N> load_next_chirotope();
    constexpr int ignored_number_of_lines() const;
    bool file_exists() const;
    bool is_exhausted() const;
    Chirotope<R,N> operator*() const;
    ReadChirotopesFromFile& operator++();
    bool operator==(const ReadChirotopesFromFile&) const;
    bool operator!=(const ReadChirotopesFromFile&) const;
    ReadChirotopesFromFile& begin();
    ReadChirotopesFromFile end() const;
};

// =============================
// ReadChirotopesFromFiles<R, N>
// =============================

// Reads all chirotopes in a given collection of files.
//
// The files must be as follows. All chirotopes in a given
// file must have the same number of bases and must occupy
// an entire row. The first few lines of a file may be junk,
// but this have to be specified (see later). Empty lines
// are permitted and ignored. A file is allowed to contain
// 0 chirotopes. There may be multiple files which contain
// chirotopes with the same number of bases.
//
// To construct an instance of `ReadChirotopesFromFiles`,
// a function pointer `std::str (*)(int n_bases, int idx)`
// must be specified, which constructs the path of the
// `idx`th file which contains chirotopes with `n_bases`
// bases. `idx` starts at 0, and if the file at the path
// generated with parameters `(n_bases, idx)` is missing,
// then all files with parameters `(n_bases, i)` will be
// ignored for `i > id`.
//
// An `int ignore_first_n_lines` must be specified: the
// first this many lines of the file are thrown out before
// looking for chirotopes.
//
// Usage:
// ```
// for (std::pair<int, Chirotope<R, N>> p : ReadChirotopesFromFiles<R,N>
// (&path_generating_function, ig)) {
//     /* USER CODE*/
// }
// ```
// where :
// - `path_generating_function` is the function mentioned above, of
//   of signature `std::string path_generating_function(int, int)`,
// - `ig` is an integer, the first this many lines will be ignored in
//   in each file,
// - `p.first` is the number of bases of `p.second`; it is non-decreasing,
// - `p.second` is the next chirotope read from the file structure.
template<int R, int N>
struct ReadChirotopesFromFiles {
    using value_type = Chirotope<R, N>;

    private:
    std::string (*path_constructor)(int, int);
    int number_of_bases;
    int idx_of_file;
    ReadChirotopesFromFile<R, N> file_reader;

    ReadChirotopesFromFiles(std::string(*)(int,int), int, int, int);

    public:
    ReadChirotopesFromFiles(const ReadChirotopesFromFiles&);

    ReadChirotopesFromFiles(std::string(*)(int, int), int);
    
    ReadChirotopesFromFile<R, N>& update_file_reader();
    std::ifstream& get_file_stream();
    constexpr bool is_exhausted() const;
    constexpr int ignored_number_of_lines() const;
    Chirotope<R, N> load_next_chirotope();
    
    std::pair<int, Chirotope<R, N>> operator*() const;
    ReadChirotopesFromFiles& operator++();
    bool operator==(const ReadChirotopesFromFiles&) const;
    bool operator!=(const ReadChirotopesFromFiles&) const;

    ReadChirotopesFromFiles& begin();
    ReadChirotopesFromFiles end() const;
};

#include "OM_IO.cpp"