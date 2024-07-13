#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include "OMs.hpp"

// ============================
// ReadOMDataFromFile<R, N>
// ============================

// Reads a collection of (whitespace separated) data entries from a file,
// assuming the typename `StreamReadable` is a type which can be read into 
// using the `>>` operator, and which can be default constructed.
// The first few lines are skipped, as specified by the user. Usage:
// ```
// for (StreamReadable p : ReadOMDataFromFile<StreamReadable>("path/to/file", ig)) {
//     /* USER CODE*/
// }
// ```
// where the first `ig` lines of the given file will be ignored.
template<typename StreamReadable>
struct ReadOMDataFromFile {
    using value_type = StreamReadable;
    
    private:
    std::string path;
    const int ignore_first_n_lines;
    std::ifstream file_stream;
    bool file_found;
    bool exhausted;
    StreamReadable current_value;

    public:
    ReadOMDataFromFile(std::string path, int ignored_lines, bool exhausted_=false);
    ReadOMDataFromFile(const ReadOMDataFromFile&);

    ReadOMDataFromFile& change_path(std::string);
    StreamReadable load_next_value();
    constexpr int ignored_number_of_lines() const;
    bool file_exists() const;
    bool is_exhausted() const;
    ReadOMDataFromFile& exhaust();
    StreamReadable operator*() const;
    ReadOMDataFromFile& operator++();
    bool operator==(const ReadOMDataFromFile&) const;
    bool operator!=(const ReadOMDataFromFile&) const;
    ReadOMDataFromFile& begin();
    ReadOMDataFromFile end() const;
};

// =============================
// ReadOMDataFromFiles<R, N>
// =============================

// THIS DESCRIPTION IS INCORRECT.
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
// To construct an instance of `ReadOMDataFromFiles`,
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
// for (std::pair<int, Chirotope<R, N>> p : ReadOMDataFromFiles<R,N>
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
template<typename StreamReadable>
struct ReadOMDataFromFiles {
    using value_type = StreamReadable;

    private:
    std::string (*path_constructor)(int, int);
    int number_of_bases;
    int idx_of_file;
    bool exhausted;
    ReadOMDataFromFile<StreamReadable> file_reader;

    ReadOMDataFromFiles(std::string(*)(int,int), int, int, int);

    public:
    ReadOMDataFromFiles(const ReadOMDataFromFiles&, bool exhausted_ = false);

    ReadOMDataFromFiles(std::string(*)(int, int), int);
    
    ReadOMDataFromFile<StreamReadable>& update_file_reader();
    std::ifstream& get_file_stream();
    constexpr bool is_exhausted() const;
    constexpr int ignored_number_of_lines() const;
    
    std::pair<int, StreamReadable> operator*() const;
    ReadOMDataFromFiles& operator++();
    bool operator==(const ReadOMDataFromFiles&) const;
    bool operator!=(const ReadOMDataFromFiles&) const;

    ReadOMDataFromFiles& begin();
    ReadOMDataFromFiles end() const;
};

#include "OM_IO.cpp"