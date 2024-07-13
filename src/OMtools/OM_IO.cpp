#include <iostream>
#include <fstream>
#include <string>
#include <format>
#include "mymath.hpp"
#include "OMs.hpp"
#include "OM_IO.hpp"

// ==================
// ReadOMDataFromFile
// ==================

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

template<typename StreamReadable>
ReadOMDataFromFile<StreamReadable>::ReadOMDataFromFile
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
        current_value = load_next_value();
    }  
}

template<typename StreamReadable>
ReadOMDataFromFile<StreamReadable>::ReadOMDataFromFile
(const ReadOMDataFromFile<StreamReadable>& other):
path(other.path), file_stream(other.path.c_str()),
file_found(file_stream.good()), ignore_first_n_lines(other.ignore_first_n_lines) {
    if (file_found != other.file_found) {
        throw std::invalid_argument("Could not copy construct ReadOMDataFromFiles, as the target file was created/destroyed since the creation of the original instance.");
    } else if (other.is_exhausted()) {
        exhausted = true;
    } else {
        exhausted = false;
        ignore_n_lines(file_stream, ignore_first_n_lines);
        current_value = load_next_value();
        while (!exhausted && current_value != other.current_value) {
            current_value = load_next_value();
        }
        if (exhausted) {
            throw std::invalid_argument("Could not copy construct ReadOMDataFromFiles, because could not find current Chirotope of the original instance in the file.");
        }
    }
}

template<typename StreamReadable>
ReadOMDataFromFile<StreamReadable>& ReadOMDataFromFile<StreamReadable>::change_path(std::string new_path) {
    path = new_path;
    file_stream.close();
    file_stream.open(path.c_str());
    file_found = file_stream.good();
    exhausted = !file_found;
    if (file_found) {
        ignore_n_lines(file_stream, ignore_first_n_lines);
    }
    current_value = load_next_value();
    return (*this);
}

template<typename StreamReadable>
StreamReadable ReadOMDataFromFile<StreamReadable>::load_next_value() {
    if (exhausted) return StreamReadable();
    if (!(file_stream >> std::ws).good()) {
        exhausted = true;
        return StreamReadable();
    } else {
        StreamReadable ret; file_stream >> ret;
        return ret;
    }
}

template<typename StreamReadable>
constexpr int ReadOMDataFromFile<StreamReadable>::ignored_number_of_lines() const {
    return ignore_first_n_lines;
}

template<typename StreamReadable>
bool ReadOMDataFromFile<StreamReadable>::file_exists() const {
    return file_found;
}

template<typename StreamReadable>
bool ReadOMDataFromFile<StreamReadable>::is_exhausted() const {
    return exhausted;
}

template<typename StreamReadable>
ReadOMDataFromFile<StreamReadable>& ReadOMDataFromFile<StreamReadable>::exhaust() {
    exhausted = true;
    return (*this);
}

template<typename StreamReadable>
StreamReadable ReadOMDataFromFile<StreamReadable>::operator*() const {
    return current_value;
}

template<typename StreamReadable>
ReadOMDataFromFile<StreamReadable>& ReadOMDataFromFile<StreamReadable>::operator++() {
    current_value = load_next_value();
    return (*this);
}

template<typename StreamReadable>
bool ReadOMDataFromFile<StreamReadable>::operator==
(const ReadOMDataFromFile<StreamReadable>& other) const {
    return (exhausted && other.exhausted
        && (path == other.path) 
    ) || (!exhausted && !other.exhausted
        && (path == other.path)
        && (current_value == other.current_value)
    );
}

template<typename StreamReadable>
bool ReadOMDataFromFile<StreamReadable>::operator!=
(const ReadOMDataFromFile<StreamReadable>& other) const {
    return (!exhausted || !other.exhausted
        || (path != other.path) 
    ) && (exhausted || other.exhausted
        || (path != other.path)
        || (current_value != other.current_value)
    );
}

template<typename StreamReadable>
ReadOMDataFromFile<StreamReadable>& ReadOMDataFromFile<StreamReadable>::begin() {
    return (*this);
}

template<typename StreamReadable>
ReadOMDataFromFile<StreamReadable> ReadOMDataFromFile<StreamReadable>::end() const {
    return ReadOMDataFromFile<StreamReadable>(path, ignore_first_n_lines, true);
}

// ===================
// ReadOMDataFromFiles
// ===================

template<typename StreamReadable>
ReadOMDataFromFiles<StreamReadable>::ReadOMDataFromFiles(
    std::string(*file_path_constructor)(int,int), 
    int ignored_number_of_lines,
    int number_of_bases_, 
    int idx_of_file_
):  number_of_bases(number_of_bases_),
    idx_of_file(idx_of_file_),
    path_constructor(file_path_constructor),
    file_reader(file_path_constructor(number_of_bases_, idx_of_file_), ignored_number_of_lines),
    exhausted(false)
{
    update_file_reader();
}

template<typename StreamReadable>
ReadOMDataFromFiles<StreamReadable>::ReadOMDataFromFiles
(const ReadOMDataFromFiles& other, bool exhausted_):
file_reader(other.path_constructor(other.number_of_bases, other.idx_of_file), other.ignored_number_of_lines()) {
    exhausted = exhausted_;
    path_constructor = other.path_constructor;
    if (other.is_exhausted() || exhausted_) {
        number_of_bases = 0;
        idx_of_file = -1;
    } else {
        number_of_bases = other.number_of_bases;
        idx_of_file = other.idx_of_file;
        update_file_reader();
    }
}

template<typename StreamReadable>
ReadOMDataFromFiles<StreamReadable>::ReadOMDataFromFiles
(std::string(*file_path_constructor)(int, int), int ignored_number_of_lines): 
file_reader(file_path_constructor(1, 0), ignored_number_of_lines) {
    exhausted = false;
    path_constructor = file_path_constructor;
    number_of_bases = 1;
    idx_of_file = 0;
    update_file_reader();
}

template<typename StreamReadable>
ReadOMDataFromFile<StreamReadable>& ReadOMDataFromFiles<StreamReadable>::update_file_reader() {
    if (exhausted) return file_reader;
    file_reader.change_path(path_constructor(number_of_bases, idx_of_file));
    if (file_reader.file_exists() || idx_of_file == 0) {
        exhausted = !file_reader.file_exists();
        return file_reader;
    }
    idx_of_file = 0;
    number_of_bases++;
    file_reader.change_path(path_constructor(number_of_bases, idx_of_file));
    while (file_reader.file_exists() && file_reader.is_exhausted()) {
        number_of_bases++;
        file_reader.change_path(path_constructor(number_of_bases, idx_of_file));
    }
    exhausted = !file_reader.file_exists();
    return file_reader;
}

template<typename StreamReadable>
constexpr bool ReadOMDataFromFiles<StreamReadable>::is_exhausted() const {
    return exhausted;
}

template<typename StreamReadable>
constexpr int ReadOMDataFromFiles<StreamReadable>::ignored_number_of_lines() const {
    return file_reader.ignored_number_of_lines();
}

template<typename StreamReadable>
std::pair<int, StreamReadable> ReadOMDataFromFiles<StreamReadable>::operator*() const {
    return std::pair<int, StreamReadable>(number_of_bases, *file_reader);
}

template<typename StreamReadable>
ReadOMDataFromFiles<StreamReadable>& ReadOMDataFromFiles<StreamReadable>::operator++() {
    ++file_reader;
    if (file_reader.is_exhausted()) {
        idx_of_file++;
        update_file_reader();
    } 
    return (*this);
}

template<typename StreamReadable>
bool ReadOMDataFromFiles<StreamReadable>::operator==(const ReadOMDataFromFiles& other) const {
    return (!is_exhausted() && !other.is_exhausted()
        && (number_of_bases == other.number_of_bases)
        && (idx_of_file == other.idx_of_file)
        && (file_reader == other.file_reader)
        && (path_constructor) == (other.path_constructor)
    ) || (is_exhausted() && other.is_exhausted()
        && path_constructor == other.path_constructor
    );
}

template<typename StreamReadable>
bool ReadOMDataFromFiles<StreamReadable>::operator!=(const ReadOMDataFromFiles& other) const {
    return (is_exhausted() || other.is_exhausted()
        || (number_of_bases != other.number_of_bases)
        || (idx_of_file != other.idx_of_file)
        || (file_reader != other.file_reader)
        || (path_constructor) != (other.path_constructor)
    ) && (!is_exhausted() || !other.is_exhausted()
        || path_constructor != other.path_constructor
    );
}

template<typename StreamReadable>
ReadOMDataFromFiles<StreamReadable>& ReadOMDataFromFiles<StreamReadable>::begin() {
    return (*this);
}

template<typename StreamReadable>
ReadOMDataFromFiles<StreamReadable> ReadOMDataFromFiles<StreamReadable>::end() const {
    return ReadOMDataFromFiles<StreamReadable>(
        *this,
        true
    );
}
