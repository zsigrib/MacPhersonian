#pragma once

#include <concepts>
#include <cstddef>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <vector>
#include <string>
#include "OMs.hpp"
#include "program_template.hpp"

namespace programs {

// Utility functions, mainly for printing debug info.
namespace utility {

// Print `R`-tuples as a table with `R` rows and
// `binomial_coefficient(N, R)` columns, each column
// being an `R`-tuple.
template<int R, int N>
void print_Rtuples() {
    for (int t = 0; t < R; ++t) {
        for (int idx = 0; idx < Chirotope<R,N>::RTUPLES::NR; ++idx) {
            std::cout << int(Chirotope<R,N>::RTUPLES::LIST::array[idx][t]);
        }
        std::cout << "\n";
    }
}

// Print an iterable producing integers as a horitontal
// right-adjusted single row table, with column width of
// `spacing` many characters.
template<typename Iterable>
void print_iterable_of_ints(const Iterable& iter, int spacing=3) {
    for (auto elem: iter) {
        std::cout << std::setw(spacing) << int(elem);
    }
}

// Print an iterable producing integers as a horizontal
// right-adjusted single row table, with column width of
// `spacing` many characters - excluding the comma.
template<typename Iterable>
void print_comma_separated_iterable_of_ints(const Iterable& iter, int spacing=3) {
    bool first = true;
    for (auto elem: iter) {
        if (!first) {
            std::cout << ",";
        }
        first = false;
        std::cout << std::setw(spacing) << int(elem);
    }
}

// A helper class for printing arrays with vertically
// aligned columns. `add_cell(...)` adds a new cell to
// the right of the last row, `newline(...)` starts a
// new row.
class ArrayFormatter {
public:
    enum alignment { LEFT, CENTER, RIGHT};
private:
    struct Cell { 
        std::string str; 
        alignment align; 
        Cell(std::string s, alignment a): str(s), align(a){}
    };
    std::vector<std::vector<Cell>> cells;
public:
    ArrayFormatter() {
        cells.push_back(std::vector<Cell>());
    }
    // Convert `value` to a string and store it in a new cell at the right
    // of the last row. The cell will have alignment `align`.
    template<typename T>
    void add_cell(const T& value, enum alignment align = RIGHT) {
        cells.back().push_back(Cell(std::to_string(value), align));
    }
    void add_cell(std::string value, enum alignment align = RIGHT) {
        cells.back().push_back(Cell(value, align));
    }
    // Start a new row in the array.
    void newline() {
        cells.push_back(std::vector<Cell>());
    }
    // Print the array to std::cout, with `separator` inserted between cells
    // in the same row, and `left_pad` appearing before every row. Each row
    // will be ended by `column_end`.
    // If the array is `uniform`, all columns will have the same width.
    void print(bool uniform = false, std::string separator = " ", std::string column_end = "", std::string left_pad = "") {
        // Calculate column widths
        std::vector<size_t> column_widths;
        for (auto row: cells) {
            for (size_t i = 0; i < row.size(); ++i) {
                if (i < column_widths.size()) {
                    if (column_widths[i] < row[i].str.size()) {
                        column_widths[i] = row[i].str.size();
                    }
                } else {
                    column_widths.push_back(row[i].str.size());
                }
            }
        }
        if (uniform) {
            size_t max_column_width = 0;
            for (size_t width: column_widths) {
                if (width > max_column_width) {
                    max_column_width = width;
                }
            }
            for (size_t i = 0; i < column_widths.size(); ++i) {
                column_widths[i] = max_column_width;
            }
        }
        // Actual print
        for (auto row: cells) {
            std::cout << left_pad;
            for (size_t i = 0; i < row.size() - 1; ++i) {
                print_cell(row[i], column_widths[i]);
                std::cout << separator;
            }
            print_cell(row.back(),column_widths[row.size() - 1]);
            std::cout << column_end << "\n";
        }
    }
private:
    void print_cell(const Cell& cell, size_t column_width) {
        size_t total_margin = column_width - cell.str.size();
        size_t left_margin, right_margin;
        switch (cell.align) {
            case LEFT:
                left_margin = 0;
                right_margin = total_margin;
                break;
            case CENTER:
                left_margin = total_margin/2;
                right_margin = total_margin - left_margin;
                break;
            case RIGHT:
                left_margin = total_margin;
                right_margin = 0;
                break;
        }
        std::cout << std::string(left_margin, ' ') << cell.str << std::string(right_margin, ' ');
    }
};


}

}