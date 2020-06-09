#pragma once

#include <vector>
#include <iostream>

template<class T>
class Matrix2D {
private:
    size_t _rows, _cols;
    std::vector<T> matrix;

public:
    Matrix2D(size_t rows, size_t cols) :
        _rows(rows), _cols(cols), matrix(_rows*_cols) {}

    Matrix2D(size_t rows, size_t cols, T init) :
        _rows(rows), _cols(cols), matrix(_rows*_cols, init) {}

    T& operator()(size_t row, size_t col) {
        // check bounds here
        return matrix[row * _cols + col];
    }

    T operator()(size_t row, size_t col) const {
        // check bounds here
        return matrix[row * _cols + col];
    }

    size_t rows() const { return _rows; }
    size_t cols() const { return _cols; }

    friend std::ostream& operator<<(std::ostream& out, const Matrix2D& in)
    {
        for (int i = 0; i < in.rows(); i++)
        {
            for (int j = 0; j < in.cols(); j++)
            {
                out << in(i, j) << ' ';
            }
            out << std::endl;
        }
        return out;
    }
};

