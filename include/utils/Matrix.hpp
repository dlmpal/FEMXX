#ifndef MATRIX_HPP 
#define MATRIX_HPP 

#include "Vector.hpp"
#include "SparseMatrix.hpp"

class Matrix{



public:
    Matrix(){};
    Matrix(int Rows , int Col);
    Matrix(const Matrix& rhs);
    Matrix& operator=(const Matrix& rhs);
    ~Matrix();
    void InsertElement(int row_idx , int col_idx , double val);
    double GetElement(int row_idx , int col_idx);
    Vector MatrixVectorMult(Vector& Vec , double factor);
private: 
    int Rows; 
    int Cols;
    SparseMatrix Mat; 
    friend class Vector;
    friend class LinearSystem;

};






#endif // MATRIX_HPP 