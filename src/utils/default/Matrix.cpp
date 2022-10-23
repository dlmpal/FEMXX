
#include <cstddef>
#include "utils/Matrix.hpp"

Matrix::Matrix(int Rows , int Cols){

    this->Rows = Rows; 
    this->Cols = Cols; 
    this->Mat = SparseMatrix(Rows,Cols);
}


Matrix::Matrix(const Matrix& rhs){
    this->Rows = rhs.Rows; 
    this->Cols = rhs.Cols; 
    this->Mat = rhs.Mat;
}

Matrix& Matrix::operator=(const Matrix& rhs){
    this->Rows = rhs.Rows; 
    this->Cols = rhs.Cols; 
    this->Mat = rhs.Mat;
    return *this;
}


Matrix::~Matrix(){
}

void Matrix::InsertElement(int row_idx , int col_idx , double val){
    this->Mat.InsertElement(row_idx , col_idx , val);
}

double Matrix::GetElement(int row_idx , int col_idx){
    return this->Mat.GetElement(row_idx , col_idx);
}

Vector Matrix::MatrixVectorMult(Vector& Vec, double factor){

    if(Vec.size != this->Cols){
        printf("Matrix and vector dimensions do not match\n");
        exit(-1);
    }
    Vector result(this->Rows);
    result.vec = Mat.MatrixVectorProduct(Vec.vec , factor);
    return result;
}
