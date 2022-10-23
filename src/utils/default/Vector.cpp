
#include "utils/Vector.hpp"


Vector::Vector(int size){
    this->size = size; 
    this->vec = new double[size];
    for(int i = 0 ; i<size;i++){
        vec[i] = 0; 
    }
}

Vector::Vector(const Vector& rhs){
    this->size = rhs.size; 
    this->vec = new double[size];
    for(int i = 0 ; i < this->size ; i++){
        this->vec[i] = rhs.vec[i];
    }

}

Vector& Vector::operator=(const Vector& rhs){
    this->size = rhs.size; 
    this->vec = new double[size];
    for(int i = 0 ; i < this->size ; i++){
        this->vec[i] = rhs.vec[i];
    }
    return *this;
}

Vector::~Vector(){
    delete[] vec; 
}

double Vector::GetElement(int idx){
    return this->vec[idx];
}

void Vector::InsertElement(int idx , double val ){

    vec[idx] = val; 
}

void Vector::SetAll(double val){
    for(int i = 0 ; i < size ; i++){
        vec[i] = val;
    }
}
