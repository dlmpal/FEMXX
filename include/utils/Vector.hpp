#ifndef VECTOR_HPP 
#define VECTOR_HPP 


class Vector{

public: 
    Vector(){};
    Vector(int size);
    Vector(const Vector& rhs);
    Vector& operator=(const Vector& rhs);
    ~Vector();
    void SetAll(double val);
    void InsertElement(int idx , double val);
    double GetElement(int idx);
private: 
    int size; 
    double* vec; 
    friend class Matrix; 
    friend class LinearSystem;

};



#endif // VECTOR_HPP