#ifndef LINEARSYSTEM_HPP 
#define LINEARSYSTEM_HPP 
#include <cstddef>

#include "utils/Matrix.hpp"
#include "utils/Vector.hpp"

class LinearSystem{


public: 
    LinearSystem(){};
    LinearSystem(int dim);
    ~LinearSystem();
    void SolveSystem();
    void Print();


private: 
    int dim; 
    double tolerance; 
    double max_iters; 
    Matrix LHS; 
    Vector RHS; 
    Vector Sol; 
    bool isReady; 
    friend class Model;
};





#endif // LINEARSYSTEM_HPP