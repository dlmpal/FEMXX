#include "LinearSystem.hpp"

#include "CG.hpp"

LinearSystem::LinearSystem(int dim ){
    this->isReady = false; 
    this->max_iters = 1e4; 
    this->tolerance = 1e-4; 
    this->dim = dim; 
    this->LHS = Matrix(dim , dim);
    this->RHS = Vector(dim);
    this->Sol = Vector(dim);    
}
LinearSystem::~LinearSystem(){

}

void LinearSystem::SolveSystem(){
  
  this->LHS.Mat.ConvertToCSR();
  CG linear_solver(false);
  linear_solver.Run(this->tolerance , this->max_iters , this->Sol.vec , this->LHS.Mat , this->RHS.vec);

}

void LinearSystem::Print(){

  printf("%d %d\n",LHS.Cols , RHS.size);


}





