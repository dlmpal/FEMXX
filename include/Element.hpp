#ifndef ELEMENT_HPP 
#define ELEMENT_HPP 

#include "Elements.hpp"

#include <vector>


class Element{

    public: 
        Element();
        Element(int dim , int NumNodes , int NumGP , int *IX , double** NodesXY , ElementPrototype& Proto);
        ~Element();
        Element(const Element& rhs); 
        Element& operator=(const Element& rhs);
        void SetVal(double val);
        double GetVal();

    private: 
        void AllocMem();
        void ComputeJac(ElementPrototype& Proto);
        void InvertJac();
        void _InvertJac2d();
        void _InvertJac3d();
        void ComputeBasisGrad(ElementPrototype& Proto);

        
    private:
        int dim , NumNodes, NumGP; 
        double Val; 
        double **NodesXYZ; 
        double ***Jac;  
        int *IX; 
        double ***JacInv; 
        double *JacDet;   
        double ***dN_dx;
        friend class ElementManager;
        friend class Solver; 
        friend class ThermalSolver;
        friend class Model;
       
    
};





#endif // ELEMENT_HPP 