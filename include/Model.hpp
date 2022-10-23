#ifndef MODEL_HPP
#define MODEL_HPP

#include "MeshManager.hpp"
#include "DofManager.hpp"
#include "ElementManager.hpp"
#include "LinearSystem.hpp"
#include "utils/Vector.hpp"
#include "utils/Matrix.hpp"

enum TypeBC
{
    Dirichlet,
    Neumann
};

class Model
{
public:
    Model(MeshManager &mesh_manager, DofManager &dof_manager, int NumGP , bool Transient);
    ~Model();
    void ApplyBC(std::string BoundaryName, TypeBC Type, double Value, int *ApplyToDof);
    void AddSrc(double (*src_func)(double*,double*),int* ApplyToDof);
    void Laplace(double coeff);
    void Elasticity(double E , double nu);
    void SolveSystem();
    void TransientSetUp(double dt , double StartVal);
    void WriteSolutionToFile(std::string time);

private:
    void SetMeshAndDofData(MeshManager &mesh_manager, DofManager &dof_manager);
    void ModelAlloc();
    void AddImplicitTransientSrc();
    void AssembleRHS();
    void CreateAxb(Matrix& LHS , Vector& RHS);
    void RecoverSolution(Vector& sol);


private:
    // Mesh data
    MeshData mesh_data;

    // Dof data
    DofData dof_data;
    bool *isDofDirichlet;
    int NumFreeDof;
    int NumNonZeroes; // the number of non-zero entries in the lhs (sparse) matrix

    // FE data
    int NumGP; // number of gauss points (internal elements)
    ElementManager elementManager; // Vector of finite elements

    // System assembly data
    Vector FI; // fixed dirichlet values
    Vector FLUXI; // total source vector
    Vector EXTSRC; // external source vector
    Vector BSRC; // neumman values 
    Vector TSRC; 
    Matrix K;
    Matrix M; 

    // Linear System Data: Ax=b
    Vector Sol; 

    // Transient
    bool Transient; 
    double dt;
   

};

#endif // MODEL_HPP
