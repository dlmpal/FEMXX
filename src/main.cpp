
#include <MeshManager.hpp>
#include <DofManager.hpp>
#include <Model.hpp>
#include <math.h>
#include <fstream>


double GaussianSource(double *ksi_eta , double *xy){

    double val = 0; 
    double x = xy[0];
    double y = xy[1];
    double dist1 = sqrt(pow(x-5,2));
    if(dist1 < 0.5){
        val = -50000;
    }
   
    return val ; 
}

// Write neumann bc 
// Implement correct timestepping
// Compile files into libs
// Add quadratic 3d ele

int main(int argc, char *argv[])
{

    std::string MeshFilePath = argv[1];
    int dim = atoi(argv[2]);
    MeshOrder order = static_cast<MeshOrder>(atoi(argv[3]));

    MeshManager mesh_manager(MeshFilePath, order, dim);
    DofManager dof_manager(3, mesh_manager);


    int choose1[] = {1, 1 ,1}; // ux , uy , uz == fixed
    int choose2[] = {0 , 1  , 0}; // ux == fixed , uy , uz == free
    Model Structural(mesh_manager , dof_manager , 1 , false);
    Structural.ApplyBC("FixedSupport",Dirichlet, 0 , choose1);
    Structural.ApplyBC("LoadedEnd",Dirichlet, 0 , choose1);
    Structural.AddSrc(GaussianSource , choose2);
    Structural.Elasticity(1e8 , 0.3);
    Structural.SolveSystem();
    Structural.WriteSolutionToFile("");
    printf("Finished.\n");

    return 0;
};