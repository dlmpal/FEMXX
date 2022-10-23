#ifndef DOF_MANAGER_HPP
#define DOF_MANAGER_HPP

#include <map>
#include <MeshManager.hpp>


struct DofData{

    int NumDofPerNode;
    int NumDof;
    std::vector<std::vector<int>> *DofConnectivityVec;
    int **NodeDofMatrix;
    int NumNonZeroes; 


};

class DofManager
{
public:
    DofManager(int NumDofPerNode, MeshManager &MshManager)
    {
        std::cout<<"Initializing Dof manager...\n";
        this->NumDofPerNode = NumDofPerNode;
        MshManager.CreateMeshData(Data);
        NumDof = NumDofPerNode * Data.NumNodes;
        InitStatus = true;
        DofManagerAlloc();
        std::cout<<"Generating connectivity...\n";
        GenerateDofConnectivity();
        Print();
        std::cout<<"Dof manager exiting succesfully.\n";
    };
    ~DofManager()
    {
        if (InitStatus)
        {
            for (int i = 0; i < Data.NumNodes; i++)
            {
                delete[] NodeDofMatrix[i];
            }
            delete[] NodeDofMatrix;
        }
    };
    void CreateDofData(DofData& Data);    
    void Print();

private:
    void GenerateDofConnectivity();
    void DofManagerAlloc();

    int NumDofPerNode;
    int NumDof;
    int NumNonZeroes = 0; 
    MeshData Data;

    bool InitStatus;

    int **NodeDofMatrix;

    std::vector<std::vector<int>> DofConnectivityVec; // NumDof x (Num adjacent dof per dof)
};

#endif // DOF_MANAGER_HPP