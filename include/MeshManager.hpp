#ifndef MESH_MANAGER_HPP
#define MESH_MANAGER_HPP

#include <string>
#include <iostream>
#include <map>
#include <vector>

enum MeshOrder{
    Linear = 1, 
    Quadratic = 2, 
    Cubic = 3
};


struct BoundaryCount
{
    int Count = 0;
};

struct MeshData
{
    MeshOrder Order;
    int Dim;
    int NumNodes = 0;
    int NumBoundaryElements = 0;
    int NumElements = 0;
    int NumElementsTotal = 0;
    int NumBoundaries = 0;
    int NodesPerBoundaryElement = 0; 
    int NodesPerElement = 0; 
    double **NodesVec=nullptr; 
    int **IX = nullptr; 
    int **IXB = nullptr; 
    std::map<std::string, BoundaryCount> *BoundaryElementCount=nullptr;
    std::map<std::string, BoundaryCount> *BoundaryNodeCount=nullptr;
    std::map<std::string , std::vector<int>> *BoundaryElementIdx; 
    std::map<std::string , std::vector<int>> *BoundaryNodeIdx; 

};

class MeshManager
{

public:
    MeshManager(std::string FilePath , MeshOrder Order , int Dim)
    {
        this->FilePath = FilePath;
        this->Order = Order; 
        this->Dim = Dim; 
        NumNodes = 0;
        NumBoundaryElements = 0;
        NumElements = 0;
        NumElementsTotal = 0;
        NumBoundaries = 0;
        InitStatus = false;
        std::cout << "Initializing mesh manager...\n";
        GetElementType();
        GetMeshSize();
        PrintData();
        MeshAlloc();
        ParseMesh();
        std::cout << "Mesh manager exiting succesfully.\n";
    }
    ~MeshManager()
    {
        if (InitStatus)
        {   
            for(int i = 0; i < NumNodes ; i++){
                delete[] NodesVec[i];
            }
            delete[] NodesVec;
            for(int i = 0 ; i < NumBoundaryElements ; i++){
                delete[] IXB[i]; 
            }
            delete[] IXB;
            for(int i = 0; i < NumElements ; i++){
                delete[] IX[i];
            }
            delete[] IX;
        }
    }
    int GetNumNodes();
    int GetNumElements();
    int GetNumBoundaryElements();
    int GetNumBoundaries();
    void PrintData();
    void PrintBoundaryData();
    void CreateMeshData(MeshData& Data);

private:
    std::string FilePath;
    MeshOrder Order; 
    int Dim; 
    int NumNodes;
    int NumElementsTotal;
    int NumElements;
    int NumBoundaryElements;
    int NodesPerElement;
    int NodesPerBoundaryElement;
    int NumBoundaries;
    int ElementType; 
    int BoundaryElementType;

    std::map<std::string , int> BoundaryNameToTag; 
    std::map<int , std::string > BoundaryTagToName;

    std::map<std::string, BoundaryCount> BoundaryElementCount;
    std::map<std::string, BoundaryCount> BoundaryNodeCount;

    bool InitStatus;
    double **NodesVec;
    int **IX;
    int **IXB;
    
    std::map<std::string , std::vector<int>> BoundaryElementIdx; 
    std::map<std::string , std::vector<int>> BoundaryNodeIdx; 


    void GetMeshSize();
    void ParseMesh();
    void MeshAlloc();
    void GetElementType();
};

#endif // MESH_MANAGER_HPP