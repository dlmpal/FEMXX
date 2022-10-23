
#include <MeshManager.hpp>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>

struct DuplicateCheck
{
    bool exists = false;
};

int MeshManager::GetNumNodes()
{
    return this->NumNodes;
}
int MeshManager::GetNumElements()
{
    return this->NumElements;
}
int MeshManager::GetNumBoundaryElements()
{
    return this->NumBoundaryElements;
}
int MeshManager::GetNumBoundaries()
{
    return this->NumBoundaries;
}
void MeshManager::PrintData()
{   if(!this->InitStatus){
    printf("MeshManager not yet initialised!\n");
    return; 
    }
    std::cout << "Number of nodes: " << GetNumNodes() << "\n";
    std::cout << "Number of boundary elements: " << GetNumBoundaryElements() << "\n";
    std::cout << "Number of elements: " << GetNumElements() << "\n";
    std::cout << "Number of boundaries: " << GetNumBoundaries() << "\n";
    std::cout << "ElementType: " <<this->ElementType <<"\n";
    std::cout << "BoundaryElementType: " <<this->BoundaryElementType <<"\n";

}

void MeshManager::PrintBoundaryData(){
    std::cout<<"Boundaries"<<"\n";
    for(auto &x : this->BoundaryNameToTag){
        std::cout << x.first << " , " <<x.second <<"\n";
    }
}
void MeshManager::GetElementType(){

    std::map<std::pair<MeshOrder , int> , std::pair<int,int>> GmshEleType; 

    GmshEleType[std::make_pair(Linear , 1)] = std::make_pair(1,2); // 2 node line
    GmshEleType[std::make_pair(Linear , 2)] = std::make_pair(2,3); // 3 node triangle
    GmshEleType[std::make_pair(Linear , 3)] = std::make_pair(4,4); // 4 node tetrahedron
    GmshEleType[std::make_pair(Quadratic , 1)] = std::make_pair(8,3); // 3 node line 
    GmshEleType[std::make_pair(Quadratic , 2)] = std::make_pair(9,6); // 6 node triangle
    GmshEleType[std::make_pair(Quadratic , 3)] = std::make_pair(11,10); // 10 node tetrahedron
    this->ElementType = GmshEleType[std::make_pair(this->Order , this->Dim)].first; 
    this->NodesPerElement = GmshEleType[std::make_pair(this->Order , this->Dim)].second; 
    this->BoundaryElementType = GmshEleType[std::make_pair(this->Order , this->Dim-1)].first;     
    this->NodesPerBoundaryElement = GmshEleType[std::make_pair(this->Order , this->Dim-1)].second; 
}

void MeshManager::GetMeshSize()
{

    // gmsh flags
    std::string PhysicalNamesFlag = "$PhysicalNames";
    std::string NodesFlag = "$Nodes";
    std::string ElementsFlag = "$Elements";
    int BoundaryElementFlag = this->BoundaryElementType;
    int ElementFlag = this->ElementType;


    std::fstream MeshFile;
    MeshFile.open(this->FilePath, std::ios::in);

    if (MeshFile.is_open())
    {
        std::string FileLine;

        double xyz;
        int empty_int;
        int NodeCount = 0;
        int BoundaryIdx;
        int CurElementType;

        while (std::getline(MeshFile, FileLine))
        {   
            if(!FileLine.compare(PhysicalNamesFlag)){
                std::getline(MeshFile,FileLine);
                this->NumBoundaries = std::stoi(FileLine)-1;
                for(int i = 0 ; i < NumBoundaries ; i++){
                    std::getline(MeshFile,FileLine); 
                    std::istringstream iss(FileLine);
                    int BoundaryTag; 
                    std::string BoundaryName; 
                    iss >> BoundaryTag; 
                    iss >> BoundaryTag; 
                    iss >> BoundaryName;
                    this->BoundaryNameToTag[BoundaryName] = BoundaryTag;
                    this->BoundaryTagToName[BoundaryTag]  = BoundaryName;
                }
            }
            if (!FileLine.compare(NodesFlag))
            {
                std::getline(MeshFile, FileLine);
                this->NumNodes = std::stoi(FileLine);
                for (int i = 0; i < this->NumNodes; i++)
                {
                    std::getline(MeshFile, FileLine);
                }
            }
            if (!FileLine.compare(ElementsFlag))
            {
                std::getline(MeshFile, FileLine);
                std::string BoundaryTag; 
                this->NumElementsTotal = std::stoi(FileLine);
                int count_boundary_ele = 0;
                int count_ele = 0;
                for (int i = 0; i < this->NumElementsTotal; i++)
                {
                    std::getline(MeshFile, FileLine);
                    std::istringstream iss(FileLine);
                    iss >> empty_int;
                    // TODO: use enum
                    iss >> CurElementType;
                    if (CurElementType == BoundaryElementFlag)
                    {
                        iss >> empty_int;
                        //iss >> empty_int;
                        iss >> BoundaryIdx;
                        BoundaryTag = this->BoundaryTagToName[BoundaryIdx];
                        BoundaryElementCount[BoundaryTag].Count++; // gmsh indexing starts at one
                        count_boundary_ele += 1; // total number of boundary elements
                    }
                    if (CurElementType == ElementFlag)
                    {
                        count_ele += 1;
                    }
                }
                this->NumBoundaryElements = count_boundary_ele;
                this->NumElements = count_ele;
               
            }
        }
        MeshFile.close();
        InitStatus = true;
    }
    else
    {
        std::cout << "Error opening mesh file , check path...\n";
        exit(-1);
    }
}

void MeshManager::MeshAlloc()
{
    if (this->InitStatus)
    {
        int dim = this->Dim;
        NodesVec = new double *[NumNodes];
        for (int i = 0; i < this->NumNodes; i++)
        {
            NodesVec[i] = new double[dim];
        }
        IX = new int *[NumElements];
        for (int i = 0; i < this->NumElements; i++)
        {
            IX[i] = new int[NodesPerElement];
        }
        IXB = new int *[NumBoundaryElements];
        for (int i = 0; i < this->NumBoundaryElements; i++)
        {
            IXB[i] = new int[NodesPerBoundaryElement];
        }
        

        for(auto &x : BoundaryElementCount){
            std::string BoundaryTag; 
            BoundaryTag = x.first; 
            int BoundaryNumEle = x.second.Count;
            BoundaryElementIdx[BoundaryTag] = std::vector<int>();
            BoundaryElementIdx[BoundaryTag].reserve(BoundaryNumEle);
            BoundaryNodeIdx[BoundaryTag] = std::vector<int>();
            BoundaryNodeIdx[BoundaryTag].reserve(NumNodes);

        }
        std::cout << "Mesh manager memory allocation complete.\n";
    }
    else
    {
        std::cout << "Mesh manager failed to allocate memory.Exiting";
        exit(-1);
    }
}

void MeshManager::ParseMesh()
{
    std::cout << "Parsing mesh file...\n";
    // gmsh flags
    std::string PhysicalNamesFlag = "$PhysicalNames";
    std::string NodesFlag = "$Nodes";
    std::string ElementsFlag = "$Elements";
    int ElementFlag = this->ElementType;
    int BoundaryElementFlag = this->BoundaryElementType;
    
    std::fstream MeshFile;
    MeshFile.open(this->FilePath, std::ios::in);

    if (MeshFile.is_open())
    {
        std::string FileLine;
        double xyz;
        int empty_int;
        int NodeCount = 0;
        int ElementType;
        int BoundaryIdx;
        int NodeIdx;
        std::map<std::pair<std::string, int>, DuplicateCheck> BoundaryNodeDuplicateCheck;

        while (std::getline(MeshFile, FileLine))
        {
            if (!FileLine.compare(NodesFlag))
            {
                std::getline(MeshFile, FileLine);
                this->NumNodes = std::stoi(FileLine);
                for (int i = 0; i < this->NumNodes; i++)
                {
                    std::getline(MeshFile, FileLine);

                    std::istringstream iss(FileLine);
                    iss >> empty_int;
                    for (int i = 0; i < this->Dim; i++)
                    {
                        iss >> this->NodesVec[NodeCount][i];
                    }
                    NodeCount++;
                }
            }
            if (!FileLine.compare(ElementsFlag))
            {
                std::getline(MeshFile, FileLine);
                this->NumElementsTotal = std::stoi(FileLine);
                int BoundaryTag; 
                std::string BoundaryName; 
                int boundary_ele_idx = 0;
                int ele_idx = 0;

                for (int i = 0; i < this->NumElementsTotal; i++)
                {
                    std::getline(MeshFile, FileLine);
                    std::istringstream iss(FileLine);
                    iss >> empty_int;
                    // TODO: use enum
                    iss >> ElementType;

                    if (ElementType == ElementFlag)
                    {
                        for (int j = 0; j < 3; j++)
                        {
                            iss >> empty_int;
                        }
                        for (int j = 0; j < NodesPerElement; j++)
                        {
                            iss >> IX[ele_idx][j];
                            IX[ele_idx][j] -= 1;
                        }
                        ele_idx += 1;
                    }
                    if (ElementType == BoundaryElementFlag)
                    {
                        iss >> empty_int; 
                        iss >> BoundaryIdx;
                        BoundaryName = this->BoundaryTagToName[BoundaryIdx];
                        // printf("face: %d ele_idx %d\n",BoundaryIdx , count_ele_2d);
                        BoundaryElementIdx[BoundaryName].push_back(boundary_ele_idx);
                        iss >> empty_int; 
                        for (int j = 0; j < NodesPerBoundaryElement; j++)
                        {
                            iss >> NodeIdx;
                            NodeIdx -= 1; // gmsh indexing starts from one
                            IXB[boundary_ele_idx][j] = NodeIdx;
                            if (!BoundaryNodeDuplicateCheck[std::make_pair(BoundaryName, NodeIdx)].exists)
                            {
                                BoundaryNodeDuplicateCheck[std::make_pair(BoundaryName, NodeIdx)].exists = true;
                                BoundaryNodeCount[BoundaryName].Count++;
                                BoundaryNodeIdx[BoundaryName].push_back(NodeIdx);
                            }
                        }
                        boundary_ele_idx += 1;
                    }
                }
                this->NumBoundaryElements = boundary_ele_idx;
                this->NumElements = ele_idx;
            }
        }
        MeshFile.close();
        // cleanup unused memory
        for(auto &x : BoundaryNameToTag){
            std::string BoundaryName = x.first; 
            BoundaryNodeIdx[BoundaryName].erase(BoundaryNodeIdx[BoundaryName].begin() + 
                            BoundaryNodeCount[BoundaryName].Count , BoundaryNodeIdx[BoundaryName].end());
        }
    }
    else
    {
        std::cout << "Error opening mesh file , check path...\n";
    }
}

void MeshManager::CreateMeshData(MeshData &Data)
{
    if (this->InitStatus)
    {
        Data.Dim = this->Dim;
        Data.Order = this->Order;
        Data.NumNodes = this->NumNodes;
        Data.NumBoundaryElements = this->NumBoundaryElements;
        Data.NumElements = this->NumElements;
        Data.NumElementsTotal = this->NumElementsTotal;
        Data.NumBoundaries = this->NumBoundaries;
        Data.NodesPerBoundaryElement = this->NodesPerBoundaryElement;
        Data.NodesPerElement = this->NodesPerElement;
        Data.NodesVec = this->NodesVec;
        Data.IX = this->IX;
        Data.IXB = this->IXB;
        Data.BoundaryElementCount = &(this->BoundaryElementCount);
        Data.BoundaryNodeCount = &(this->BoundaryNodeCount);
        Data.BoundaryElementIdx = &(this->BoundaryElementIdx);
        Data.BoundaryNodeIdx = &(this->BoundaryNodeIdx);
    }
    else
    {
        std::cout << "Mesh manager not initialized.\n";
    }
}
