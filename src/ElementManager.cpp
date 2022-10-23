#include "ElementManager.hpp"


ElementManager::ElementManager(int NGP , MeshData& meshData){

    this->NGP = NGP; 
    this->NumEle = meshData.NumElements;
    int dim = meshData.Dim; 
    int numNodes = meshData.NodesPerElement;
    this->FE.reserve(this->NumEle);
    this->FE_PROTO = ElementPrototype(dim ,numNodes , NGP , meshData.Order);
    this->NGP = FE_PROTO.NumGaussPoints;
    BuildElements(meshData);

}

ElementManager::ElementManager(const ElementManager& rhs){

    this->NGP = rhs.NGP; 
    this->NumEle = rhs.NumEle;
    this->FE_PROTO = rhs.FE_PROTO;
    this->FE.reserve(this->NumEle);
    for(int i = 0 ; i < this->NumEle ; i++){
        this->FE.push_back(rhs.FE[i]);
    }
}


ElementManager& ElementManager::operator=(const ElementManager& rhs){

    this->NGP = rhs.NGP; 
    this->NumEle = rhs.NumEle;
    this->FE_PROTO = rhs.FE_PROTO;
    this->FE.reserve(this->NumEle);
    for(int i = 0 ; i < this->NumEle ; i++){
        this->FE.push_back(rhs.FE[i]);
    }
    return *this;
}

ElementManager::~ElementManager(){

    
}


void ElementManager::BuildElements(MeshData& meshData){

    double **NodesXYZ = new double *[this->FE_PROTO.NumNodes];
    int *IX = new int[this->FE_PROTO.NumNodes];

    for (int i = 0; i < this->FE_PROTO.NumNodes; i++)
    {
        NodesXYZ[i] = new double[this->FE_PROTO.Dim];
    }

    for (int ele_idx = 0; ele_idx < meshData.NumElements; ele_idx++)
    {
        for (int i = 0; i < this->FE_PROTO.NumNodes; i++)
        {
            IX[i] = meshData.IX[ele_idx][i];
            for (int j = 0; j < FE_PROTO.Dim ;  j++)
            {
                NodesXYZ[i][j] = meshData.NodesVec[IX[i]][j];
            }
        }
        FE.push_back(Element(this->FE_PROTO.Dim , this->FE_PROTO.NumNodes , this->NGP , IX , NodesXYZ , this->FE_PROTO));
    }
    delete[] IX;
    for (int i = 0; i < this->FE_PROTO.NumNodes; i++)
    {
        delete[] NodesXYZ[i];
    }
    delete[] NodesXYZ;

}