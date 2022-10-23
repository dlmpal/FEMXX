#include <DofManager.hpp>

struct check
{

    bool val = false;
};

void DofManager::DofManagerAlloc()
{   
    if (InitStatus)
    {
        int DofIdx=0;
        NodeDofMatrix = new int *[Data.NumNodes];
        for (int i = 0; i < Data.NumNodes; i++)
        {
            NodeDofMatrix[i] = new int[NumDofPerNode];
            for(int j = 0; j < NumDofPerNode ; j++){
                NodeDofMatrix[i][j] = DofIdx;
                DofIdx++;
            }
        }

        DofConnectivityVec.reserve(this->NumDof);

        for (int i = 0; i < NumDof; i++)
        {
            this->DofConnectivityVec.push_back(std::vector<int>());
            this->DofConnectivityVec[i].reserve(150); /////??????????????//
        }
    }
}

void DofManager::GenerateDofConnectivity()
{

    std::map<std::pair<int , int>, check> DuplicateCheck;
    int ElementNodes[this->Data.NodesPerElement];
    int NodesPerElement = this->Data.NodesPerElement;
    int NodeIdx;
    int AdjNodeIdx;
    int DofIdx;
    int AdjDofIdx;
    // Iterate through all elements
    for (int ele_idx = 0; ele_idx < this->Data.NumElements; ele_idx++)
    {
        // Get element nodes and store
        for (int j = 0; j < NodesPerElement; j++)
        {
            ElementNodes[j] = this->Data.IX[ele_idx][j];
        }
        // Iterate through the dofs of each node
        for (int i = 0; i < NodesPerElement; i++)
        {
            NodeIdx = ElementNodes[i];
            for (int j = 0; j < this->NumDofPerNode; j++)
            {
                DofIdx = NodeDofMatrix[NodeIdx][j];
                // Iterate through adjacent node dofs
                for (int ii = 0; ii < NodesPerElement; ii++)
                {
                    AdjNodeIdx = ElementNodes[ii];
                    for (int jj = 0; jj < this->NumDofPerNode; jj++)
                    {
                        AdjDofIdx = NodeDofMatrix[AdjNodeIdx][jj];
                        if (!DuplicateCheck[std::make_pair(DofIdx, AdjDofIdx)].val && AdjDofIdx != DofIdx)
                        {   
                            DuplicateCheck[std::make_pair(DofIdx, AdjDofIdx)].val = true;
                            DofConnectivityVec[DofIdx].push_back(AdjDofIdx);
                            NumNonZeroes++; 
                        }
                    }
                }
            }
        }
    }
}

void DofManager::CreateDofData(DofData& Data){
    if(InitStatus){
    Data.NumDof = this->NumDof; 
    Data.NumDofPerNode = this->NumDofPerNode; 
    Data.DofConnectivityVec = &(this->DofConnectivityVec);
    Data.NodeDofMatrix = this->NodeDofMatrix;
    Data.NumNonZeroes = this->NumNonZeroes;
    }
    else{
        std::cout<<"Dof manager not initialized.\n";
    }

}

void DofManager::Print(){
    std::cout << "Number of DOFs: " << this->NumDof << "\n";
    std::cout << "Number of non zero entries: " << this->NumNonZeroes << "\n";
}