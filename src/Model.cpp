

#include "Model.hpp"
#include "Elements.hpp"
#include "SparseMatrix.hpp"
#include "CG.hpp"
#include "GaussSeidel.hpp"

#include <stdio.h>
#include <math.h>
#include <fstream>

////////////////////////Memory Allocation and initialization////////////////////////////////////////////

Model::Model(MeshManager &mesh_manager, DofManager &dof_manager, int NumGP, bool Transient)
{
    this->NumGP = NumGP;
    this->Transient = Transient;
    SetMeshAndDofData(mesh_manager, dof_manager);
    NumFreeDof = dof_data.NumDof;
    NumNonZeroes = dof_data.NumNonZeroes;
    ModelAlloc();
    elementManager = ElementManager(NumGP, mesh_data);
};

void Model::SetMeshAndDofData(MeshManager &mesh_manager, DofManager &dof_manager)
{
    mesh_manager.CreateMeshData(mesh_data);
    dof_manager.CreateDofData(dof_data);
}
void Model::TransientSetUp(double dt, double StartVal)
{
    this->dt = dt;
    this->Sol.SetAll(StartVal);
}

void Model::ModelAlloc()
{

    FI = Vector(dof_data.NumDof);
    FLUXI = Vector(dof_data.NumDof);
    EXTSRC = Vector(dof_data.NumDof);
    BSRC = Vector(dof_data.NumDof);
    TSRC = Vector(dof_data.NumDof);
    Sol = Vector(dof_data.NumDof);
    isDofDirichlet = new bool[dof_data.NumDof];
    for (int i = 0; i < dof_data.NumDof; i++)
    {
        isDofDirichlet[i] = false;
    }

    K = Matrix(dof_data.NumDof, dof_data.NumDof);
    M = Matrix(dof_data.NumDof, dof_data.NumDof);
}

Model::~Model()
{

    delete[] isDofDirichlet;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////Boundary conditions and sources////////////////////////////////////////////

void Model::ApplyBC(std::string BoundaryName, TypeBC type, double Value, int *ApplyToDof)
{

    int node_idx;
    int dof_idx;
    if (type == Dirichlet)
    {
        for (int i = 0; i < mesh_data.BoundaryNodeCount[0][BoundaryName].Count; i++)
        {
            node_idx = mesh_data.BoundaryNodeIdx[0][BoundaryName][i];
            for (int j = 0; j < dof_data.NumDofPerNode; j++)
            {

                dof_idx = dof_data.NodeDofMatrix[node_idx][j];
                this->FI.InsertElement(dof_idx, Value * ((double)ApplyToDof[j]));
                if (!this->isDofDirichlet[dof_idx])
                {
                    NumFreeDof -= 1 * ApplyToDof[j];
                }
                isDofDirichlet[dof_idx] = (bool)ApplyToDof[j];
            }
        }
    }
    if (type == Neumann)
    { // TODO: need to loop over elements and integrate...
        for (int i = 0; i < mesh_data.BoundaryNodeCount[0][BoundaryName].Count; i++)
        {

            node_idx = mesh_data.BoundaryNodeIdx[0][BoundaryName][i];
            for (int j = 0; j < dof_data.NumDofPerNode; j++)
            {

                dof_idx = dof_data.NodeDofMatrix[node_idx][j];
                BSRC.InsertElement(dof_idx, Value * ((double)ApplyToDof[j]));
            }
        }
    }
}

void Model::AddSrc(double (*func)(double *, double *), int *ApplyToDof)
{

    int NumNodes = mesh_data.NodesPerElement;
    int NumDof = dof_data.NumDofPerNode;
    int NumGP = elementManager.NGP;
    int dim = mesh_data.Dim;

    for (auto &fe : this->elementManager.FE)
    {
        for (int i = 0; i < NumNodes; i++)
        {
            int node_idx = fe.IX[i];
            double node_xyz[dim];
            for (int j = 0; j < dim; j++)
            {
                node_xyz[j] = mesh_data.NodesVec[node_idx][j];
            }

            for (int j = 0; j < NumDof; j++)
            {
                int dof_idx = this->dof_data.NodeDofMatrix[node_idx][j];
                double temp = 0;
                for (int igp = 0; igp < NumGP; igp++)
                {
                    temp += ApplyToDof[j] * func(this->elementManager.FE_PROTO.GP[igp], node_xyz) * fe.JacDet[igp] * this->elementManager.FE_PROTO.GP[igp][dim];
                }
                EXTSRC.InsertElement(dof_idx, temp);
            }
        }
    }
}

void Model::AssembleRHS()
{

    this->FLUXI.SetAll(0.0);
    for (int i = 0; i < this->dof_data.NumDof; i++)
    {
        this->FLUXI.InsertElement(i, this->EXTSRC.GetElement(i) + BSRC.GetElement(i) + this->TSRC.GetElement(i));
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////System Assembly and Solution////////////////////////////////////////////

void Model::Laplace(double coeff)
{
    int dim = this->mesh_data.Dim;
    int NumNodes = this->mesh_data.NodesPerElement;
    int NGP = this->elementManager.NGP;

    double rho_cp = 8000 * 400;

    for (int ele_idx = 0; ele_idx < this->elementManager.NumEle; ele_idx++)
    {
        Element fe = this->elementManager.FE[ele_idx];
        fe.SetVal(coeff);
        double KELE[NumNodes][NumNodes];
        double MELE[NumNodes][NumNodes];
        for (int i = 0; i < NumNodes; i++)
        {
            for (int j = 0; j < NumNodes; j++)
            {
                KELE[i][j] = 0.0;
                MELE[i][j] = 0.0;
            }
        }
        for (int igp = 0; igp < NGP; igp++)
        {
            for (int i = 0; i < NumNodes; i++)
            {
                for (int j = 0; j < NumNodes; j++)
                {
                    double tempK = 0.0;
                    double tempM = 0.0;
                    for (int k = 0; k < dim; k++)
                    {
                        tempK += fe.dN_dx[igp][i][k] * fe.dN_dx[igp][j][k];
                    }
                    tempK *= fe.JacDet[igp] * this->elementManager.FE_PROTO.GP[igp][dim] * fe.GetVal();
                    if (Transient)
                    {
                        tempM = this->elementManager.FE_PROTO.N[igp][i] * this->elementManager.FE_PROTO.N[igp][j] * rho_cp * fe.JacDet[igp] * elementManager.FE_PROTO.GP[igp][dim];
                        MELE[i][j] += tempM;
                    }
                    KELE[i][j] += tempK;
                }
            }
        }
        double tempM;
        double tempK;
        for (int i = 0; i < NumNodes; i++)
        {
            for (int j = 0; j < NumNodes; j++)
            {
                tempK = K.GetElement(fe.IX[i], fe.IX[j]) + KELE[i][j]; //+ MELE[i][j]/dt; // not general , fix later
                K.InsertElement(fe.IX[i], fe.IX[j], tempK);
                tempM = M.GetElement(fe.IX[i], fe.IX[j]) + MELE[i][j];
                M.InsertElement(fe.IX[i], fe.IX[j], tempM);
            }
        }
    }
};

void Model::Elasticity(double E , double nu)
{

    int NumNodes = this->mesh_data.NodesPerElement;
    int NumDof = this->dof_data.NumDofPerNode * NumNodes;
    int dim = this->mesh_data.Dim;

    double **D = new double *[3 * (dim - 1)];
    for (int i = 0; i < 3 * (dim - 1); i++)
    {
        D[i] = new double[3 * (dim - 1)];
    }
    if (dim == 2)
    {
        double D2[3][3] = {{1, nu, 0}, {nu, 1, 0}, {0, 0, (1 - nu) / 2}};
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                D[i][j] = D2[i][j] * E / (1 - nu * nu);
            }
        }
    }
    if (dim == 3)
    {   
        double d11 = nu;
        double d33 = 1 - nu;
        double d22 = (1 - 2 * nu) / 2;
        double D3[6][6] = {{d33, d11, d11, 0, 0, 0}, {d11, d33, d11, 0, 0, 0}, {d11, d11, d33, 0, 0, 0}, {0, 0, 0, d22, 0, 0}, {0, 0, 0, 0, d22, 0}, {0, 0, 0, 0, 0, d22}};
        for (int i = 0; i < 6; i++)
        {
            for (int j = 0; j < 6; j++)
            {
                D[i][j] = D3[i][j] * E / ((1 + nu) * (1 - 2 * nu));
            }
        }
    }

    for (int ele_idx = 0; ele_idx < this->elementManager.NumEle; ele_idx++)
    {
        Element fe = this->elementManager.FE[ele_idx];

        double KELE[NumDof][NumDof];
        for (int i = 0; i < NumDof; i++)
        {
            for (int j = 0; j < NumDof; j++)
            {
                KELE[i][j] = 0.0;
            }
        }

        for (int igp = 0; igp < this->elementManager.NGP; igp++)
        {   

            double DB[3 * (dim - 1)][NumDof];
            double B[3 * (dim - 1)][NumDof];

            for (int i = 0; i < 3 * (dim - 1); i++)
            {
                for (int j = 0; j < NumDof; j++)
                {
                    B[i][j] = 0;
                    DB[i][j] = 0;
                }
            }

            for (int i = 0; i < dim; i++)
            {
                int b_col_idx = i;
                for (int j = 0; j < NumNodes; j++)
                {
                    B[i][b_col_idx] = fe.dN_dx[igp][j][i];
                    b_col_idx += dim;
                }
            }

            int dx_inds[][2] = {{1, 0}, {2, 1}, {2, 0}};
            int local_inds[][2] = {{0, 1}, {1, 2}, {0, 2}};

            for (int i = dim; i < 3 * (dim - 1); i++)
            {
                int col_ind = 0;
                for (int j = 0; j < NumNodes; j++)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        B[i][local_inds[i - dim][k] + col_ind] = fe.dN_dx[igp][j][dx_inds[i - dim][k]];
                    }
                    col_ind += dim;
                }
            }
            
            for (int i = 0; i < 3 * (dim - 1); i++)
            {
                for (int j = 0; j < NumDof; j++)
                {
                    double temp = 0;
                    for (int k = 0; k < 3 * (dim - 1); k++)
                    {
                        temp += D[i][k] * B[k][j];
                    }
                    DB[i][j] = temp;
                }
            }
            for (int i = 0; i < NumDof; i++)
            {
                for (int j = 0; j < NumDof; j++)
                {
                    double temp = 0;
                    for (int k = 0; k < 3 * (dim - 1); k++)
                    {
                        temp += B[k][i] * DB[k][j];
                    }
                    KELE[i][j] += temp * fe.JacDet[igp] * this->elementManager.FE_PROTO.GP[igp][this->mesh_data.Dim];
                }
            }
        }

        int EleDof[NumDof];
        int count = 0;
        for (int i = 0; i < NumNodes; i++)
        {
            for (int j = 0; j < dim; j++)
            {
                EleDof[count] = this->dof_data.NodeDofMatrix[fe.IX[i]][j];
                count++;
            }
        }

        for (int i = 0; i < NumDof; i++)
        {
            for (int j = 0; j < NumDof; j++)
            {
                double temp = this->K.GetElement(EleDof[i], EleDof[j]);
                this->K.InsertElement(EleDof[i], EleDof[j], temp + KELE[i][j]);
            }
        }
    }
    
    for(int i = 0 ; i < 3*(dim-1) ; i++){
    delete[] D[i];
    }
    delete[] D;
}

void Model::AddImplicitTransientSrc()
{

    if (!this->Transient)
    {
        return;
    }
    this->TSRC = M.MatrixVectorMult(Sol, 1.0 / this->dt);
}

void Model::CreateAxb(Matrix &LHS, Vector &RHS)
{

    int *DofEqIdx = new int[dof_data.NumDof];
    int eq_idx = 0;
    for (int i = 0; i < dof_data.NumDof; i++)
    {
        DofEqIdx[i] = -1;
        if (!isDofDirichlet[i])
        {
            DofEqIdx[i] = eq_idx;
            eq_idx++;
        }
    }

    for (int dof_idx = 0; dof_idx < dof_data.NumDof; dof_idx++)
    {
        int dof_eq_idx = DofEqIdx[dof_idx];
        if (!isDofDirichlet[dof_idx])
        {
            // RHS[dof_eq_idx] = FLUXI[dof_idx];
            LHS.InsertElement(dof_eq_idx, dof_eq_idx, K.GetElement(dof_idx, dof_idx));
            RHS.InsertElement(dof_eq_idx, FLUXI.GetElement(dof_idx));

            for (auto &adj_dof_idx : dof_data.DofConnectivityVec[0][dof_idx])
            {
                if (!isDofDirichlet[adj_dof_idx])
                {
                    int adj_dof_eq_idx = DofEqIdx[adj_dof_idx];
                    LHS.InsertElement(dof_eq_idx, adj_dof_eq_idx, K.GetElement(dof_idx, adj_dof_idx));
                }
                else
                {
                    // RHS[dof_eq_idx] -= K.GetElement(dof_idx, adj_dof_idx) * FI[adj_dof_idx];
                    double temp = RHS.GetElement(dof_eq_idx) - K.GetElement(dof_idx, adj_dof_idx) * FI.GetElement(adj_dof_idx);
                    RHS.InsertElement(dof_eq_idx, temp);
                }
            }
        }
    }

    delete[] DofEqIdx;
}

void Model::RecoverSolution(Vector &sol)
{

    int eq_idx = 0;
    for (int i = 0; i < dof_data.NumDof; i++)
    {
        double temp = isDofDirichlet[i] ? FI.GetElement(i) : sol.GetElement(eq_idx);
        this->Sol.InsertElement(i, temp);
        if (!isDofDirichlet[i])
        {
            eq_idx++;
        }
    }
}

void Model::SolveSystem()
{
    // AddImplicitTransientSrc();
    AssembleRHS();
    LinearSystem System(this->NumFreeDof);
    CreateAxb(System.LHS, System.RHS);
    System.SolveSystem();
    RecoverSolution(System.Sol);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

void Model::WriteSolutionToFile(std::string time)
{

    std::string NamesU[] = {"Ux", "Uy", "Uz"};
    std::string NamesXYZ[] = {"x", "y", "z"};

    std::fstream SolutionFile;
    std::string FilePath = "U" + time + ".csv";
    SolutionFile.open(FilePath, std::ios::out);
    for (int i = 0; i < mesh_data.Dim; i++)
    {
        SolutionFile << NamesXYZ[i] + ",";
    }

    int i = 0;
    for (i = 0; i < dof_data.NumDofPerNode - 1; i++)
    {
        SolutionFile << NamesU[i] + ",";
    }
    SolutionFile << NamesU[i];
    SolutionFile << "\n";
    int count = 0;
    for (int j = 0; j < mesh_data.NumNodes; j++)
    {
        for (int k = 0; k < mesh_data.Dim; k++)
        {
            SolutionFile << mesh_data.NodesVec[j][k] << ",";
        }
        for (int k = 0; k < dof_data.NumDofPerNode; k++)
        {
            SolutionFile << this->Sol.GetElement(count) << ",";
            count++;
        }
        SolutionFile << "\n";
    }
    SolutionFile.close();
}
