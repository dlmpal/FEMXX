

#include "Element.hpp"

Element::Element()
{
}

Element::Element(int dim, int NumNodes ,  int NumGP, int *IX , double** NodesXY, ElementPrototype &Proto)
{
    this->dim = dim;
    this->NumNodes = NumNodes;
    this->NumGP = NumGP;

    AllocMem();

    for(int i = 0 ; i < NumNodes ; i++){
        this->IX[i] = IX[i];
        for(int j = 0 ; j < dim ; j++){
            this->NodesXYZ[i][j] = NodesXY[i][j];
        }
    }

    
    ComputeJac(Proto);
    InvertJac();
    ComputeBasisGrad(Proto);

}

void Element::AllocMem(){
    this->IX = new int[NumNodes];
    this->NodesXYZ = new double *[NumNodes];
    this->Jac = new double **[NumGP];
    this->JacInv = new double **[NumGP];
    this->JacDet = new double[NumGP];
    this->dN_dx  = new double**[NumGP];
    for (int i = 0; i < this->NumNodes; i++)
    {
        this->NodesXYZ[i] = new double[dim];
    }

    for (int i = 0; i < NumGP; i++)
    {
        this->Jac[i] = new double *[dim];
        this->JacInv[i] = new double *[dim];
        this->dN_dx[i] = new double*[NumNodes];
        for (int j = 0; j < dim; j++)
        {
            this->Jac[i][j] = new double[dim];
            this->JacInv[i][j] = new double[dim];
        }
        for(int j = 0; j < NumNodes ; j++){
            this->dN_dx[i][j] = new double[dim];
        }
    }

}


Element::Element(const Element& rhs){
    this->dim = rhs.dim; 
    this->NumNodes = rhs.NumNodes; 
    this->NumGP = rhs.NumGP;
    this->Val = rhs.Val;
    AllocMem();
    for(int i = 0 ; i < NumNodes ; i++){
        this->IX[i] = rhs.IX[i];
        for(int j = 0 ; j < dim ; j++){
            this->NodesXYZ[i][j] = rhs.NodesXYZ[i][j];
        }
    }
    for(int igp = 0 ; igp < NumGP ; igp++){
        this->JacDet[igp] = rhs.JacDet[igp];
        for(int i = 0 ; i < dim ; i++){
            for(int j = 0 ; j < dim ; j++){
                this->Jac[igp][i][j] = rhs.Jac[igp][i][j];
                this->JacInv[igp][i][j] = rhs.JacInv[igp][i][j];
            }
        }
        for(int i = 0; i < NumNodes; i++){
            for(int j = 0 ; j < dim ; j++){
                this->dN_dx[igp][i][j] = rhs.dN_dx[igp][i][j];
            }
        }
    }
}

Element& Element::operator=(const Element& rhs){
    this->dim = rhs.dim; 
    this->NumNodes = rhs.NumNodes; 
    this->NumGP = rhs.NumGP;
    this->Val = rhs.Val;
    AllocMem();
    for(int i = 0 ; i < NumNodes ; i++){
        this->IX[i] = rhs.IX[i];
        for(int j = 0 ; j < dim ; j++){
            this->NodesXYZ[i][j] = rhs.NodesXYZ[i][j];
        }
    }
    for(int igp = 0 ; igp < NumGP ; igp++){
        this->JacDet[igp] = rhs.JacDet[igp];
        for(int i = 0 ; i < dim ; i++){
            for(int j = 0 ; j < dim ; j++){
                this->Jac[igp][i][j] = rhs.Jac[igp][i][j];
                this->JacInv[igp][i][j] = rhs.JacInv[igp][i][j];
            }
        }
        for(int i = 0; i < NumNodes; i++){
            for(int j = 0 ; j < dim ; j++){
                this->dN_dx[igp][i][j] = rhs.dN_dx[igp][i][j];
            }
        }
    }
    return *this;
}



Element::~Element(){
    delete[] IX; 
    for(int i = 0 ; i < NumNodes ; i++){
        delete[] NodesXYZ[i];
    }
    delete[] NodesXYZ;
    delete[] JacDet;

    for(int i = 0 ; i < NumGP ; i++){
        for(int j = 0 ; j < NumNodes ; j++){
            delete[]  dN_dx[i][j];
        }
        for(int j = 0 ; j < dim; j++){
            delete[] Jac[i][j];
            delete[] JacInv[i][j];
        }
        delete[] dN_dx[i];
        delete[] Jac[i];
        delete[] JacInv[i];
    }
    delete[] Jac;
    delete[] JacInv;
    delete[] dN_dx;

}

void Element::ComputeJac(ElementPrototype &Proto)
{

    for (int GP_idx = 0; GP_idx < NumGP; GP_idx++)
    {
        for (int i = 0; i < dim; i++)
        {
            for (int j = 0; j < dim; j++)
            {
                this->Jac[GP_idx][i][j] = 0;
                for (int k = 0; k < NumNodes; k++)
                {
                    this->Jac[GP_idx][i][j] += this->NodesXYZ[k][i] * Proto.Nxi[GP_idx][k][j]; // J[i][j] = dx_i / dksi_j = Sum dNk/dksi_j * x_k
                }
            }
        }
    }
}

void Element::InvertJac()
{

    if (this->dim == 2)
    {
        this->_InvertJac2d();
    }
    else if (this->dim == 3)
        this->_InvertJac3d();
}

void Element::_InvertJac2d()
{

    for (int GP_idx = 0; GP_idx < NumGP; GP_idx++)
    {

        double det_J = this->Jac[GP_idx][0][0] * this->Jac[GP_idx][1][1] - this->Jac[GP_idx][1][0] * this->Jac[GP_idx][0][1];
        double det_J_inv = 1.0 / det_J;
        this->JacInv[GP_idx][0][0] = this->Jac[GP_idx][1][1] * det_J_inv;
        this->JacInv[GP_idx][0][1] = -this->Jac[GP_idx][0][1] * det_J_inv;
        this->JacInv[GP_idx][1][0] = -this->Jac[GP_idx][1][0] * det_J_inv;
        this->JacInv[GP_idx][1][1] = this->Jac[GP_idx][0][0] * det_J_inv;
        this->JacDet[GP_idx] = det_J;

    }
}

void Element::_InvertJac3d(){

    for(int igp = 0; igp < NumGP ; igp++){
     double det_J = Jac[igp][0][0] * (Jac[igp][1][1] * Jac[igp][2][2] - Jac[igp][2][1] * Jac[igp][1][2]) - 
     Jac[igp][0][1] * (Jac[igp][1][0] * Jac[igp][2][2] - Jac[igp][2][0] * Jac[igp][1][2]) + 
     Jac[igp][0][2] * (Jac[igp][1][0] * Jac[igp][2][1] - Jac[igp][1][1] * Jac[igp][2][0]);

     double det_J_inv = 1.0 / det_J;

    double J_inv_3d[3][3] = {{Jac[igp][1][1] * Jac[igp][2][2] - Jac[igp][1][2] * Jac[igp][2][1], -(Jac[igp][0][1] * Jac[igp][2][2] - Jac[igp][0][2] * Jac[igp][2][1]), Jac[igp][0][1] * Jac[igp][1][2] - Jac[igp][0][2] * Jac[igp][1][1]},
                                         {-(Jac[igp][1][0] * Jac[igp][2][2] - Jac[igp][1][2] * Jac[igp][2][0]), Jac[igp][0][0] * Jac[igp][2][2] - Jac[igp][0][2] * Jac[igp][2][0], -(Jac[igp][0][0] * Jac[igp][1][2] - Jac[igp][0][2] * Jac[igp][1][0])},
                                         {Jac[igp][1][0] * Jac[igp][2][1] - Jac[igp][1][1] * Jac[igp][2][0], -(Jac[igp][0][0] * Jac[igp][2][1] - Jac[igp][0][1] * Jac[igp][2][0]), Jac[igp][0][0] * Jac[igp][1][1] - Jac[igp][0][1] * Jac[igp][1][0]}};
                for (int i = 0; i < dim; i++)
                {
                    for (int j = 0; j < dim; j++)
                    {
                        JacInv[igp][i][j] = J_inv_3d[i][j] * det_J_inv;
                    }
                }
    this->JacDet[igp] = det_J;
    }

}

void Element::ComputeBasisGrad(ElementPrototype& Proto){
    for(int igp=0 ;igp<NumGP;igp++){
     for (int i = 0; i < this->NumNodes; i++)

            {

                for (int j = 0; j < dim; j++)
                {
                    dN_dx[igp][i][j] = 0.0;
                    
                    for (int k = 0; k < dim; k++)
                    {
                        dN_dx[igp][i][j] += Proto.Nxi[igp][i][k]* this->JacInv[igp][k][j];

                    }
                    
                }

            }

    }


}


void Element::SetVal(double val ){
    this->Val = val; 
}

double Element::GetVal(){
    return this->Val;
}