#ifndef ELEMENTS_HPP
#define ELEMENTS_HPP

#include "MeshManager.hpp"

class ElementPrototype
{
public:
    int Dim;
    int NumNodes;
    int NumGaussPoints;
    double **GP;   // Gauss points: xi , eta , zeta (if Dim=3) , weight
    double **N;    // Value of shape functions at each Gauss point
    double ***Nxi; //  Value of shape function derivative at each Gauss point
    ElementPrototype()
    {
    }
    ElementPrototype(int Dim, int NumNodes, int NumGaussPoints, MeshOrder order)
    {
        this->Dim = Dim;
        this->NumNodes = NumNodes;
        this->NumGaussPoints = NumGaussPoints;
        GP = new double *[this->NumGaussPoints];
        N = new double *[this->NumGaussPoints];
        Nxi = new double **[this->NumGaussPoints];
        for (int i = 0; i < this->NumGaussPoints; i++)
        {
            GP[i] = new double[this->Dim + 1];
            N[i] = new double[this->NumNodes];
            Nxi[i] = new double *[this->NumNodes];
            for (int j = 0; j < NumNodes; j++)
            {
                Nxi[i][j] = new double[this->Dim];
            }
        }
        GetGP(order);
        GetBasis(order);
    }

    ElementPrototype(const ElementPrototype &rhs)
    {
        this->Dim = rhs.Dim;
        this->NumNodes = rhs.NumNodes;
        this->NumGaussPoints = rhs.NumGaussPoints;
        GP = new double *[this->NumGaussPoints];
        N = new double *[this->NumGaussPoints];
        Nxi = new double **[this->NumGaussPoints];

        for (int i = 0; i < this->NumGaussPoints; i++)
        {
            this->GP[i] = new double[this->Dim + 1];
            this->N[i] = new double[this->NumNodes];
            this->Nxi[i] = new double *[this->NumNodes];
            for (int j = 0; j < this->Dim + 1; j++)
            {
                this->GP[i][j] = rhs.GP[i][j];
            }
            for (int j = 0; j < NumNodes; j++)
            {
                this->N[i][j] = rhs.N[i][j];
                Nxi[i][j] = new double[this->Dim];
                for (int k = 0; k < this->Dim; k++)
                {
                    this->Nxi[i][j][k] = rhs.Nxi[i][j][k];
                }
            }
        }
    }

    ElementPrototype& operator = (const ElementPrototype& rhs)
    {

        this->Dim = rhs.Dim;
        this->NumNodes = rhs.NumNodes;
        this->NumGaussPoints = rhs.NumGaussPoints;
        GP = new double *[this->NumGaussPoints];
        N = new double *[this->NumGaussPoints];
        Nxi = new double **[this->NumGaussPoints];

        for (int i = 0; i < this->NumGaussPoints; i++)
        {
            this->GP[i] = new double[this->Dim + 1];
            this->N[i] = new double[this->NumNodes];
            this->Nxi[i] = new double *[this->NumNodes];
            for (int j = 0; j < this->Dim + 1; j++)
            {
                this->GP[i][j] = rhs.GP[i][j];
            }
            for (int j = 0; j < NumNodes; j++)
            {
                this->N[i][j] = rhs.N[i][j];
                Nxi[i][j] = new double[this->Dim];
                for (int k = 0; k < this->Dim; k++)
                {
                    this->Nxi[i][j][k] = rhs.Nxi[i][j][k];
                }
            }
        }

        return *this;
    }

    ~ElementPrototype()
    {
        for (int i = 0; i < this->NumGaussPoints; i++)
        {
            delete[] GP[i];
            delete[] N[i];

            for (int j = 0; j < NumNodes; j++)
            {
                delete[] Nxi[i][j];
            }
            delete[] Nxi[i];
        }
        delete[] GP;
        delete[] N;
        delete[] Nxi;
    }
    void Print()
    {
    }

private:
    void CopyGP2D(int Rows, double (*src)[3], double **dest)
    {
        for (int i = 0; i < Rows; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                dest[i][j] = src[i][j];
            }
        }
    }
    void CopyGP3D(int Rows, double (*src)[4], double **dest)
    {
        for (int i = 0; i < Rows; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                dest[i][j] = src[i][j];
            }
        }
    }

    void GetGP(MeshOrder order)
    {
        if (this->Dim == 2)
        { /* DIM2 */

            switch (order)
            {
            case (Linear):
            {
                double gp1[1][3] = {{1 / 3.0, 1 / 3.0, 1.0}};
                this->NumGaussPoints = 1;
                CopyGP2D(this->NumGaussPoints, gp1, GP);
                break;
            }

            case (Quadratic):
            {
                switch (NumGaussPoints)
                {

                case 7:
                {
                    double gp2[7][3] = {{0.3333333333, 0.3333333333, 0.2250000000},
                                        {0.1012865073, 0.7974269854, 0.1259391805},
                                        {0.1012865073, 0.1012865073, 0.1259391805},
                                        {0.7974269854, 0.1012865073, 0.1259391805},
                                        {0.0597158718, 0.4701420641, 0.1323941528},
                                        {0.4701420641, 0.4701420641, 0.1323941528},
                                        {0.4701420641, 0.0597158718, 0.1323941528}};
                    CopyGP2D(this->NumGaussPoints, gp2, GP);
                    break;
                }
                default:
                {
                    double gp3[3][3] = {{2.0 / 3.0, 1.0 / 6.0, 1.0 / 3.0},
                                        {1.0 / 6.0, 2.0 / 3.0, 1.0 / 3.0},
                                        {1.0 / 6.0, 1.0 / 6.0, 1.0 / 3.0}};
                    this->NumGaussPoints = 3;
                    CopyGP2D(this->NumGaussPoints, gp3, GP);
                    break;
                }
                }
            }
            break;
            default:
            {
                printf("Element order (MeshOrder) not supported\n");
                exit(-1);
            }
            }
        } /*END DIM2 */

        else if (this->Dim == 3)
        { /* DIM3 */
            switch (order)
            {
            case (Linear):
            {
                double gp4[1][4] = {{0.25, 0.25, 0.25, 0.166666667}};
                this->NumGaussPoints = 1;
                CopyGP3D(this->NumGaussPoints, gp4, GP);
                break;
            }
            case (Quadratic):
            {
                switch (NumGaussPoints)
                {

                case 8:
                {
                    double gp5[8][4] = {{0.01583591, 0.328054697, 0.328054697, 0.023087995},
                                        {0.328054697, 0.01583591, 0.328054697, 0.023087995},
                                        {0.328054697, 0.328054697, 0.01583591, 0.023087995},
                                        {0.328054697, 0.328054697, 0.328054697, 0.023087995},
                                        {0.679143178, 0.106952274, 0.106952274, 0.018578672},
                                        {0.106952274, 0.679143178, 0.106952274, 0.018578672},
                                        {0.106952274, 0.106952274, 0.679143178, 0.018578672},
                                        {0.106952274, 0.106952274, 0.106952274, 0.018578672}};
                    CopyGP3D(this->NumGaussPoints, gp5, GP);
                    break;
                }

                default:
                {
                    double gp6[4][4] = {{0.13819660, 0.13819660, 0.13819660, 0.041666667},
                                        {0.58541020, 0.13819660, 0.13819660, 0.041666667},
                                        {0.13819660, 0.58541020, 0.13819660, 0.041666667},
                                        {0.13819660, 0.13819660, 0.58541020, 0.041666667}};
                    this->NumGaussPoints = 4;
                    CopyGP3D(this->NumGaussPoints, gp6, GP);
                    break;
                }
                }
                break;
            }
            default:
            {
                printf("Element order (MeshOrder) not supported\n");
                exit(-1);
            }
            }
        } /*END DIM2 */
        else
        {
            printf("Dimension error");
            exit(-1);
        }
    }

    void GetBasis(MeshOrder order)
    {

        if (this->Dim == 2)
        {
            ComputeBasisGrad2D(order);
        }
        else if (this->Dim == 3)
        {
            ComputeBasisGrad3D(order);
        }
    }
    void ComputeBasisGrad2D(MeshOrder order)
    {
        switch (order)
        {
        case Linear:
            for (int GP_idx = 0; GP_idx < NumGaussPoints; GP_idx++)
            {
                double pt[2];
                pt[0] = GP[GP_idx][0];
                pt[1] = GP[GP_idx][1];
                double n[3];
                double nxi[6];

                n[0] = 1.0 - pt[0] - pt[1];
                n[1] = pt[0];
                n[2] = pt[1];

                nxi[0] = -1.0;
                nxi[1] = -1.0;
                nxi[2] = 1.0;
                nxi[3] = 0.0;
                nxi[4] = 0.0;
                nxi[5] = 1.0;
                int count = 0;
                for (int i = 0; i < this->NumNodes; i++)
                {
                    this->N[GP_idx][i] = n[i];
                    for (int j = 0; j < this->Dim; j++)
                    {
                        this->Nxi[GP_idx][i][j] = nxi[count];
                        count++;
                    }
                }
            }
            break;

        case Quadratic:
            for (int GP_idx = 0; GP_idx < NumGaussPoints; GP_idx++)
            {
                double pt[2];
                pt[0] = GP[GP_idx][0];
                pt[1] = GP[GP_idx][1];
                double n[6];
                double nxi[12];

                n[0] = (1.0 - pt[0] - pt[1]) * (1.0 - 2.0 * pt[0] - 2.0 * pt[1]);
                n[1] = pt[0] * (2.0 * pt[0] - 1.0);
                n[2] = pt[1] * (2.0 * pt[1] - 1.0);
                n[3] = 4.0 * pt[0] * (1.0 - pt[0] - pt[1]);
                n[4] = 4.0 * pt[0] * pt[1];
                n[5] = 4.0 * pt[1] * (1.0 - pt[0] - pt[1]);

                nxi[0] = 4.0 * pt[0] + 4.0 * pt[1] - 3.0;
                nxi[1] = 4.0 * pt[0] + 4.0 * pt[1] - 3.0;
                nxi[2] = 4.0 * pt[0] - 1.0;
                nxi[3] = 0.0;
                nxi[4] = 0.0;
                nxi[5] = 4.0 * pt[1] - 1.0;
                nxi[6] = 4.0 - 8.0 * pt[0] - 4.0 * pt[1];
                nxi[7] = -4.0 * pt[0];
                nxi[8] = 4.0 * pt[1];
                nxi[9] = 4.0 * pt[0];
                nxi[10] = -4.0 * pt[1];
                nxi[11] = 4.0 - 4.0 * pt[0] - 8.0 * pt[1];
                int count = 0;
                for (int i = 0; i < this->NumNodes; i++)
                {
                    this->N[GP_idx][i] = n[i];
                    for (int j = 0; j < this->Dim; j++)
                    {
                        this->Nxi[GP_idx][i][j] = nxi[count];
                        count++;
                    }
                }
            }
            break;

        default:
            printf("Error computing basis (2D) \n");
            exit(-1);
        }
    }
    void ComputeBasisGrad3D(MeshOrder order)
    {
        switch (order)
        {
        case Linear:
        {
            for (int GP_idx = 0; GP_idx < NumGaussPoints; GP_idx++)
            {
                double pt[3];
                pt[0] = GP[GP_idx][0];
                pt[1] = GP[GP_idx][1];
                pt[2] = GP[GP_idx][2];

                double n[4];
                double nxi[12];
                n[0] = 1 - pt[0] - pt[1] - pt[2];
                n[1] = pt[0];
                n[2] = pt[1];
                n[3] = pt[2];

                nxi[0] = -1.0;
                nxi[1] = -1.0;
                nxi[2] = -1.0;
                nxi[3] = 1.0;
                nxi[4] = 0.0;
                nxi[5] = 0.0;
                nxi[6] = 0.0;
                nxi[7] = 1.0;
                nxi[8] = 0.0;
                nxi[9] = 0.0;
                nxi[10] = 0.0;
                nxi[11] = 1.0;

                int count = 0;
                for (int i = 0; i < this->NumNodes; i++)
                {
                    this->N[GP_idx][i] = n[i];
                    for (int j = 0; j < this->Dim; j++)
                    {
                        this->Nxi[GP_idx][i][j] = nxi[count];
                        count++;
                    }
                }
            }
            break;
        }
        case Quadratic:
        {
            for (int GP_idx = 0; GP_idx < NumGaussPoints; GP_idx++)
            {
                double pt[3];
                pt[0] = GP[GP_idx][0];
                pt[1] = GP[GP_idx][1];
                pt[2] = GP[GP_idx][2];
                double n[10];
                double nxi[30];

                double l0 = 1 - pt[0] - pt[1] - pt[2];
                double l1 = pt[0];
                double l2 = pt[1];
                double l3 = pt[2];

                // Corner nodes
                n[0] = l0 * (2 * l0 - 1);
                n[1] = l1 * (2 * l1 - 1);
                n[2] = l2 * (2 * l2 - 1);
                n[3] = l3 * (2 * l3 - 1);

                // Mid-side nodes
                n[4] = 4 * l1 * l0;
                n[5] = 4 * l1 * l2;
                n[6] = 4 * l2 * l0;
                n[7] = 4 * l3 * l0;
                n[8] = 4 * l1 * l3;
                n[9] = 4 * l3 * l2;

                // Corner node derivatives
                nxi[0] = 4 * pt[0] + 4 * pt[1] + 4 * pt[2] - 3;
                nxi[1] = 4 * pt[0] + 4 * pt[1] + 4 * pt[2] - 3;
                nxi[2] = 4 * pt[0] + 4 * pt[1] + 4 * pt[2] - 3;
                nxi[3] = 4 * pt[0] - 1;
                nxi[4] = 0;
                nxi[5] = 0;
                nxi[6] = 0;
                nxi[7] = 4 * pt[1] - 1;
                nxi[8] = 0;
                nxi[9] = 0;
                nxi[10] = 0;
                nxi[11] = 4 * pt[2] - 1;

                // Mid node derivatives
                nxi[12] = -4 * (2 * pt[0] + pt[1] + pt[2] - 1);
                nxi[13] = -4 * pt[0];
                nxi[14] = -4 * pt[0];

                nxi[15] = 4 * pt[1];
                nxi[16] = 4 * pt[0];
                nxi[17] = 0.0;

                nxi[18] = -4 * pt[1];
                nxi[19] = -4 * (pt[0] + 2 * pt[1] + pt[2] - 1);
                nxi[20] = -4 * pt[1];

                nxi[21] = -4 * pt[2];
                nxi[22] = -4 * pt[2];
                nxi[23] = -4 * (pt[0] + pt[1] + 2 * pt[2] - 1);

                nxi[24] = 4 * pt[2];
                nxi[25] = 0.0;
                nxi[26] = 4 * pt[0];

                nxi[27] = 0.0;
                nxi[28] = 4 * pt[2];
                nxi[29] = 4 * pt[1];

                int count = 0;
                for (int i = 0; i < this->NumNodes; i++)
                {
                    this->N[GP_idx][i] = n[i];
                    for (int j = 0; j < this->Dim; j++)
                    {
                        this->Nxi[GP_idx][i][j] = nxi[count];
                        count++;
                    }
                }
            }
            break;
        }
        default:
            printf("Error computing basis (3D) \n");
            exit(-1);
        }
    }
};

#endif // ELEMENTS_HPP
