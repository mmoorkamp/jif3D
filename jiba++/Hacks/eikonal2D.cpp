//============================================================================
// Name        : eikonal.cpp
// Author      : Aug 13, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================

#include <iostream>
#include <iomanip>
#include <boost/multi_array.hpp>
#include "../Global/convert.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include "../Tomo/ReadWriteTomographyData.h"

double SolveU(const double a, const double b, const double f)
  {
    double u = b + 1.0 / f;
    if (u <= a)
      return u;
    u = (a + b + sqrt(-a * a - b * b + 2 * a * b + 2.0 / (f * f))) / 2.0;
    return u;
  }

void FindMinDir(const size_t currx, const size_t curry,
    const boost::multi_array<double, 2> &U, double &a, double &b)
  {
    double minx, miny;

    minx = 1000;
    if (currx > 0)
      {
        minx = std::min(minx, U[currx - 1][curry]);
      }
    if (currx < U.shape()[0] - 1)
      {
        minx = std::min(minx, U[currx + 1][curry]);
      }

    miny = 1000;
    if (curry > 0)
      {
        miny = std::min(miny, U[currx][curry - 1]);
      }
    if (curry < U.shape()[1] - 1)
      {
        miny = std::min(miny, U[currx][curry + 1]);
      }

    b = std::min(minx, miny);
    a = std::max(minx, miny);

  }

void UpdateOne(size_t i, size_t j, boost::multi_array<double, 2> &U,
    std::vector<bool> &L, std::vector<bool> &Lnew, const boost::multi_array<
        double, 2> &Vel)
  {
    size_t offset = i * U.shape()[1] + j;
    if (!L[offset])
      {
        double p = U[i][j];
        double a, b;
        FindMinDir(i, j, U, a, b);
        double q = SolveU(a, b, Vel[i][j]);
        if (p > q)
          {
            U[i][j] = q;
            Lnew[offset] = true;

          }
      }
  }

void UpdateNeighbors(size_t i, size_t j, boost::multi_array<double, 2> &U,
    std::vector<bool> &L, std::vector<bool> &Lnew, const boost::multi_array<
        double, 2> &Vel)
  {
    if (i < U.shape()[0] - 1)
      {
        UpdateOne(i + 1, j, U, L, Lnew, Vel);
      }
    if (i > 0)
      {
        UpdateOne(i - 1, j, U, L, Lnew, Vel);
      }
    if (j < U.shape()[1] - 1)
      {
        UpdateOne(i, j + 1, U, L, Lnew, Vel);
      }
    if (j > 0)
      {
        UpdateOne(i, j - 1, U, L, Lnew, Vel);
      }

  }

int main()
  {
    const size_t nx = 7;
    const size_t ny = 7;

    const int ngrid = nx * ny;
    boost::multi_array<double, 2> U(boost::extents[nx][ny]);
    boost::multi_array<double, 2> Vel(boost::extents[nx][ny]);

    const double sourcex_index = 1;
    const double sourcey_index = 1;

    const double deltax = 1;
    const double InfValue = 1000;
    const double eps = 1e-10;
    std::fill_n(U.origin(), ngrid, InfValue);
    std::fill_n(Vel.origin(), ngrid, 1.0);
    std::vector<bool> L(ngrid, false);

    U[sourcex_index][sourcey_index] = 0.0;
    U[5][4] = 0.0;
    if (sourcex_index < nx - 1)
      {
        L[U.shape()[0] * (sourcex_index + 1) + sourcey_index] = true;

      }
    if (sourcex_index > 0)
      {
        L[U.shape()[0] * (sourcex_index - 1) + sourcey_index] = true;
      }
    if (sourcey_index < ny - 1)
      {
        L[U.shape()[0] * (sourcex_index) + sourcey_index + 1] = true;
      }
    if (sourcey_index > 0)
      {
        L[U.shape()[0] * (sourcex_index) + sourcey_index - 1] = true;
      }


    if (sourcex_index < nx - 1)
      {
        L[U.shape()[0] * (5+1) + 4] = true;

      }
    if (sourcex_index > 0)
      {
        L[U.shape()[0] * (5 - 1) + 4] = true;
      }
    if (sourcey_index < ny - 1)
      {
        L[U.shape()[0] * (5) + 4 + 1] = true;
      }
    if (sourcey_index > 0)
      {
        L[U.shape()[0] * (5) + 4 - 1] = true;
      }


    size_t nactive = 4;
    size_t iteration = 1;
    do
      {
        std::vector<bool> Lnew(L);
        for (int index = 0; index < ngrid; ++index)
          {
            if (L[index])
              {

                int j = index % U.shape()[1];
                int i = (index - j) / U.shape()[1];
                ;
                double p = U[i][j];
                double a, b;
                FindMinDir(i, j, U, a, b);
                double q = SolveU(a, b, Vel[i][j]);
                U[i][j] = q;
                if (fabs(p - q) < eps)
                  {
                    UpdateNeighbors(i, j, U, L, Lnew, Vel);

                    Lnew[index] = false;

                  }
              }
          }
        L = Lnew;
        nactive = std::accumulate(L.begin(), L.end(), 0);
        std::cout << nactive << std::endl;
        std::vector<double> Active(ngrid);
        std::copy(L.begin(), L.end(), Active.begin());
        jiba::PlotTimeField("front" + jiba::stringify(iteration) + ".vtk",
            &Active[0], deltax, nx, ny, 1);
        iteration++;
      } while (nactive > 0);
    for (size_t i = 0; i < nx ; ++i)
      {
      for (size_t j = 0; j < ny; ++j)
        {
        std::cout << std::setw(8) << std::setprecision(4) <<  U[i][j] << " ";
        }
      std::cout << std::endl;
      }
    jiba::PlotTimeField("times.vtk", U.origin(), deltax, nx, ny, 1);
    //std::cout << "U[49][49]: " << U[49][49] << std::endl;
    //std::cout << "U[4][1]: " << U[4][sourcey_index] << std::endl;
    /*std::cout << U[11][10][10] << std::endl;
     std::cout << U[10][11][10] << std::endl;
     std::cout << U[10][10][11] << std::endl;
     std::cout << U[11][11][10] << std::endl;
     std::cout << U[11][10][11] << std::endl;
     std::cout << U[10][11][11] << std::endl;
     std::cout << U[11][11][11] << std::endl;*/
  }
