//============================================================================
// Name        : eikonal.cpp
// Author      : Aug 13, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================

#include <iostream>
#include <boost/multi_array.hpp>
#include "../Global/convert.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include "../Tomo/ReadWriteTomographyData.h"

double SolveU(const double a, const double b, const double c, const double f)
  {
    double u = c + 1.0 / f;
    if (u <= b)
      return u;
    u = (b + c + sqrt(-b * b - c * c + 2 * b * c + 2.0 / (f * f))) / 2.0;
    if (u <= a)
      return u;
    u = (2 * (a + b + c) + sqrt(4 * pow(a + b + c, 2) - 12 * (a * a + b * b + c
        * c - 1.0 / (f * f)))) / 6.0;
    return u;
  }

void FindMinDir(const size_t currx, const size_t curry, const size_t currz,
    const boost::multi_array<double, 3> &U, double &a, double &b, double &c)
  {
    double minx, miny, minz;

    minx = 1000;
    if (currx > 0)
      {
        minx = std::min(minx, U[currx - 1][curry][currz]);
      }
    if (currx < U.shape()[0] - 1)
      {
        minx = std::min(minx, U[currx + 1][curry][currz]);
      }

    miny = 1000;
    if (curry > 0)
      {
        miny = std::min(miny, U[currx][curry - 1][currz]);
      }
    if (curry < U.shape()[1] - 1)
      {
        miny = std::min(miny, U[currx][curry + 1][currz]);
      }

    minz = 1000;
    if (currz > 0)
      {
        minz = std::min(minz, U[currx][curry][currz - 1]);
      }
    if (currz < U.shape()[2] - 1)
      {
        minz = std::min(minz, U[currx][curry][currz + 1]);
      }

    std::vector<double> vals(3, 0.0);
    vals[0] = minx;
    vals[1] = miny;
    vals[2] = minz;
    std::sort(vals.begin(), vals.end());
    a = vals[2];
    b = vals[1];
    c = vals[0];
    /*c = std::min(minx, std::min(miny, minz));
     a = std::max(minx, std::max(miny, minz));
     if (minx != a && minx != c)
     b = minx;
     if (miny != a && miny != c)
     b = miny;
     if (minz != a && minz != c)
     b = minz;*/

  }

void UpdateOne(size_t i, size_t j, size_t k, boost::multi_array<double, 3> &U,
    std::vector<bool> &L, std::vector<bool> &Lnew, const boost::multi_array<
        double, 3> &Vel)
  {
    size_t offset = U.shape()[2] * (i * U.shape()[1] + j) + k;
    if (!L[offset])
      {
        double p = U[i][j][k];
        double a, b, c;
        FindMinDir(i, j, k, U, a, b, c);
        double q = SolveU(a, b, c, Vel[i][j][k]);
        if (p > q)
          {
            U[i][j][k] = q;
            std::copy(U.origin(), U.origin() + U.num_elements(),
                std::ostream_iterator<double>(std::cout, " "));
            std::cout << std::endl;

            Lnew[offset] = true;

          }
      }
  }

void UpdateNeighbors(size_t i, size_t j, size_t k,
    boost::multi_array<double, 3> &U, std::vector<bool> &L,
    std::vector<bool> &Lnew, const boost::multi_array<double, 3> &Vel)
  {
    if (i < U.shape()[0] - 1)
      {
        UpdateOne(i + 1, j, k, U, L, Lnew, Vel);
      }
    if (i > 0)
      {
        UpdateOne(i - 1, j, k, U, L, Lnew, Vel);
      }
    if (j < U.shape()[1] - 1)
      {
        UpdateOne(i, j + 1, k, U, L, Lnew, Vel);
      }
    if (j > 0)
      {
        UpdateOne(i, j - 1, k, U, L, Lnew, Vel);
      }
    if (k < U.shape()[2] - 1)
      {
        UpdateOne(i, j, k + 1, U, L, Lnew, Vel);
      }
    if (k > 0)
      {
        UpdateOne(i, j, k - 1, U, L, Lnew, Vel);
      }
  }

int main()
  {
    const size_t nx = 10;
    const size_t ny = 10;
    const size_t nz = 5;
    const size_t ngrid = nx * ny * nz;
    boost::multi_array<double, 3> U(boost::extents[nx][ny][nz]);
    jiba::ThreeDSeismicModel Model;
    Model.SetSlownesses().resize(boost::extents[nx][ny][nz]);

    const double sourcex_index = 5;
    const double sourcey_index = 5;
    const double sourcez_index = 2;
    const double deltax = 5;
    const double InfValue = 1000;
    const double eps = 1e-10;
    std::fill_n(U.origin(), ngrid, InfValue);
    Model.SetCellSize(deltax, nx, ny, nz);
    std::fill_n(Model.SetSlownesses().origin(), ngrid, 1.0);
    std::vector<bool> L(ngrid, false);

    U[sourcex_index][sourcey_index][sourcez_index] = 0.0;
    L[Model.IndexToOffset(sourcex_index + 1, sourcey_index, sourcez_index)]
        = true;
    L[Model.IndexToOffset(sourcex_index - 1, sourcey_index, sourcez_index)]
        = true;
    L[Model.IndexToOffset(sourcex_index, sourcey_index + 1, sourcez_index)]
        = true;
    L[Model.IndexToOffset(sourcex_index, sourcey_index - 1, sourcez_index)]
        = true;
    L[Model.IndexToOffset(sourcex_index, sourcey_index, sourcez_index + 1)]
        = true;
    L[Model.IndexToOffset(sourcex_index, sourcey_index, sourcez_index - 1)]
        = true;
    size_t nactive = 6;
    size_t iteration = 1;
    do
      {
        std::vector<bool> Lnew(L);
        for (int index = 0; index < ngrid; ++index)
          {
            if (L[index])
              {
                int i, j, k;
                Model.OffsetToIndex(index, i, j, k);
                double p = U[i][j][k];
                double a, b, c;
                FindMinDir(i, j, k, U, a, b, c);
                double q = SolveU(a, b, c, Model.GetSlownesses()[i][j][k]);
                U[i][j][k] = q;
                if (fabs(p - q) < eps)
                  {
                    UpdateNeighbors(i, j, k, U, L, Lnew, Model.GetSlownesses());

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
            &Active[0], deltax, nx, ny, nz);
        iteration++;
        std::copy(U.origin(), U.origin() + U.num_elements(),
            std::ostream_iterator<double>(std::cout, " "));
        std::cout << std::endl;
      } while (nactive > 0);
    jiba::PlotTimeField("times.vtk", U.origin(), deltax, nx, ny, nz);
     std::cout << "U[8][8][4]: " <<  U[8][8][4] << std::endl;
     std::cout << "U[8][5][2]: " <<  U[8][sourcey_index][sourcez_index] << std::endl;
     /*std::cout << U[11][10][10] << std::endl;
     std::cout << U[10][11][10] << std::endl;
     std::cout << U[10][10][11] << std::endl;
     std::cout << U[11][11][10] << std::endl;
     std::cout << U[11][10][11] << std::endl;
     std::cout << U[10][11][11] << std::endl;
     std::cout << U[11][11][11] << std::endl;*/
  }
