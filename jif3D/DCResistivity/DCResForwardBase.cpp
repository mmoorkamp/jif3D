//============================================================================
// Name        : DCResForwardBase.cpp
// Author      : Zhanjie Shi and Richard.W Hobbs
// Version     : April 2014
// Copyright   : 2014, Zhanjie Shi and Richard.W Hobbs
//============================================================================

#include "../Global/FatalException.h"
#include "../Global/convert.h"
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/LU>
#include <Eigen/Core>
#include <boost/config.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp>
#include <boost/multi_array.hpp>
#include <omp.h>
#include <limits>
#include "DCResForwardBase.h"

using namespace Eigen;

namespace jif3D
  {
    /****** Interpolation function used to interpolate source and receiver position to cell centre
     *@param x  The centre coordinate of each cell in x direction.
     *@param y  The centre coordinate of each cell in y direction.
     *@param z  The centre coordinate of each cell in z direction.
     *@param xr The source/receiver coordinate in x direction.
     *@param yr The source/receiver coordinate in y direction.
     *@param zr The source/receiver coordinate in z direction.
     * */
    std::vector<double> Linint(const std::vector<double> &x, const std::vector<double> &y,
        const std::vector<double> &z, double xr, double yr, double zr)
      {

        double minx;
        size_t imx = 0;
        std::vector<size_t> ind_x(2, 0);
        std::vector<double> dx(2, 0.0);
        minx = fabs(xr - x[0]);
        for (size_t i = 1; i < x.size(); i++)
          {
            if (fabs(xr - x[i]) < minx)
              {
                minx = fabs(xr - x[i]);
                imx = i;
              }
          }
        if ((xr - x[imx]) >= 0.0)
          {
            ind_x[0] = imx;
            ind_x[1] = imx + 1;
          }
        else if ((xr - x[imx]) < 0.0)
          {
            ind_x[0] = imx - 1;
            ind_x[1] = imx;
          }
        dx[0] = xr - x[ind_x[0]];
        dx[1] = x[ind_x[1]] - xr;

        double miny;
        size_t imy = 0;
        std::vector<size_t> ind_y(2, 0);
        std::vector<double> dy(2, 0.0);
        miny = fabs(yr - y[0]);
        for (size_t i = 1; i < y.size(); i++)
          {
            if (fabs(yr - y[i]) < miny)
              {
                miny = fabs(yr - y[i]);
                imy = i;
              }
          }
        if ((yr - y[imy]) >= 0.0)
          {
            ind_y[0] = imy;
            ind_y[1] = imy + 1;
          }
        else if ((yr - y[imy]) < 0.0)
          {
            ind_y[0] = imy - 1;
            ind_y[1] = imy;
          }
        dy[0] = yr - y[ind_y[0]];
        dy[1] = y[ind_y[1]] - yr;

        double minz;
        size_t imz = 0;
        std::vector<size_t> ind_z(2, 0);
        std::vector<double> dz(2, 0.0);
        minz = fabs(zr - z[0]);
        for (size_t i = 1; i < z.size(); i++)
          {
            if (fabs(zr - z[i]) < minz)
              {
                minz = fabs(zr - z[i]);
                imz = i;
              }
          }
        if ((zr - z[imz]) >= 0.0)
          {
            ind_z[0] = imz;
            ind_z[1] = imz + 1;
          }
        else if ((zr - z[imz]) < 0.0)
          {
            ind_z[0] = imz - 1;
            ind_z[1] = imz;
          }
        dz[0] = zr - z[ind_z[0]];
        dz[1] = z[ind_z[1]] - zr;

        double Dx = (x[ind_x[1]] - x[ind_x[0]]);
        double Dy = (y[ind_y[1]] - y[ind_y[0]]);
        double Dz = (z[ind_z[1]] - z[ind_z[0]]);

        std::vector<double> v(x.size() * y.size() * z.size(), 0);

        v[ind_x[0] + ind_y[0] * x.size() + ind_z[0] * x.size() * y.size()] = (1
            - dx[0] / Dx) * (1 - dy[0] / Dy) * (1 - dz[0] / Dz);
        v[ind_x[0] + ind_y[1] * x.size() + ind_z[0] * x.size() * y.size()] = (1
            - dx[0] / Dx) * (1 - dy[1] / Dy) * (1 - dz[0] / Dz);
        v[ind_x[1] + ind_y[0] * x.size() + ind_z[0] * x.size() * y.size()] = (1
            - dx[1] / Dx) * (1 - dy[0] / Dy) * (1 - dz[0] / Dz);
        v[ind_x[1] + ind_y[1] * x.size() + ind_z[0] * x.size() * y.size()] = (1
            - dx[1] / Dx) * (1 - dy[1] / Dy) * (1 - dz[0] / Dz);
        v[ind_x[0] + ind_y[0] * x.size() + ind_z[1] * x.size() * y.size()] = (1
            - dx[0] / Dx) * (1 - dy[0] / Dy) * (1 - dz[1] / Dz);
        v[ind_x[0] + ind_y[1] * x.size() + ind_z[1] * x.size() * y.size()] = (1
            - dx[0] / Dx) * (1 - dy[1] / Dy) * (1 - dz[1] / Dz);
        v[ind_x[1] + ind_y[0] * x.size() + ind_z[1] * x.size() * y.size()] = (1
            - dx[1] / Dx) * (1 - dy[0] / Dy) * (1 - dz[1] / Dz);
        v[ind_x[1] + ind_y[1] * x.size() + ind_z[1] * x.size() * y.size()] = (1
            - dx[1] / Dx) * (1 - dy[1] / Dy) * (1 - dz[1] / Dz);

        return v;

      }

    /****** Basic forward modelling function */
    jif3D::rvec ResForward(const GEOMETRY_RES &geo, const GRID_STRUCT_RES &grid,
        size_t ndata_res)
      {
        std::vector<double> cellwidth_x(grid.nx, grid.dx);
        std::vector<double> cellwidth_y(grid.ny, grid.dy);
        std::vector<double> cellwidth_z(grid.dz);

        //*******Generate index and value of non-zero elements of the divergence matrix D which is a sparse matrix.
        std::vector<size_t> lxd, lyd, lzd, jxd, jyd, jzd; //the index of non-zero elements in sparse matrix D.
        std::vector<double> kxd, kyd, kzd; //the value of non-zero elements in sparse matrix D.
        /*****Divergence matrix D[np][nax+nay+naz], np=grid.nx*grid.ny*grid.nz, nax=(grid.nx-1)*grid.ny*grid.nz,
         * nay=grid.nx*(grid.ny-1)*grid.nz, naz=grid.nx*grid.ny*(grid.nz-1).
         *for (size_t i=0; i<2*nax; i++), D[lxd[i]][jxd[i]]=kxd[i]
         *for (size_t i=0; i<2*nay; i++), D[lyd[i]][jyd[i]+nax]=kyd[i]
         *for (size_t i=0; i<2*naz; i++), D[lzd[i]][jzd[i]+nax+nay]=kzd[i]
         *****/

        //***generate d/dx
        //*Entries(l,j,k)
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 1; k < grid.nx - 1; k++)
              {
                lxd.push_back(k + j * grid.nx + i * grid.nx * grid.ny);
              }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 0; k < grid.nx - 2; k++)
              {
                jxd.push_back(k + j * (grid.nx - 1) + i * (grid.nx - 1) * grid.ny);
              }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 1; k < grid.nx - 1; k++)
              {
                kxd.push_back(-1.0 / (cellwidth_x[k] * 1.0 * 1.0));
              }
        //*Entries(l+1,j,k)
        for (size_t i = 0; i < (grid.nx - 2) * grid.ny * grid.nz; i++)
          {
            lxd.push_back(lxd[i]);
          }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 1; k < grid.nx - 1; k++)
              {
                jxd.push_back(k + j * (grid.nx - 1) + i * (grid.nx - 1) * grid.ny);
              }
        for (size_t i = 0; i < (grid.nx - 2) * grid.ny * grid.nz; i++)
          {
            kxd.push_back(-kxd[i]);
          }
        //*BC at x=0
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            {
              lxd.push_back(j * grid.nx + i * grid.nx * grid.ny);
            }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            {
              jxd.push_back(j * (grid.nx - 1) + i * (grid.nx - 1) * grid.ny);
            }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            {
              kxd.push_back(1.0 / (cellwidth_x[0] * 1.0 * 1.0));
            }
        //*BC at x=end
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            {
              lxd.push_back(grid.nx - 1 + j * grid.nx + i * grid.nx * grid.ny);
            }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            {
              jxd.push_back(
                  grid.nx - 2 + j * (grid.nx - 1) + i * (grid.nx - 1) * grid.ny);
            }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            {
              kxd.push_back(-1.0 / (cellwidth_x[grid.nx - 1] * 1.0 * 1.0));
            }
        //***generate d/dy
        //*Entries(l,j,k)
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 1; j < grid.ny - 1; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                lyd.push_back(k + j * grid.nx + i * grid.nx * grid.ny);
              }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny - 2; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                jyd.push_back(k + j * grid.nx + i * grid.nx * (grid.ny - 1));
              }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 1; j < grid.ny - 1; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                kyd.push_back(-1.0 / (1.0 * cellwidth_y[j] * 1.0));
              }
        //*Entries(l+1,j,k)
        for (size_t i = 0; i < grid.nx * (grid.ny - 2) * grid.nz; i++)
          {
            lyd.push_back(lyd[i]);
          }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 1; j < grid.ny - 1; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                jyd.push_back(k + j * grid.nx + i * grid.nx * (grid.ny - 1));
              }
        for (size_t i = 0; i < grid.nx * (grid.ny - 2) * grid.nz; i++)
          {
            kyd.push_back(-kyd[i]);
          }
        //*BC on y=0
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t k = 0; k < grid.nx; k++)
            {
              lyd.push_back(k + i * grid.nx * grid.ny);
            }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t k = 0; k < grid.nx; k++)
            {
              jyd.push_back(k + i * grid.nx * (grid.ny - 1));
            }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t k = 0; k < grid.nx; k++)
            {
              kyd.push_back(1.0 / (1.0 * cellwidth_y[0] * 1.0));
            }
        //*BC on y=end
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t k = 0; k < grid.nx; k++)
            {
              lyd.push_back(k + grid.nx * (grid.ny - 1) + i * grid.nx * grid.ny);
            }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t k = 0; k < grid.nx; k++)
            {
              jyd.push_back(k + grid.nx * (grid.ny - 2) + i * grid.nx * (grid.ny - 1));
            }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t k = 0; k < grid.nx; k++)
            {
              kyd.push_back(-1.0 / (1.0 * cellwidth_y[grid.ny - 1] * 1.0));
            }
        //***generate d/dz
        //*Entries(l,j,k)
        for (size_t i = 1; i < grid.nz - 1; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                lzd.push_back(k + j * grid.nx + i * grid.nx * grid.ny);
              }
        for (size_t i = 0; i < grid.nz - 2; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                jzd.push_back(k + j * grid.nx + i * grid.nx * grid.ny);
              }
        for (size_t i = 1; i < grid.nz - 1; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                kzd.push_back(-1.0 / (1.0 * 1.0 * cellwidth_z[i]));
              }
        //*Entries(l+1,j,k)
        for (size_t i = 0; i < grid.nx * grid.ny * (grid.nz - 2); i++)
          {
            lzd.push_back(lzd[i]);
          }
        for (size_t i = 1; i < grid.nz - 1; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                jzd.push_back(k + j * grid.nx + i * grid.nx * grid.ny);
              }
        for (size_t i = 0; i < grid.nx * grid.ny * (grid.nz - 2); i++)
          {
            kzd.push_back(-kzd[i]);
          }
        //*BC on z=0
        for (size_t j = 0; j < grid.ny; j++)
          for (size_t k = 0; k < grid.nx; k++)
            {
              lzd.push_back(k + j * grid.nx);
            }
        for (size_t j = 0; j < grid.ny; j++)
          for (size_t k = 0; k < grid.nx; k++)
            {
              jzd.push_back(k + j * grid.nx);
            }
        for (size_t j = 0; j < grid.ny; j++)
          for (size_t k = 0; k < grid.nx; k++)
            {
              kzd.push_back(1.0 / (1.0 * 1.0 * cellwidth_z[0]));
            }
        //*BC on z=end
        for (size_t j = 0; j < grid.ny; j++)
          for (size_t k = 0; k < grid.nx; k++)
            {
              lzd.push_back(k + j * grid.nx + grid.nx * grid.ny * (grid.nz - 1));
            }
        for (size_t j = 0; j < grid.ny; j++)
          for (size_t k = 0; k < grid.nx; k++)
            {
              jzd.push_back(k + j * grid.nx + grid.nx * grid.ny * (grid.nz - 2));
            }
        for (size_t j = 0; j < grid.ny; j++)
          for (size_t k = 0; k < grid.nx; k++)
            {
              kzd.push_back(-1.0 / (1.0 * 1.0 * cellwidth_z[grid.nz - 1]));
            }

        //*******Generate index and value of non-zero elements of the gradient matrix G which is a sparse matrix.
        std::vector<size_t> lxg, lyg, lzg, jxg, jyg, jzg; //the index of non-zero elements in sparse matrix G.
        std::vector<double> kxg, kyg, kzg; //the value of non-zero elements in sparse matrix G.

        /*****Gradient matrix G[nax+nay+naz][np], nax=(grid.nx-1)*grid.ny*grid.nz,
         * nay=grid.nx*(grid.ny-1)*grid.nz, naz=grid.nx*grid.ny*(grid.nz-1),np=grid.nx*grid.ny*grid.nz.
         *for (size_t i=0; i<2*nax; i++), G[lxg[i]][jxg[i]]=kxg[i];
         *for (size_t i=0; i<2*nay; i++), G[lyg[i]+nax][jyg[i]]=kyg[i]
         *for (size_t i=0; i<2*naz; i++), G[lzg[i]+nax+nay][jzg[i]]=kzg[i]
         *****/

        //***generate d/dx
        //*Entries (l,j,k)
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 0; k < grid.nx - 1; k++)
              {
                lxg.push_back(k + j * (grid.nx - 1) + i * (grid.nx - 1) * grid.ny);
              }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 0; k < grid.nx - 1; k++)
              {
                jxg.push_back(k + j * grid.nx + i * grid.nx * grid.ny);
              }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 0; k < grid.nx - 1; k++)
              {
                kxg.push_back(
                    -2.0 / (cellwidth_x[k] * 1.0 * 1.0 + cellwidth_x[k + 1] * 1.0 * 1.0));
              }
        //*Entries (l+1,j,k)
        for (size_t i = 0; i < (grid.nx - 1) * grid.ny * grid.nz; i++)
          {
            lxg.push_back(lxg[i]);
          }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 1; k < grid.nx; k++)
              {
                jxg.push_back(k + j * grid.nx + i * grid.nx * grid.ny);
              }
        for (size_t i = 0; i < (grid.nx - 1) * grid.ny * grid.nz; i++)
          {
            kxg.push_back(-kxg[i]);
          }
        //***generate d/dy
        //*Entries (l,j,k)
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny - 1; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                lyg.push_back(k + j * grid.nx + i * grid.nx * (grid.ny - 1));
              }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny - 1; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                jyg.push_back(k + j * grid.nx + i * grid.nx * grid.ny);
              }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny - 1; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                kyg.push_back(
                    -2.0 / (1.0 * cellwidth_y[j] * 1.0 + 1.0 * cellwidth_y[j + 1] * 1.0));
              }
        //*Entries (l+1,j,k)
        for (size_t i = 0; i < grid.nx * (grid.ny - 1) * grid.nz; i++)
          {
            lyg.push_back(lyg[i]);
          }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 1; j < grid.ny; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                jyg.push_back(k + j * grid.nx + i * grid.nx * grid.ny);
              }
        for (size_t i = 0; i < grid.nx * (grid.ny - 1) * grid.nz; i++)
          {
            kyg.push_back(-kyg[i]);
          }
        //***generate d/dz
        //*Entries (l,j,k)
        for (size_t i = 0; i < grid.nz - 1; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                lzg.push_back(k + j * grid.nx + i * grid.nx * grid.ny);
              }
        for (size_t i = 0; i < grid.nz - 1; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                jzg.push_back(k + j * grid.nx + i * grid.nx * grid.ny);
              }
        for (size_t i = 0; i < grid.nz - 1; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                kzg.push_back(
                    -2.0 / (1.0 * 1.0 * cellwidth_z[i] + 1.0 * 1.0 * cellwidth_z[i + 1]));
              }
        //*Entries (l+1,j,k)
        for (size_t i = 0; i < grid.nx * grid.ny * (grid.nz - 1); i++)
          {
            lzg.push_back(lzg[i]);
          }
        for (size_t i = 1; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                jzg.push_back(k + j * grid.nx + i * grid.nx * grid.ny);
              }
        for (size_t i = 0; i < grid.nx * grid.ny * (grid.nz - 1); i++)
          {
            kzg.push_back(-kzg[i]);
          }

        //*******Generate index and value of non-zero elements of the diagonal matrix S(rho) which is a sparse matrix.
        std::vector<size_t> lxs, lys, lzs, jxs, jys, jzs; //the index of non-zero elements in sparse matrix S.
        std::vector<double> kxs, kys, kzs; //the value of non-zero elements in sparse matrix S.

        /*****Diagonal matrix S[nax+nay+naz][nax+nay+naz], nax=(grid.nx-1)*grid.ny*grid.nz,
         * nay=grid.nx*(grid.ny-1)*grid.nz, naz=grid.nx*grid.ny*(grid.nz-1).
         *for (size_t i=0; i<nax; i++), S[lxs[i]][jxs[i]]=kxs[i];
         *for (size_t i=0; i<nay; i++), S[lys[i]+nax][jys[i]+nax]=kys[i]
         *for (size_t i=0; i<naz; i++), S[lzs[i]+nax+nay][jzs[i]+nax+nay]=kzs[i]
         *****/

        std::vector<double> rhof_x, rhof_y, rhof_z;
        //***generate x coefficients
        //*Average rho on x face
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 1; k < grid.nx; k++)
              {
                rhof_x.push_back(
                    (cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[i]
                        * grid.rho[k + j * grid.nx + i * grid.nx * grid.ny]
                        + cellwidth_x[k - 1] * cellwidth_y[j] * cellwidth_z[i]
                            * grid.rho[k - 1 + j * grid.nx + i * grid.nx * grid.ny])
                        / (cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[i]
                            + cellwidth_x[k - 1] * cellwidth_y[j] * cellwidth_z[i]));
              }
        for (size_t i = 0; i < (grid.nx - 1) * grid.ny * grid.nz; i++)
          {
            lxs.push_back(i);
            jxs.push_back(i);
            kxs.push_back(rhof_x[i]);
          }
        //***generate y coefficients
        //*Average rho on y face
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 1; j < grid.ny; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                rhof_y.push_back(
                    (cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[i]
                        * grid.rho[k + j * grid.nx + i * grid.nx * grid.ny]
                        + cellwidth_x[k] * cellwidth_y[j - 1] * cellwidth_z[i]
                            * grid.rho[k + (j - 1) * grid.nx + i * grid.nx * grid.ny])
                        / (cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[i]
                            + cellwidth_x[k] * cellwidth_y[j - 1] * cellwidth_z[i]));
              }
        for (size_t i = 0; i < grid.nx * (grid.ny - 1) * grid.nz; i++)
          {
            lys.push_back(i);
            jys.push_back(i);
            kys.push_back(rhof_y[i]);
          }
        //***generate z coefficients
        //*Average rho on z face
        for (size_t i = 1; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                rhof_z.push_back(
                    (cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[i]
                        * grid.rho[k + j * grid.nx + i * grid.nx * grid.ny]
                        + cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[i - 1]
                            * grid.rho[k + j * grid.nx + (i - 1) * grid.nx * grid.ny])
                        / (cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[i]
                            + cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[i - 1]));
              }
        for (size_t i = 0; i < grid.nx * grid.ny * (grid.nz - 1); i++)
          {
            lzs.push_back(i);
            jzs.push_back(i);
            kzs.push_back(rhof_z[i]);
          }

        //*******Generate Forward operator A=D*S*G.
        size_t np = grid.nx * grid.ny * grid.nz;
        size_t nax = (grid.nx - 1) * grid.ny * grid.nz;
        size_t nay = grid.nx * (grid.ny - 1) * grid.nz;
        size_t naz = grid.nx * grid.ny * (grid.nz - 1);

        //generate sparse matrix D,
        typedef Eigen::Triplet<double> DT;
        std::vector<DT> DV;
        DV.reserve(2 * (nax + nay + naz));
        for (size_t i = 0; i < 2 * nax; i++)
          {
            DV.push_back(DT(lxd[i], jxd[i], kxd[i]));
          }
        for (size_t i = 0; i < 2 * nay; i++)
          {
            DV.push_back(DT(lyd[i], jyd[i] + nax, kyd[i]));
          }
        for (size_t i = 0; i < 2 * naz; i++)
          {
            DV.push_back(DT(lzd[i], jzd[i] + nax + nay, kzd[i]));
          }
        Eigen::SparseMatrix<double> D(np, nax + nay + naz);
        D.setFromTriplets(DV.begin(), DV.end());

        //generate sparse matrix S,
        typedef Eigen::Triplet<double> ST;
        std::vector<ST> SV;
        SV.reserve(nax + nay + naz);
        for (size_t i = 0; i < nax; i++)
          {
            SV.push_back(ST(lxs[i], jxs[i], 1.0 / kxs[i]));
          }
        for (size_t i = 0; i < nay; i++)
          {
            SV.push_back(ST(lys[i] + nax, jys[i] + nax, 1.0 / kys[i]));
          }
        for (size_t i = 0; i < naz; i++)
          {
            SV.push_back(ST(lzs[i] + nax + nay, jzs[i] + nax + nay, 1.0 / kzs[i]));
          }
        Eigen::SparseMatrix<double> S(nax + nay + naz, nax + nay + naz);
        S.setFromTriplets(SV.begin(), SV.end());

        //generate sparse matrix G,
        typedef Eigen::Triplet<double> GT;
        std::vector<GT> GV;
        GV.reserve(2 * (nax + nay + naz));
        for (size_t i = 0; i < 2 * nax; i++)
          {
            GV.push_back(GT(lxg[i], jxg[i], kxg[i]));
          }
        for (size_t i = 0; i < 2 * nay; i++)
          {
            GV.push_back(GT(lyg[i] + nax, jyg[i], kyg[i]));
          }
        for (size_t i = 0; i < 2 * naz; i++)
          {
            GV.push_back(GT(lzg[i] + nax + nay, jzg[i], kzg[i]));
          }
        Eigen::SparseMatrix<double> G(nax + nay + naz, np);
        G.setFromTriplets(GV.begin(), GV.end());

        //generate sparse matrix, forward operator A
        Eigen::SparseMatrix<double> A1(np, nax + nay + naz), A(np, np);
        A1 = D * S;
        A = A1 * G;
        A.coeffRef(0, 0) = 1.0 / (cellwidth_x[0] * cellwidth_y[0] * cellwidth_z[0]);

        //*****Permutation of forward operator A using Cuthill-Mckee Ordering Method*********************************//
        /* Because Eigen library has not Cuthill-Mckee Ordering function, but Boost library has, so we used Boost uBlas library to finish Cuthill-Mckee Ordering of forward operator A.
         * When implementing Cuthill-Mckee Ordering using Boost, the forward operator A which is Eigen SparseMatrix need to be input into a Boost SparseMatrix firstly. We use a loop to
         * input the non-zero element of A one by one, so this probably can cause program running slower. But it is difficult to generate forward operator A using Boost C++ library due
         * to no proper multiplication function of sparse matrix can be used. So we use Eigen to generate forward operator A from D*S*G as above. In the future, if there is a new function
         *  which can be used to implement sparse matrix product in Boost uBlas library, it is a better method to directly use Boost to generate forward operator A from D*S*G.
         */
        typedef boost::numeric::ublas::compressed_matrix<double> MatType;
        typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
            boost::property<boost::vertex_color_t, boost::default_color_type,
                boost::property<boost::vertex_degree_t, size_t> > > GraphType;
        typedef boost::graph_traits<GraphType>::vertex_descriptor VertexType;
        typedef boost::graph_traits<GraphType>::vertices_size_type size_type;
        GraphType GFOR;
        //declare a Boost sparse matrix FOR, which is equal to forward operator A which is Eigen Sparse Matrix.
        MatType FOR(np, np);
        //order is corresponding to 'p' in Matlab and is used to rearrange sequence of FOR for permutation, inv_order is the same with 'up' in Matlab and is used to regain original sequence.
        boost::numeric::ublas::vector<size_t> order(np), inv_order(order.size());
        //declare row, col number and value of non-zero element of permutation of Sparse Matrix FOR.
        std::vector<size_t> rowpFOR, colpFOR;
        std::vector<double> valuepFOR;
        //Input non-zero element of Eigen Sparse Matrix A into Boost sparse matrix FOR.
        for (int k = 0; k < A.outerSize(); ++k)
          {
            for (SparseMatrix<double>::InnerIterator it(A, k); it; ++it)
              {
                FOR(it.row(), it.col()) = it.value();
              }
          }
        //Input FOR into GFOR, GFOR is a Graph Class. In Boost uBlas C++ library, it is using Graph theory to implement Cuthill-Mckee Ordering of Sparse Matrix.
        for (MatType::const_iterator1 row_it = FOR.begin1(); row_it != FOR.end1();
            ++row_it)
          {
            for (MatType::const_iterator2 col_it = row_it.begin(); col_it != row_it.end();
                ++col_it)
              {
                boost::add_edge(col_it.index1(), col_it.index2(), GFOR);
              }
          }
        // Get order
        boost::property_map<GraphType, boost::vertex_index_t>::type index_map = get(
            boost::vertex_index, GFOR);
        // Use the boost::graph cuthill mckee algorithm
        std::vector<VertexType> inv_perm(boost::num_vertices(GFOR));
        cuthill_mckee_ordering(GFOR, inv_perm.rbegin(), get(boost::vertex_color, GFOR),
            make_degree_map(GFOR));
        for (size_t i = 0; i < inv_perm.size(); i++)
          {
            order[i] = index_map[inv_perm[i]];
          }
        // Get inv_order
        for (size_t i = 0; i < order.size(); i++)
          {
            inv_order[order[i]] = i;
          }
        //Get row,col and value of non-zero of permutation Sparse Matrix of FOR.
        for (MatType::const_iterator1 row_it = FOR.begin1(); row_it != FOR.end1();
            ++row_it)
          {
            for (MatType::const_iterator2 col_it = row_it.begin(); col_it != row_it.end();
                ++col_it)
              {
                size_t row_ind = inv_order(col_it.index1());
                size_t col_ind = inv_order(col_it.index2());
                rowpFOR.push_back(row_ind);
                colpFOR.push_back(col_ind);
                valuepFOR.push_back(*col_it);
              }
          }
        // Then generate permutation sparse matrix of FOR, pFOR and implement Incomplete LU (ILU) of pFOR.
        //generate permutation sparse matrix of FOR, pFOR
        typedef Eigen::Triplet<double> pFORT;
        std::vector<pFORT> pFORV;
        pFORV.reserve(valuepFOR.size());
        for (size_t i = 0; i < valuepFOR.size(); i++)
          {
            pFORV.push_back(pFORT(rowpFOR[i], colpFOR[i], valuepFOR[i]));
          }
        Eigen::SparseMatrix<double> pFOR(np, np);
        pFOR.setFromTriplets(pFORV.begin(), pFORV.end());
        //ILU of pFOR
        Eigen::BiCGSTAB<SparseMatrix<double>, Eigen::IncompleteLUT<double> > pFORBCGST;
        pFORBCGST.preconditioner().setDroptol(0.00001);
        pFORBCGST.compute(pFOR);
        //A permutation of the forward operator has been finished.




        //********************************Boundary correction and source term improvement*********************************//
        /* In order to improve the BCs and have an accurate solution near a source electrode, we apply a correction to the source
         * term q that both improves the BCs and reduces errors near the source locations. The method is as follows,
         * Firstly, define a homogeneous conductivity field Con_H that is the best initial guess of the average conductivity in the model space.
         * Secondly, calculate the analytical solution for the half-space potential field U_H resulting from the conductivity structure Con_H from a source term q.
         * Note, the analytical solution considers a half space that goes to infinity and therefore has the "correct" boundary conditions for our problem.
         * Then, construct the forward operator A(Con_H).
         * Evaluate q_corr = A(Con_H)*U_H.
         * Once we have the vector q_corr, we substitute it for the source term q and solve the system,
         * U_true = A(Con_true)-1*q_corr.
         * Where U_true is the potential field everywhere in our model space and Con_true is the true conductivity field.
         *
         */
        /*Because U_H is calculated for each source pair, we implement it in the loop of calculation potential for true model.
         * Now, we firstly construct the forward operator A(Con_H).
         */
        //*******Generate index and value of non-zero elements of the diagonal matrix S(1/avg_cond) which is a sparse matrix.
        //***We add "_BC" after the similar variable name with that used during constructing Sparse matrix "S" in the above code.
        std::vector<size_t> lxs_BC, lys_BC, lzs_BC, jxs_BC, jys_BC, jzs_BC; //the index of non-zero elements in sparse matrix S(1/avg_cond).
        std::vector<double> kxs_BC, kys_BC, kzs_BC; //the value of non-zero elements in sparse matrix S(1/avg_cond).

        /*****Diagonal matrix S[nax+nay+naz][nax+nay+naz], nax=(grid.nx-1)*grid.ny*grid.nz,
         * nay=grid.nx*(grid.ny-1)*grid.nz, naz=grid.nx*grid.ny*(grid.nz-1).
         *for (size_t i=0; i<nax; i++), S[lxs[i]][jxs[i]]=kxs[i];
         *for (size_t i=0; i<nay; i++), S[lys[i]+nax][jys[i]+nax]=kys[i]
         *for (size_t i=0; i<naz; i++), S[lzs[i]+nax+nay][jzs[i]+nax+nay]=kzs[i]
         *****/

        std::vector<double> rhof_x_BC, rhof_y_BC, rhof_z_BC;
        //***generate x coefficients
        //*Average rho on x face
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 1; k < grid.nx; k++)
              {
                rhof_x_BC.push_back(
                    (cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[i]
                        * (1.0/grid.avg_cond)
                        + cellwidth_x[k - 1] * cellwidth_y[j] * cellwidth_z[i]
                            * (1.0/grid.avg_cond))
                        / (cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[i]
                            + cellwidth_x[k - 1] * cellwidth_y[j] * cellwidth_z[i]));
              }
        for (size_t i = 0; i < (grid.nx - 1) * grid.ny * grid.nz; i++)
          {
            lxs_BC.push_back(i);
            jxs_BC.push_back(i);
            kxs_BC.push_back(rhof_x_BC[i]);
          }
        //***generate y coefficients
        //*Average rho on y face
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 1; j < grid.ny; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                rhof_y_BC.push_back(
                    (cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[i]
                        * (1.0/grid.avg_cond)
                        + cellwidth_x[k] * cellwidth_y[j - 1] * cellwidth_z[i]
                            * (1.0/grid.avg_cond))
                        / (cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[i]
                            + cellwidth_x[k] * cellwidth_y[j - 1] * cellwidth_z[i]));
              }
        for (size_t i = 0; i < grid.nx * (grid.ny - 1) * grid.nz; i++)
          {
            lys_BC.push_back(i);
            jys_BC.push_back(i);
            kys_BC.push_back(rhof_y_BC[i]);
          }
        //***generate z coefficients
        //*Average rho on z face
        for (size_t i = 1; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                rhof_z_BC.push_back(
                    (cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[i]
                        * (1.0/grid.avg_cond)
                        + cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[i - 1]
                            * (1.0/grid.avg_cond))
                        / (cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[i]
                            + cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[i - 1]));
              }
        for (size_t i = 0; i < grid.nx * grid.ny * (grid.nz - 1); i++)
          {
            lzs_BC.push_back(i);
            jzs_BC.push_back(i);
            kzs_BC.push_back(rhof_z_BC[i]);
          }

        //*******construct the forward operator A(Con_H) = D*S(1/avg_cond)*G.
        //generate sparse matrix S(1/avg_cond),
        typedef Eigen::Triplet<double> ST_BC;
        std::vector<ST_BC> SV_BC;
        SV_BC.reserve(nax + nay + naz);
        for (size_t i = 0; i < nax; i++)
          {
            SV_BC.push_back(ST_BC(lxs_BC[i], jxs_BC[i], 1.0 / kxs_BC[i]));
          }
        for (size_t i = 0; i < nay; i++)
          {
            SV_BC.push_back(ST_BC(lys_BC[i] + nax, jys_BC[i] + nax, 1.0 / kys_BC[i]));
          }
        for (size_t i = 0; i < naz; i++)
          {
            SV_BC.push_back(ST_BC(lzs_BC[i] + nax + nay, jzs_BC[i] + nax + nay, 1.0 / kzs_BC[i]));
          }
        Eigen::SparseMatrix<double> S_BC(nax + nay + naz, nax + nay + naz);
        S_BC.setFromTriplets(SV_BC.begin(), SV_BC.end());
        //generate sparse matrix, forward operator A(Con_H) = D*S(1/avg_cond)*G.
        Eigen::SparseMatrix<double> A1_BC(np, nax + nay + naz), A_BC(np, np);
        A1_BC = D * S_BC;
        A_BC = A1_BC * G;
        A_BC.coeffRef(0, 0) = 1.0 / (cellwidth_x[0] * cellwidth_y[0] * cellwidth_z[0]);





        //*****Generate source correction operator q_corr, substitute source term q using q_corr and calculate forward response data*********************************//
        /* Firstly, calculate each cell center's coordinate used in interpolation operation to interpolate source/receiver's into cell's centre.
         */
        std::vector<double> zpos(grid.nz + 1, 0), xpos(grid.nx + 1, 0), ypos(grid.ny + 1,
            0);
        std::vector<double> centerzeroxpos(grid.nx + 1, 0), centerzeroypos(grid.ny + 1,
            0);
        std::vector<double> cellcenterzpos(grid.nz, 0), cellcenterxpos(grid.nx, 0),
            cellcenterypos(grid.ny, 0);
        //build the 3d grid - numbered from 0 to maximum extent,
        for (size_t i = 0; i < grid.nz; i++)
          {
            zpos[i + 1] = zpos[i] + cellwidth_z[i];
          }
        for (size_t i = 0; i < grid.nx; i++)
          {
            xpos[i + 1] = xpos[i] + cellwidth_x[i];
          }
        for (size_t i = 0; i < grid.ny; i++)
          {
            ypos[i + 1] = ypos[i] + cellwidth_y[i];
          }
        //center the grid about zero. x and y's zero point is at the centre of 3d volume, z' zero at surface.
        for (size_t i = 0; i < grid.nx + 1; i++)
          {
            centerzeroxpos[i] = xpos[i] - xpos[grid.nx] / 2.0;
          }
        for (size_t i = 0; i < grid.ny + 1; i++)
          {
            centerzeroypos[i] = ypos[i] - ypos[grid.ny] / 2.0;
          }
        //set cellcenter's position,
        for (size_t i = 0; i < grid.nx; i++)
          {
            cellcenterxpos[i] = centerzeroxpos[i] + cellwidth_x[i] / 2.0;
          }
        for (size_t i = 0; i < grid.ny; i++)
          {
            cellcenterypos[i] = centerzeroypos[i] + cellwidth_y[i] / 2.0;
          }
        for (size_t i = 0; i < grid.nz; i++)
          {
            cellcenterzpos[i] = zpos[i] + cellwidth_z[i] / 2.0;
          }

        //calculate how many data points we will have
	size_t ndata = 0;
	std::vector<int> startindices;
	startindices.push_back(0);
	for (size_t i = 0; i < geo.nsource; i++)
	{
		for (size_t Qi = 0; Qi < ndata_res; Qi++)
		{
			if (i == geo.sno[Qi]) {
				ndata++;
			}
		}
		startindices.push_back(ndata);
	}

        /* The follows are the main loop for calculating forward response for a model and corresponding geometry of sources/receivers.
         */
        jif3D::rvec result(ndata);

#pragma omp parallel for default(shared)
        for (size_t i = 0; i < geo.nsource; i++)
          {
            std::vector<double> pvel(np, 0), nvel(np, 0), pveimagl(np, 0), nveimagl(np, 0), U_H(np, 0), q_corr(np, 0), Q1(np, 0), Q2(np, 0), Q;
            Eigen::VectorXd permutationeq(np), permutationeu, eu(np), eforwarddata, eU_H(np), eq_corr(np);

            //generate source term q,
            //q1 = Linint(cellcenterxpos, cellcenterypos, cellcenterzpos, geo.PosSx[i],
                //geo.PosSy[i], geo.PosSz[i]);
            //q2 = Linint(cellcenterxpos, cellcenterypos, cellcenterzpos, geo.NegSx[i],
                //geo.NegSy[i], geo.NegSz[i]);
            //q.reserve(np);
            //for (size_t j = 0; j < grid.nz; j++)
              //for (size_t k = 0; k < grid.ny; k++)
                //for (size_t l = 0; l < grid.nx; l++)
                  //{
                    //q.push_back(
                        //(q1[l + k * grid.nx + j * grid.nx * grid.ny]
                            //- q2[l + k * grid.nx + j * grid.nx * grid.ny])
                           /// (cellwidth_x[l] * cellwidth_y[k] * cellwidth_z[j]));
                  //}

            //norm of current electrode and 1st potential electrode and generate potential of homogeneous model, U_H
            for (size_t l=0; l<grid.nz; l++)
            	for (size_t m=0; m<grid.ny; m++)
            		for (size_t n=0; n<grid.nx; n++)
            		{
            			//norm of current electrode and 1st potential electrode
            			pvel[n+m*grid.nx+l*grid.nx*grid.ny] = std::sqrt((cellcenterxpos[n]-geo.PosSx[i])*(cellcenterxpos[n]-geo.PosSx[i])
            					+(cellcenterypos[m]-geo.PosSy[i])*(cellcenterypos[m]-geo.PosSy[i])+(cellcenterzpos[l]-geo.PosSz[i])*(cellcenterzpos[l]-geo.PosSz[i]));
            			nvel[n+m*grid.nx+l*grid.nx*grid.ny] = std::sqrt((cellcenterxpos[n]-geo.NegSx[i])*(cellcenterxpos[n]-geo.NegSx[i])
            					+(cellcenterypos[m]-geo.NegSy[i])*(cellcenterypos[m]-geo.NegSy[i])+(cellcenterzpos[l]-geo.NegSz[i])*(cellcenterzpos[l]-geo.NegSz[i]));
            			//norm of imaginary current electrode and 1st potential electrode
            			pveimagl[n+m*grid.nx+l*grid.nx*grid.ny] = std::sqrt((cellcenterxpos[n]-geo.PosSx[i])*(cellcenterxpos[n]-geo.PosSx[i])
            					+(cellcenterypos[m]-geo.PosSy[i])*(cellcenterypos[m]-geo.PosSy[i])+(cellcenterzpos[l]+geo.PosSz[i])*(cellcenterzpos[l]+geo.PosSz[i]));
            			nveimagl[n+m*grid.nx+l*grid.nx*grid.ny] = std::sqrt((cellcenterxpos[n]-geo.NegSx[i])*(cellcenterxpos[n]-geo.NegSx[i])
            					+(cellcenterypos[m]-geo.NegSy[i])*(cellcenterypos[m]-geo.NegSy[i])+(cellcenterzpos[l]+geo.NegSz[i])*(cellcenterzpos[l]+geo.NegSz[i]));
            			//generate potential of homogeneous model, U_H
            			U_H[n+m*grid.nx+l*grid.nx*grid.ny] = (1.0/(grid.avg_cond*4.0*3.1415926))*(1.0/pvel[n+m*grid.nx+l*grid.nx*grid.ny]-1.0/nvel[n+m*grid.nx+l*grid.nx*grid.ny]
            					                              +1.0/pveimagl[n+m*grid.nx+l*grid.nx*grid.ny]-1.0/nveimagl[n+m*grid.nx+l*grid.nx*grid.ny]);
            		}
            //now check for singularities due to the source being on a node
            double singularityvalue = std::numeric_limits<double>::infinity();
            for (size_t Uid=0; Uid<np; Uid++)
            {
            	if (U_H[Uid]==singularityvalue||U_H[Uid]==-singularityvalue)
            	{
            		size_t idx, idy, idz;
            		idx = Uid % grid.nx;
            		idz = (Uid - idx) / grid.nx;
            		idy = idz % grid.ny;
            		idz = (idz - idy) / grid.ny;
            		if (idz==0)
            			U_H[Uid] = (U_H[idx+1+idy*grid.nx+idz*grid.nx*grid.ny]+U_H[idx+(idy+1)*grid.nx+idz*grid.nx*grid.ny]
            			               +U_H[idx+idy*grid.nx+(idz+1)*grid.nx*grid.ny]+U_H[idx-1+idy*grid.nx+idz*grid.nx*grid.ny]
            			                   +U_H[idx+(idy-1)*grid.nx+idz*grid.nx*grid.ny])/5.0;
            		else
            			U_H[Uid] = (U_H[idx+1+idy*grid.nx+idz*grid.nx*grid.ny]+U_H[idx+(idy+1)*grid.nx+idz*grid.nx*grid.ny]
            			               +U_H[idx+idy*grid.nx+(idz+1)*grid.nx*grid.ny]+U_H[idx-1+idy*grid.nx+idz*grid.nx*grid.ny]
            			                   +U_H[idx+(idy-1)*grid.nx+idz*grid.nx*grid.ny]+U_H[idx+idy*grid.nx+(idz-1)*grid.nx*grid.ny])/6.0;
            	}
            }
            //calculate correction operator q_corr
            for (size_t cid=0; cid<np; cid++)
            {
            	eU_H(cid) = U_H[cid];
            }
            eq_corr=A_BC*eU_H;//source term correction operator q_corr expressed as sparse vector eq_corr using multiplication of sparse matrix A_BC and sparse vector eU_H.
            for (size_t cid=0; cid<np; cid++)
            {
            	q_corr[cid] = eq_corr(cid);
            }
            //substitute source term q using correction operator q_corr and calculate potential volume corresponding to i source
            for (size_t m = 0; m < np; m++)
              {
                permutationeq(m) = q_corr[order[m]];
              }
            permutationeu = pFORBCGST.solve(permutationeq);
            for (size_t n = 0; n < np; n++)
              {
                eu(n) = permutationeu(inv_order[n]);
              }
            //select related receivers of i source and calculate forward response,
            size_t count = 0;
            std::vector<double> GeometryFactor;//geometry factor for each electrode pair
            for (size_t Qi = 0; Qi < ndata_res; Qi++)
              {
                if (i == geo.sno[Qi])
                  {
                    Q1 = Linint(cellcenterxpos, cellcenterypos, cellcenterzpos,
                        geo.rx1[Qi], geo.ry1[Qi], geo.rz1[Qi]);
                    Q2 = Linint(cellcenterxpos, cellcenterypos, cellcenterzpos,
                        geo.rx2[Qi], geo.ry2[Qi], geo.rz2[Qi]);
                    Q.reserve(np);
                    for (size_t j = 0; j < grid.nz; j++)
                      for (size_t k = 0; k < grid.ny; k++)
                        for (size_t l = 0; l < grid.nx; l++)
                          {
                            Q.push_back(
                                Q1[l + k * grid.nx + j * grid.nx * grid.ny]
                                    - Q2[l + k * grid.nx + j * grid.nx * grid.ny]);
                          }
                    count++;

                    double R1S1 = std::sqrt((geo.rx1[Qi]-geo.PosSx[i])*(geo.rx1[Qi]-geo.PosSx[i])+(geo.ry1[Qi]-geo.PosSy[i])*(geo.ry1[Qi]-geo.PosSy[i]));
                    double R2S1 = std::sqrt((geo.rx2[Qi]-geo.PosSx[i])*(geo.rx2[Qi]-geo.PosSx[i])+(geo.ry2[Qi]-geo.PosSy[i])*(geo.ry2[Qi]-geo.PosSy[i]));
                    double R1S2 = std::sqrt((geo.rx1[Qi]-geo.NegSx[i])*(geo.rx1[Qi]-geo.NegSx[i])+(geo.ry1[Qi]-geo.NegSy[i])*(geo.ry1[Qi]-geo.NegSy[i]));
                    double R2S2 = std::sqrt((geo.rx2[Qi]-geo.NegSx[i])*(geo.rx2[Qi]-geo.NegSx[i])+(geo.ry2[Qi]-geo.NegSy[i])*(geo.ry2[Qi]-geo.NegSy[i]));
                    GeometryFactor.push_back(2*3.1415926/(1/R1S1-1/R2S1+1/R2S2-1/R1S2));

                  }
              }

            Eigen::MatrixXd eQ(count, np);

            for (size_t Qc = 0; Qc < count; Qc++)
              for (size_t h = 0; h < np; h++)
                {
                  eQ(Qc, h) = Q[Qc * np + h];
                }
            eforwarddata = eQ * eu;

            for (size_t ndatatemp = 0; ndatatemp < count; ndatatemp++)
              {
            	if (GeometryFactor[ndatatemp]==singularityvalue||GeometryFactor[ndatatemp]==-singularityvalue)
            	{
            		result(startindices.at(i) + ndatatemp) = 1/grid.avg_cond;
            		//datatemp.push_back(1/grid.avg_cond);
            	}
            	else
            	{
            		//result(startindices.at(i) + ndatatemp) = eforwarddata(ndatatemp);
            		result(startindices.at(i) + ndatatemp) = eforwarddata(ndatatemp)*GeometryFactor[ndatatemp];
            		//datatemp.push_back(eforwarddata(ndatatemp)*GeometryFactor[ndatatemp]);
                	//datatemp.push_back(eforwarddata(ndatatemp));
            	}
              }
            //GeometryFactor.clear();
          }


        //std::copy(datatemp.begin(), datatemp.end(), result.begin());

        return result;

      }

    /****** Basic gradient calculation function used to get gradient of object function */
    jif3D::rvec ResGradient(const GEOMETRY_RES &geo, const GRID_STRUCT_RES &grid,
        const jif3D::rvec &wdwmisfit)
      {
        std::vector<double> cellwidth_x(grid.nx, grid.dx);
        std::vector<double> cellwidth_y(grid.ny, grid.dy);
        std::vector<double> cellwidth_z(grid.dz);

        //*******Generate index and value of non-zero elements of the divergence matrix D which is a sparse matrix.
        std::vector<size_t> lxd, lyd, lzd, jxd, jyd, jzd; //the index of non-zero elements in sparse matrix D.
        std::vector<double> kxd, kyd, kzd; //the value of non-zero elements in sparse matrix D.
        /*****Divergence matrix D[np][nax+nay+naz], np=grid.nx*grid.ny*grid.nz, nax=(grid.nx-1)*grid.ny*grid.nz,
         * nay=grid.nx*(grid.ny-1)*grid.nz, naz=grid.nx*grid.ny*(grid.nz-1).
         *for (size_t i=0; i<2*nax; i++), D[lxd[i]][jxd[i]]=kxd[i]
         *for (size_t i=0; i<2*nay; i++), D[lyd[i]][jyd[i]+nax]=kyd[i]
         *for (size_t i=0; i<2*naz; i++), D[lzd[i]][jzd[i]+nax+nay]=kzd[i]
         *****/

        //***generate d/dx
        //*Entries(l,j,k)
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 1; k < grid.nx - 1; k++)
              {
                lxd.push_back(k + j * grid.nx + i * grid.nx * grid.ny);
              }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 0; k < grid.nx - 2; k++)
              {
                jxd.push_back(k + j * (grid.nx - 1) + i * (grid.nx - 1) * grid.ny);
              }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 1; k < grid.nx - 1; k++)
              {
                kxd.push_back(-1.0 / (cellwidth_x[k] * 1.0 * 1.0));
              }
        //*Entries(l+1,j,k)
        for (size_t i = 0; i < (grid.nx - 2) * grid.ny * grid.nz; i++)
          {
            lxd.push_back(lxd[i]);
          }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 1; k < grid.nx - 1; k++)
              {
                jxd.push_back(k + j * (grid.nx - 1) + i * (grid.nx - 1) * grid.ny);
              }
        for (size_t i = 0; i < (grid.nx - 2) * grid.ny * grid.nz; i++)
          {
            kxd.push_back(-kxd[i]);
          }
        //*BC at x=0
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            {
              lxd.push_back(j * grid.nx + i * grid.nx * grid.ny);
            }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            {
              jxd.push_back(j * (grid.nx - 1) + i * (grid.nx - 1) * grid.ny);
            }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            {
              kxd.push_back(1.0 / (cellwidth_x[0] * 1.0 * 1.0));
            }
        //*BC at x=end
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            {
              lxd.push_back(grid.nx - 1 + j * grid.nx + i * grid.nx * grid.ny);
            }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            {
              jxd.push_back(
                  grid.nx - 2 + j * (grid.nx - 1) + i * (grid.nx - 1) * grid.ny);
            }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            {
              kxd.push_back(-1.0 / (cellwidth_x[grid.nx - 1] * 1.0 * 1.0));
            }
        //***generate d/dy
        //*Entries(l,j,k)
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 1; j < grid.ny - 1; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                lyd.push_back(k + j * grid.nx + i * grid.nx * grid.ny);
              }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny - 2; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                jyd.push_back(k + j * grid.nx + i * grid.nx * (grid.ny - 1));
              }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 1; j < grid.ny - 1; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                kyd.push_back(-1.0 / (1.0 * cellwidth_y[j] * 1.0));
              }
        //*Entries(l+1,j,k)
        for (size_t i = 0; i < grid.nx * (grid.ny - 2) * grid.nz; i++)
          {
            lyd.push_back(lyd[i]);
          }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 1; j < grid.ny - 1; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                jyd.push_back(k + j * grid.nx + i * grid.nx * (grid.ny - 1));
              }
        for (size_t i = 0; i < grid.nx * (grid.ny - 2) * grid.nz; i++)
          {
            kyd.push_back(-kyd[i]);
          }
        //*BC on y=0
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t k = 0; k < grid.nx; k++)
            {
              lyd.push_back(k + i * grid.nx * grid.ny);
            }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t k = 0; k < grid.nx; k++)
            {
              jyd.push_back(k + i * grid.nx * (grid.ny - 1));
            }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t k = 0; k < grid.nx; k++)
            {
              kyd.push_back(1.0 / (1.0 * cellwidth_y[0] * 1.0));
            }
        //*BC on y=end
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t k = 0; k < grid.nx; k++)
            {
              lyd.push_back(k + grid.nx * (grid.ny - 1) + i * grid.nx * grid.ny);
            }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t k = 0; k < grid.nx; k++)
            {
              jyd.push_back(k + grid.nx * (grid.ny - 2) + i * grid.nx * (grid.ny - 1));
            }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t k = 0; k < grid.nx; k++)
            {
              kyd.push_back(-1.0 / (1.0 * cellwidth_y[grid.ny - 1] * 1.0));
            }
        //***generate d/dz
        //*Entries(l,j,k)
        for (size_t i = 1; i < grid.nz - 1; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                lzd.push_back(k + j * grid.nx + i * grid.nx * grid.ny);
              }
        for (size_t i = 0; i < grid.nz - 2; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                jzd.push_back(k + j * grid.nx + i * grid.nx * grid.ny);
              }
        for (size_t i = 1; i < grid.nz - 1; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                kzd.push_back(-1.0 / (1.0 * 1.0 * cellwidth_z[i]));
              }
        //*Entries(l+1,j,k)
        for (size_t i = 0; i < grid.nx * grid.ny * (grid.nz - 2); i++)
          {
            lzd.push_back(lzd[i]);
          }
        for (size_t i = 1; i < grid.nz - 1; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                jzd.push_back(k + j * grid.nx + i * grid.nx * grid.ny);
              }
        for (size_t i = 0; i < grid.nx * grid.ny * (grid.nz - 2); i++)
          {
            kzd.push_back(-kzd[i]);
          }
        //*BC on z=0
        for (size_t j = 0; j < grid.ny; j++)
          for (size_t k = 0; k < grid.nx; k++)
            {
              lzd.push_back(k + j * grid.nx);
            }
        for (size_t j = 0; j < grid.ny; j++)
          for (size_t k = 0; k < grid.nx; k++)
            {
              jzd.push_back(k + j * grid.nx);
            }
        for (size_t j = 0; j < grid.ny; j++)
          for (size_t k = 0; k < grid.nx; k++)
            {
              kzd.push_back(1.0 / (1.0 * 1.0 * cellwidth_z[0]));
            }
        //*BC on z=end
        for (size_t j = 0; j < grid.ny; j++)
          for (size_t k = 0; k < grid.nx; k++)
            {
              lzd.push_back(k + j * grid.nx + grid.nx * grid.ny * (grid.nz - 1));
            }
        for (size_t j = 0; j < grid.ny; j++)
          for (size_t k = 0; k < grid.nx; k++)
            {
              jzd.push_back(k + j * grid.nx + grid.nx * grid.ny * (grid.nz - 2));
            }
        for (size_t j = 0; j < grid.ny; j++)
          for (size_t k = 0; k < grid.nx; k++)
            {
              kzd.push_back(-1.0 / (1.0 * 1.0 * cellwidth_z[grid.nz - 1]));
            }

        //*******Generate index and value of non-zero elements of the gradient matrix G which is a sparse matrix.
        std::vector<size_t> lxg, lyg, lzg, jxg, jyg, jzg; //the index of non-zero elements in sparse matrix G.
        std::vector<double> kxg, kyg, kzg; //the value of non-zero elements in sparse matrix G.

        /*****Gradient matrix G[nax+nay+naz][np], nax=(grid.nx-1)*grid.ny*grid.nz,
         * nay=grid.nx*(grid.ny-1)*grid.nz, naz=grid.nx*grid.ny*(grid.nz-1),np=grid.nx*grid.ny*grid.nz.
         *for (size_t i=0; i<2*nax; i++), G[lxg[i]][jxg[i]]=kxg[i];
         *for (size_t i=0; i<2*nay; i++), G[lyg[i]+nax][jyg[i]]=kyg[i]
         *for (size_t i=0; i<2*naz; i++), G[lzg[i]+nax+nay][jzg[i]]=kzg[i]
         *****/

        //***generate d/dx
        //*Entries (l,j,k)
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 0; k < grid.nx - 1; k++)
              {
                lxg.push_back(k + j * (grid.nx - 1) + i * (grid.nx - 1) * grid.ny);
              }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 0; k < grid.nx - 1; k++)
              {
                jxg.push_back(k + j * grid.nx + i * grid.nx * grid.ny);
              }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 0; k < grid.nx - 1; k++)
              {
                kxg.push_back(
                    -2.0 / (cellwidth_x[k] * 1.0 * 1.0 + cellwidth_x[k + 1] * 1.0 * 1.0));
              }
        //*Entries (l+1,j,k)
        for (size_t i = 0; i < (grid.nx - 1) * grid.ny * grid.nz; i++)
          {
            lxg.push_back(lxg[i]);
          }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 1; k < grid.nx; k++)
              {
                jxg.push_back(k + j * grid.nx + i * grid.nx * grid.ny);
              }
        for (size_t i = 0; i < (grid.nx - 1) * grid.ny * grid.nz; i++)
          {
            kxg.push_back(-kxg[i]);
          }
        //***generate d/dy
        //*Entries (l,j,k)
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny - 1; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                lyg.push_back(k + j * grid.nx + i * grid.nx * (grid.ny - 1));
              }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny - 1; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                jyg.push_back(k + j * grid.nx + i * grid.nx * grid.ny);
              }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny - 1; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                kyg.push_back(
                    -2.0 / (1.0 * cellwidth_y[j] * 1.0 + 1.0 * cellwidth_y[j + 1] * 1.0));
              }
        //*Entries (l+1,j,k)
        for (size_t i = 0; i < grid.nx * (grid.ny - 1) * grid.nz; i++)
          {
            lyg.push_back(lyg[i]);
          }
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 1; j < grid.ny; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                jyg.push_back(k + j * grid.nx + i * grid.nx * grid.ny);
              }
        for (size_t i = 0; i < grid.nx * (grid.ny - 1) * grid.nz; i++)
          {
            kyg.push_back(-kyg[i]);
          }
        //***generate d/dz
        //*Entries (l,j,k)
        for (size_t i = 0; i < grid.nz - 1; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                lzg.push_back(k + j * grid.nx + i * grid.nx * grid.ny);
              }
        for (size_t i = 0; i < grid.nz - 1; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                jzg.push_back(k + j * grid.nx + i * grid.nx * grid.ny);
              }
        for (size_t i = 0; i < grid.nz - 1; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                kzg.push_back(
                    -2.0 / (1.0 * 1.0 * cellwidth_z[i] + 1.0 * 1.0 * cellwidth_z[i + 1]));
              }
        //*Entries (l+1,j,k)
        for (size_t i = 0; i < grid.nx * grid.ny * (grid.nz - 1); i++)
          {
            lzg.push_back(lzg[i]);
          }
        for (size_t i = 1; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                jzg.push_back(k + j * grid.nx + i * grid.nx * grid.ny);
              }
        for (size_t i = 0; i < grid.nx * grid.ny * (grid.nz - 1); i++)
          {
            kzg.push_back(-kzg[i]);
          }

        //*******Generate index and value of non-zero elements of the diagonal matrix S(rho) which is a sparse matrix.
        std::vector<size_t> lxs, lys, lzs, jxs, jys, jzs; //the index of non-zero elements in sparse matrix S.
        std::vector<double> kxs, kys, kzs; //the value of non-zero elements in sparse matrix S.

        /*****Diagonal matrix S[nax+nay+naz][nax+nay+naz], nax=(grid.nx-1)*grid.ny*grid.nz,
         * nay=grid.nx*(grid.ny-1)*grid.nz, naz=grid.nx*grid.ny*(grid.nz-1).
         *for (size_t i=0; i<nax; i++), S[lxs[i]][jxs[i]]=kxs[i];
         *for (size_t i=0; i<nay; i++), S[lys[i]+nax][jys[i]+nax]=kys[i]
         *for (size_t i=0; i<naz; i++), S[lzs[i]+nax+nay][jzs[i]+nax+nay]=kzs[i]
         *****/

        std::vector<double> rhof_x, rhof_y, rhof_z;
        //***generate x coefficients
        //*Average rho on x face
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 1; k < grid.nx; k++)
              {
                rhof_x.push_back(
                    (cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[i]
                        * grid.rho[k + j * grid.nx + i * grid.nx * grid.ny]
                        + cellwidth_x[k - 1] * cellwidth_y[j] * cellwidth_z[i]
                            * grid.rho[k - 1 + j * grid.nx + i * grid.nx * grid.ny])
                        / (cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[i]
                            + cellwidth_x[k - 1] * cellwidth_y[j] * cellwidth_z[i]));
              }
        for (size_t i = 0; i < (grid.nx - 1) * grid.ny * grid.nz; i++)
          {
            lxs.push_back(i);
            jxs.push_back(i);
            kxs.push_back(rhof_x[i]);
          }
        //***generate y coefficients
        //*Average rho on y face
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 1; j < grid.ny; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                rhof_y.push_back(
                    (cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[i]
                        * grid.rho[k + j * grid.nx + i * grid.nx * grid.ny]
                        + cellwidth_x[k] * cellwidth_y[j - 1] * cellwidth_z[i]
                            * grid.rho[k + (j - 1) * grid.nx + i * grid.nx * grid.ny])
                        / (cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[i]
                            + cellwidth_x[k] * cellwidth_y[j - 1] * cellwidth_z[i]));
              }
        for (size_t i = 0; i < grid.nx * (grid.ny - 1) * grid.nz; i++)
          {
            lys.push_back(i);
            jys.push_back(i);
            kys.push_back(rhof_y[i]);
          }
        //***generate z coefficients
        //*Average rho on z face
        for (size_t i = 1; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                rhof_z.push_back(
                    (cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[i]
                        * grid.rho[k + j * grid.nx + i * grid.nx * grid.ny]
                        + cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[i - 1]
                            * grid.rho[k + j * grid.nx + (i - 1) * grid.nx * grid.ny])
                        / (cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[i]
                            + cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[i - 1]));
              }
        for (size_t i = 0; i < grid.nx * grid.ny * (grid.nz - 1); i++)
          {
            lzs.push_back(i);
            jzs.push_back(i);
            kzs.push_back(rhof_z[i]);
          }

        //*******Generate Forward operator A=D*S*G.
        size_t np = grid.nx * grid.ny * grid.nz;
        size_t nax = (grid.nx - 1) * grid.ny * grid.nz;
        size_t nay = grid.nx * (grid.ny - 1) * grid.nz;
        size_t naz = grid.nx * grid.ny * (grid.nz - 1);

        //generate sparse matrix D,
        typedef Eigen::Triplet<double> DT;
        std::vector<DT> DV;
        DV.reserve(2 * (nax + nay + naz));
        for (size_t i = 0; i < 2 * nax; i++)
          {
            DV.push_back(DT(lxd[i], jxd[i], kxd[i]));
          }
        for (size_t i = 0; i < 2 * nay; i++)
          {
            DV.push_back(DT(lyd[i], jyd[i] + nax, kyd[i]));
          }
        for (size_t i = 0; i < 2 * naz; i++)
          {
            DV.push_back(DT(lzd[i], jzd[i] + nax + nay, kzd[i]));
          }
        Eigen::SparseMatrix<double> D(np, nax + nay + naz);
        D.setFromTriplets(DV.begin(), DV.end());

        //generate sparse matrix S,
        typedef Eigen::Triplet<double> ST;
        std::vector<ST> SV;
        SV.reserve(nax + nay + naz);
        for (size_t i = 0; i < nax; i++)
          {
            SV.push_back(ST(lxs[i], jxs[i], 1.0 / kxs[i]));
          }
        for (size_t i = 0; i < nay; i++)
          {
            SV.push_back(ST(lys[i] + nax, jys[i] + nax, 1.0 / kys[i]));
          }
        for (size_t i = 0; i < naz; i++)
          {
            SV.push_back(ST(lzs[i] + nax + nay, jzs[i] + nax + nay, 1.0 / kzs[i]));
          }
        Eigen::SparseMatrix<double> S(nax + nay + naz, nax + nay + naz);
        S.setFromTriplets(SV.begin(), SV.end());

        //generate sparse matrix G,
        typedef Eigen::Triplet<double> GT;
        std::vector<GT> GV;
        GV.reserve(2 * (nax + nay + naz));
        for (size_t i = 0; i < 2 * nax; i++)
          {
            GV.push_back(GT(lxg[i], jxg[i], kxg[i]));
          }
        for (size_t i = 0; i < 2 * nay; i++)
          {
            GV.push_back(GT(lyg[i] + nax, jyg[i], kyg[i]));
          }
        for (size_t i = 0; i < 2 * naz; i++)
          {
            GV.push_back(GT(lzg[i] + nax + nay, jzg[i], kzg[i]));
          }
        Eigen::SparseMatrix<double> G(nax + nay + naz, np);
        G.setFromTriplets(GV.begin(), GV.end());

        //generate sparse matrix, forward operator A
        Eigen::SparseMatrix<double> A1(np, nax + nay + naz), A(np, np);
        A1 = D * S;
        A = A1 * G;
        A.coeffRef(0, 0) = 1.0 / (cellwidth_x[0] * cellwidth_y[0] * cellwidth_z[0]);

        //*****Permutation of forward operator A using Cuthill-Mckee Ordering Method*********************************//
        /* Because Eigen library has not Cuthill-Mckee Ordering function, but Boost library has, so we used Boost uBlas library to finish Cuthill-Mckee Ordering of forward operator A.
         * When implementing Cuthill-Mckee Ordering using Boost, the forward operator A which is Eigen SparseMatrix need to be input into a Boost SparseMatrix firstly. We use a loop to
         * input the non-zero element of A one by one, so this probably can cause program running slower. But it is difficult to generate forward operator A using Boost C++ library due
         * to no proper multiplication function of sparse matrix can be used. So we use Eigen to generate forward operator A from D*S*G as above. In the future, if there is a new function
         *  which can be used to implement sparse matrix product in Boost uBlas library, it is a probably better method to directly use Boost to generate forward operator A from D*S*G.
         */
        typedef boost::numeric::ublas::compressed_matrix<double> MatType;
        typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
            boost::property<boost::vertex_color_t, boost::default_color_type,
                boost::property<boost::vertex_degree_t, size_t> > > GraphType;
        typedef boost::graph_traits<GraphType>::vertex_descriptor VertexType;
        typedef boost::graph_traits<GraphType>::vertices_size_type size_type;
        GraphType GFOR;
        //declare a Boost sparse matrix FOR, which is equal to forward operator A which is Eigen Sparse Matrix.
        MatType FOR(np, np);
        //order is corresponding to 'p' in Matlab and is used to rearrange sequence of FOR for permutation, inv_order is the same with 'up' in Matlab and is used to regain original sequence.
        boost::numeric::ublas::vector<size_t> order(np), inv_order(order.size());
        //declare row, col number and value of non-zero element of permutation of Sparse Matrix FOR.
        std::vector<size_t> rowpFOR, colpFOR;
        std::vector<double> valuepFOR;
        //Input non-zero element of Eigen Sparse Matrix A into Boost sparse matrix FOR.
        for (int k = 0; k < A.outerSize(); ++k)
          {
            for (SparseMatrix<double>::InnerIterator it(A, k); it; ++it)
              {
                FOR(it.row(), it.col()) = it.value();
              }
          }
        //Input FOR into GFOR, GFOR is a Graph Class. In Boost uBlas C++ library, it is using Graph theory to implement Cuthill-Mckee Ordering of Sparse Matrix.
        for (MatType::const_iterator1 row_it = FOR.begin1(); row_it != FOR.end1();
            ++row_it)
          {
            for (MatType::const_iterator2 col_it = row_it.begin(); col_it != row_it.end();
                ++col_it)
              {
                boost::add_edge(col_it.index1(), col_it.index2(), GFOR);
              }
          }
        // Get order
        boost::property_map<GraphType, boost::vertex_index_t>::type index_map = get(
            boost::vertex_index, GFOR);
        // Use the boost::graph cuthill mckee algorithm
        std::vector<VertexType> inv_perm(boost::num_vertices(GFOR));
        cuthill_mckee_ordering(GFOR, inv_perm.rbegin(), get(boost::vertex_color, GFOR),
            make_degree_map(GFOR));
        for (size_t i = 0; i < inv_perm.size(); i++)
          {
            order[i] = index_map[inv_perm[i]];
          }
        // Get inv_order
        for (size_t i = 0; i < order.size(); i++)
          {
            inv_order[order[i]] = i;
          }
        //Get row,col and value of non-zero of permutation Sparse Matrix of FOR.
        for (MatType::const_iterator1 row_it = FOR.begin1(); row_it != FOR.end1();
            ++row_it)
          {
            for (MatType::const_iterator2 col_it = row_it.begin(); col_it != row_it.end();
                ++col_it)
              {
                size_t row_ind = inv_order(col_it.index1());
                size_t col_ind = inv_order(col_it.index2());
                rowpFOR.push_back(row_ind);
                colpFOR.push_back(col_ind);
                valuepFOR.push_back(*col_it);
              }
          }
        // Then generate permutation sparse matrix of FOR, pFOR and implement Incomplete LU (ILU) of pFOR.
        //generate permutation sparse matrix of FOR, pFOR
        typedef Eigen::Triplet<double> pFORT;
        std::vector<pFORT> pFORV;
        pFORV.reserve(valuepFOR.size());
        for (size_t i = 0; i < valuepFOR.size(); i++)
          {
            pFORV.push_back(pFORT(rowpFOR[i], colpFOR[i], valuepFOR[i]));
          }
        Eigen::SparseMatrix<double> pFOR(np, np);
        pFOR.setFromTriplets(pFORV.begin(), pFORV.end());
        //ILU of pFOR
        Eigen::BiCGSTAB<SparseMatrix<double>, Eigen::IncompleteLUT<double> > pFORBCGST;
        pFORBCGST.preconditioner().setDroptol(0.00001);
        pFORBCGST.compute(pFOR);


        //********************************Boundary correction and source term improvement*********************************//
        /* In order to improve the BCs and have an accurate solution near a source electrode, we apply a correction to the source
         * term q that both improves the BCs and reduces errors near the source locations. The method is as follows,
         * Firstly, define a homogeneous conductivity field Con_H that is the best initial guess of the average conductivity in the model space.
         * Secondly, calculate the analytical solution for the half-space potential field U_H resulting from the conductivity structure Con_H from a source term q.
         * Note, the analytical solution considers a half space that goes to infinity and therefore has the "correct" boundary conditions for our problem.
         * Then, construct the forward operator A(Con_H).
         * Evaluate q_corr = A(Con_H)*U_H.
         * Once we have the vector q_corr, we substitute it for the source term q and solve the system,
         * U_true = A(Con_true)-1*q_corr.
         * Where U_true is the potential field everywhere in our model space and Con_true is the true conductivity field.
         *
         */
        /*Because U_H is calculated for each source pair, we implement it in the loop of calculation potential for true model.
         * Now, we firstly construct the forward operator A(Con_H).
         */
        //*******Generate index and value of non-zero elements of the diagonal matrix S(1/avg_cond) which is a sparse matrix.
        //***We add "_BC" after the similar variable name with that used during constructing Sparse matrix "S" in the above code.
        std::vector<size_t> lxs_BC, lys_BC, lzs_BC, jxs_BC, jys_BC, jzs_BC; //the index of non-zero elements in sparse matrix S(1/avg_cond).
        std::vector<double> kxs_BC, kys_BC, kzs_BC; //the value of non-zero elements in sparse matrix S(1/avg_cond).

        /*****Diagonal matrix S[nax+nay+naz][nax+nay+naz], nax=(grid.nx-1)*grid.ny*grid.nz,
         * nay=grid.nx*(grid.ny-1)*grid.nz, naz=grid.nx*grid.ny*(grid.nz-1).
         *for (size_t i=0; i<nax; i++), S[lxs[i]][jxs[i]]=kxs[i];
         *for (size_t i=0; i<nay; i++), S[lys[i]+nax][jys[i]+nax]=kys[i]
         *for (size_t i=0; i<naz; i++), S[lzs[i]+nax+nay][jzs[i]+nax+nay]=kzs[i]
         *****/

        std::vector<double> rhof_x_BC, rhof_y_BC, rhof_z_BC;
        //***generate x coefficients
        //*Average rho on x face
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 1; k < grid.nx; k++)
              {
                rhof_x_BC.push_back(
                    (cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[i]
                        * (1.0/grid.avg_cond)
                        + cellwidth_x[k - 1] * cellwidth_y[j] * cellwidth_z[i]
                            * (1.0/grid.avg_cond))
                        / (cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[i]
                            + cellwidth_x[k - 1] * cellwidth_y[j] * cellwidth_z[i]));
              }
        for (size_t i = 0; i < (grid.nx - 1) * grid.ny * grid.nz; i++)
          {
            lxs_BC.push_back(i);
            jxs_BC.push_back(i);
            kxs_BC.push_back(rhof_x_BC[i]);
          }
        //***generate y coefficients
        //*Average rho on y face
        for (size_t i = 0; i < grid.nz; i++)
          for (size_t j = 1; j < grid.ny; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                rhof_y_BC.push_back(
                    (cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[i]
                        * (1.0/grid.avg_cond)
                        + cellwidth_x[k] * cellwidth_y[j - 1] * cellwidth_z[i]
                            * (1.0/grid.avg_cond))
                        / (cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[i]
                            + cellwidth_x[k] * cellwidth_y[j - 1] * cellwidth_z[i]));
              }
        for (size_t i = 0; i < grid.nx * (grid.ny - 1) * grid.nz; i++)
          {
            lys_BC.push_back(i);
            jys_BC.push_back(i);
            kys_BC.push_back(rhof_y_BC[i]);
          }
        //***generate z coefficients
        //*Average rho on z face
        for (size_t i = 1; i < grid.nz; i++)
          for (size_t j = 0; j < grid.ny; j++)
            for (size_t k = 0; k < grid.nx; k++)
              {
                rhof_z_BC.push_back(
                    (cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[i]
                        * (1.0/grid.avg_cond)
                        + cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[i - 1]
                            * (1.0/grid.avg_cond))
                        / (cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[i]
                            + cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[i - 1]));
              }
        for (size_t i = 0; i < grid.nx * grid.ny * (grid.nz - 1); i++)
          {
            lzs_BC.push_back(i);
            jzs_BC.push_back(i);
            kzs_BC.push_back(rhof_z_BC[i]);
          }

        //*******construct the forward operator A(Con_H) = D*S(1/avg_cond)*G.
        //generate sparse matrix S(1/avg_cond),
        typedef Eigen::Triplet<double> ST_BC;
        std::vector<ST_BC> SV_BC;
        SV_BC.reserve(nax + nay + naz);
        for (size_t i = 0; i < nax; i++)
          {
            SV_BC.push_back(ST_BC(lxs_BC[i], jxs_BC[i], 1.0 / kxs_BC[i]));
          }
        for (size_t i = 0; i < nay; i++)
          {
            SV_BC.push_back(ST_BC(lys_BC[i] + nax, jys_BC[i] + nax, 1.0 / kys_BC[i]));
          }
        for (size_t i = 0; i < naz; i++)
          {
            SV_BC.push_back(ST_BC(lzs_BC[i] + nax + nay, jzs_BC[i] + nax + nay, 1.0 / kzs_BC[i]));
          }
        Eigen::SparseMatrix<double> S_BC(nax + nay + naz, nax + nay + naz);
        S_BC.setFromTriplets(SV_BC.begin(), SV_BC.end());
        //generate sparse matrix, forward operator A(Con_H) = D*S(1/avg_cond)*G.
        Eigen::SparseMatrix<double> A1_BC(np, nax + nay + naz), A_BC(np, np);
        A1_BC = D * S_BC;
        A_BC = A1_BC * G;
        A_BC.coeffRef(0, 0) = 1.0 / (cellwidth_x[0] * cellwidth_y[0] * cellwidth_z[0]);





        //Convert model vector grid.rho into a diagonal sparse matrix Dm
        typedef Eigen::Triplet<double> DmT;
        std::vector<DmT> DmV;
        DmV.reserve(grid.rho.size());
        for (size_t i = 0; i < grid.rho.size(); i++)
          {
            DmV.push_back(DmT(i, i, grid.rho[i]));
          }
        Eigen::SparseMatrix<double> Dm(grid.rho.size(), grid.rho.size());
        Dm.setFromTriplets(DmV.begin(), DmV.end());





        //*****Calculate gradient of objective function*********************************//
        /* Firstly, calculate centercell's coordinate used in interpolation operation to interpolate source/receiver's into cell's centre.
         */
        std::vector<double> zpos(grid.nz + 1, 0), xpos(grid.nx + 1, 0), ypos(grid.ny + 1,
            0);
        std::vector<double> centerzeroxpos(grid.nx + 1, 0), centerzeroypos(grid.ny + 1,
            0);
        std::vector<double> cellcenterzpos(grid.nz, 0), cellcenterxpos(grid.nx, 0),
            cellcenterypos(grid.ny, 0);
        //build the 3d grid - numbered from 0 to maximum extent,
        for (size_t i = 0; i < grid.nz; i++)
          {
            zpos[i + 1] = zpos[i] + cellwidth_z[i];
          }
        for (size_t i = 0; i < grid.nx; i++)
          {
            xpos[i + 1] = xpos[i] + cellwidth_x[i];
          }
        for (size_t i = 0; i < grid.ny; i++)
          {
            ypos[i + 1] = ypos[i] + cellwidth_y[i];
          }
        //center the grid about zero. x and y's zero point is at the centre of 3d volume, z' zero at surface.
        for (size_t i = 0; i < grid.nx + 1; i++)
          {
            centerzeroxpos[i] = xpos[i] - xpos[grid.nx] / 2.0;
          }
        for (size_t i = 0; i < grid.ny + 1; i++)
          {
            centerzeroypos[i] = ypos[i] - ypos[grid.ny] / 2.0;
          }
        //set cellcenter's position,
        for (size_t i = 0; i < grid.nx; i++)
          {
            cellcenterxpos[i] = centerzeroxpos[i] + cellwidth_x[i] / 2.0;
          }
        for (size_t i = 0; i < grid.ny; i++)
          {
            cellcenterypos[i] = centerzeroypos[i] + cellwidth_y[i] / 2.0;
          }
        for (size_t i = 0; i < grid.nz; i++)
          {
            cellcenterzpos[i] = zpos[i] + cellwidth_z[i] / 2.0;
          }
        /* The follows are the main loop for calculating gradient of objective function
         */

        omp_lock_t lck;
        omp_init_lock(&lck);
        Eigen::MatrixXd eGuii(np, geo.nsource);

#pragma omp parallel for default(shared)
        for (size_t i = 0; i < geo.nsource; i++)
          {
            std::vector<double> pvel(np, 0), nvel(np, 0), pveimagl(np, 0), nveimagl(np, 0), U_H(np, 0), q_corr(np, 0), Q1(np, 0), Q2(np, 0), Q;
            Eigen::VectorXd permutationeq(np), permutationeu, eu(np), eU_H(np), eq_corr(np);

            Eigen::VectorXd eGu(nax + nay + naz);
            std::vector<double> Gux(nax), Guy(nay), Guz(naz);
            std::vector<size_t> lGcx, jGcx, lGcy, jGcy, lGcz, jGcz;
            std::vector<double> kGcx, kGcy, kGcz;
            typedef Eigen::Triplet<double> eGcT;
            std::vector<eGcT> eGcV;
            eGcV.reserve(2 * (nax + nay + naz));
            Eigen::SparseMatrix<double> eGc(nax + nay + naz, np);

            std::vector<double> wdwmisfitforproj;
            Eigen::VectorXd Qtd(np);
            Eigen::VectorXd permutationeQtd(np), permutationelm, elm(np);
            //ILU of pFOR transpose
            Eigen::BiCGSTAB<SparseMatrix<double>, Eigen::IncompleteLUT<double> > pFORBCGSTtrans;
            pFORBCGSTtrans.preconditioner().setDroptol(0.00001);
            Eigen::SparseMatrix<double> pFORtrans(np, np);
            pFORtrans = pFOR.transpose();
            pFORBCGSTtrans.compute(pFORtrans);
            Eigen::VectorXd eGui(np);
            //generate source term q,
            //q1 = Linint(cellcenterxpos, cellcenterypos, cellcenterzpos, geo.PosSx[i],
                //geo.PosSy[i], geo.PosSz[i]);
            //q2 = Linint(cellcenterxpos, cellcenterypos, cellcenterzpos, geo.NegSx[i],
                //geo.NegSy[i], geo.NegSz[i]);
            //for (size_t j = 0; j < grid.nz; j++)
              //for (size_t k = 0; k < grid.ny; k++)
                //for (size_t l = 0; l < grid.nx; l++)
                  //{
                    //q.push_back(
                        //(q1[l + k * grid.nx + j * grid.nx * grid.ny]
                            //- q2[l + k * grid.nx + j * grid.nx * grid.ny])
                            /// (cellwidth_x[l] * cellwidth_y[k] * cellwidth_z[j]));
                  //}
            //norm of current electrode and 1st potential electrode and generate potential of homogeneous model, U_H
            for (size_t l=0; l<grid.nz; l++)
            	for (size_t m=0; m<grid.ny; m++)
            		for (size_t n=0; n<grid.nx; n++)
            		{
            			//norm of current electrode and 1st potential electrode
            			pvel[n+m*grid.nx+l*grid.nx*grid.ny] = std::sqrt((cellcenterxpos[n]-geo.PosSx[i])*(cellcenterxpos[n]-geo.PosSx[i])
            					+(cellcenterypos[m]-geo.PosSy[i])*(cellcenterypos[m]-geo.PosSy[i])+(cellcenterzpos[l]-geo.PosSz[i])*(cellcenterzpos[l]-geo.PosSz[i]));
            			nvel[n+m*grid.nx+l*grid.nx*grid.ny] = std::sqrt((cellcenterxpos[n]-geo.NegSx[i])*(cellcenterxpos[n]-geo.NegSx[i])
            					+(cellcenterypos[m]-geo.NegSy[i])*(cellcenterypos[m]-geo.NegSy[i])+(cellcenterzpos[l]-geo.NegSz[i])*(cellcenterzpos[l]-geo.NegSz[i]));
            			//norm of imaginary current electrode and 1st potential electrode
            			pveimagl[n+m*grid.nx+l*grid.nx*grid.ny] = std::sqrt((cellcenterxpos[n]-geo.PosSx[i])*(cellcenterxpos[n]-geo.PosSx[i])
            					+(cellcenterypos[m]-geo.PosSy[i])*(cellcenterypos[m]-geo.PosSy[i])+(cellcenterzpos[l]+geo.PosSz[i])*(cellcenterzpos[l]+geo.PosSz[i]));
            			nveimagl[n+m*grid.nx+l*grid.nx*grid.ny] = std::sqrt((cellcenterxpos[n]-geo.NegSx[i])*(cellcenterxpos[n]-geo.NegSx[i])
            					+(cellcenterypos[m]-geo.NegSy[i])*(cellcenterypos[m]-geo.NegSy[i])+(cellcenterzpos[l]+geo.NegSz[i])*(cellcenterzpos[l]+geo.NegSz[i]));
            			//generate potential of homogeneous model, U_H
            			U_H[n+m*grid.nx+l*grid.nx*grid.ny] = (1.0/(grid.avg_cond*4.0*3.1415926))*(1.0/pvel[n+m*grid.nx+l*grid.nx*grid.ny]-1.0/nvel[n+m*grid.nx+l*grid.nx*grid.ny]
            					                              +1.0/pveimagl[n+m*grid.nx+l*grid.nx*grid.ny]-1.0/nveimagl[n+m*grid.nx+l*grid.nx*grid.ny]);
            		}
            //now check for singularities due to the source being on a node
            double singularityvalue = std::numeric_limits<double>::infinity();
            for (size_t Uid=0; Uid<np; Uid++)
            {
            	if (U_H[Uid]==singularityvalue||U_H[Uid]==-singularityvalue)
            	{
            		size_t idx, idy, idz;
            		idx = Uid % grid.nx;
            		idz = (Uid - idx) / grid.nx;
            		idy = idz % grid.ny;
            		idz = (idz - idy) / grid.ny;
            		if (idz==0)
            			U_H[Uid] = (U_H[idx+1+idy*grid.nx+idz*grid.nx*grid.ny]+U_H[idx+(idy+1)*grid.nx+idz*grid.nx*grid.ny]
            			               +U_H[idx+idy*grid.nx+(idz+1)*grid.nx*grid.ny]+U_H[idx-1+idy*grid.nx+idz*grid.nx*grid.ny]
            			                   +U_H[idx+(idy-1)*grid.nx+idz*grid.nx*grid.ny])/5.0;
            		else
            			U_H[Uid] = (U_H[idx+1+idy*grid.nx+idz*grid.nx*grid.ny]+U_H[idx+(idy+1)*grid.nx+idz*grid.nx*grid.ny]
            			               +U_H[idx+idy*grid.nx+(idz+1)*grid.nx*grid.ny]+U_H[idx-1+idy*grid.nx+idz*grid.nx*grid.ny]
            			                   +U_H[idx+(idy-1)*grid.nx+idz*grid.nx*grid.ny]+U_H[idx+idy*grid.nx+(idz-1)*grid.nx*grid.ny])/6.0;
            	}
            }
            //calculate correction operator q_corr
            for (size_t cid=0; cid<np; cid++)
            {
            	eU_H(cid) = U_H[cid];
            }
            eq_corr=A_BC*eU_H;//source term correction operator q_corr expressed as sparse vector eq_corr using multiplication of sparse matrix A_BC and sparse vector eU_H.
            for (size_t cid=0; cid<np; cid++)
            {
            	q_corr[cid] = eq_corr(cid);
            }
            //substitute source term q using correction operator q_corr and calculate potential volume corresponding to i source
            for (size_t m = 0; m < np; m++)
              {
                permutationeq(m) = q_corr[order[m]];
              }
            permutationeu = pFORBCGST.solve(permutationeq);
            for (size_t n = 0; n < np; n++)
              {
                eu(n) = permutationeu(inv_order[n]);
              }

            //Creat a sparse matrix eGc containing the derivative of the face matrix
            //Generate x coefficients, x part of non-zero terms for eGc
            eGu = G * eu;
            for (size_t l = 0; l < nax; l++)
              {
                Gux[l] = eGu(l);
              }
            for (size_t l = 0; l < nay; l++)
              {
                Guy[l] = eGu(l + nax);
              }
            for (size_t l = 0; l < naz; l++)
              {
                Guz[l] = eGu(l + nax + nay);
              }
            for (size_t l = 0; l < grid.nz; l++)
              for (size_t j = 0; j < grid.ny; j++)
                for (size_t k = 0; k < grid.nx - 1; k++)
                  {
                    lGcx.push_back(k + j * (grid.nx - 1) + l * (grid.nx - 1) * grid.ny);
                  }
            for (size_t l = 0; l < grid.nz; l++)
              for (size_t j = 0; j < grid.ny; j++)
                for (size_t k = 0; k < grid.nx - 1; k++)
                  {
                    lGcx.push_back(k + j * (grid.nx - 1) + l * (grid.nx - 1) * grid.ny);
                  }
            for (size_t l = 0; l < grid.nz; l++)
              for (size_t j = 0; j < grid.ny; j++)
                for (size_t k = 0; k < grid.nx - 1; k++)
                  {
                    jGcx.push_back(k + j * grid.nx + l * grid.nx * grid.ny);
                  }
            for (size_t l = 0; l < grid.nz; l++)
              for (size_t j = 0; j < grid.ny; j++)
                for (size_t k = 1; k < grid.nx; k++)
                  {
                    jGcx.push_back(k + j * grid.nx + l * grid.nx * grid.ny);
                  }
            for (size_t l = 0; l < grid.nz; l++)
              for (size_t j = 0; j < grid.ny; j++)
                for (size_t k = 1; k < grid.nx; k++)
                  {
                    kGcx.push_back(
                        ((cellwidth_x[k - 1] * cellwidth_y[j] * cellwidth_z[l])
                            / (cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[l]
                                + cellwidth_x[k - 1] * cellwidth_y[j] * cellwidth_z[l]))
                            * Gux[k - 1 + j * (grid.nx - 1) + l * (grid.nx - 1) * grid.ny]);
                  }
            for (size_t l = 0; l < grid.nz; l++)
              for (size_t j = 0; j < grid.ny; j++)
                for (size_t k = 1; k < grid.nx; k++)
                  {
                    kGcx.push_back(
                        ((cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[l])
                            / (cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[l]
                                + cellwidth_x[k - 1] * cellwidth_y[j] * cellwidth_z[l]))
                            * Gux[k - 1 + j * (grid.nx - 1) + l * (grid.nx - 1) * grid.ny]);
                  }
            //Generate y coefficients, y part of non-zero terms for eGc
            for (size_t l = 0; l < grid.nz; l++)
              for (size_t j = 0; j < grid.ny - 1; j++)
                for (size_t k = 0; k < grid.nx; k++)
                  {
                    lGcy.push_back(k + j * grid.nx + l * grid.nx * (grid.ny - 1));
                  }
            for (size_t l = 0; l < grid.nz; l++)
              for (size_t j = 0; j < grid.ny - 1; j++)
                for (size_t k = 0; k < grid.nx; k++)
                  {
                    lGcy.push_back(k + j * grid.nx + l * grid.nx * (grid.ny - 1));
                  }
            for (size_t l = 0; l < grid.nz; l++)
              for (size_t j = 0; j < grid.ny - 1; j++)
                for (size_t k = 0; k < grid.nx; k++)
                  {
                    jGcy.push_back(k + j * grid.nx + l * grid.nx * grid.ny);
                  }
            for (size_t l = 0; l < grid.nz; l++)
              for (size_t j = 1; j < grid.ny; j++)
                for (size_t k = 0; k < grid.nx; k++)
                  {
                    jGcy.push_back(k + j * grid.nx + l * grid.nx * grid.ny);
                  }
            for (size_t l = 0; l < grid.nz; l++)
              for (size_t j = 1; j < grid.ny; j++)
                for (size_t k = 0; k < grid.nx; k++)
                  {
                    kGcy.push_back(
                        ((cellwidth_x[k] * cellwidth_y[j - 1] * cellwidth_z[l])
                            / (cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[l]
                                + cellwidth_x[k] * cellwidth_y[j - 1] * cellwidth_z[l]))
                            * Guy[k + (j - 1) * grid.nx + l * grid.nx * (grid.ny - 1)]);
                  }
            for (size_t l = 0; l < grid.nz; l++)
              for (size_t j = 1; j < grid.ny; j++)
                for (size_t k = 0; k < grid.nx; k++)
                  {
                    kGcy.push_back(
                        ((cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[l])
                            / (cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[l]
                                + cellwidth_x[k] * cellwidth_y[j - 1] * cellwidth_z[l]))
                            * Guy[k + (j - 1) * grid.nx + l * grid.nx * (grid.ny - 1)]);
                  }
            //Generate z coefficients, z part of non-zero terms for eGc
            for (size_t l = 0; l < grid.nz - 1; l++)
              for (size_t j = 0; j < grid.ny; j++)
                for (size_t k = 0; k < grid.nx; k++)
                  {
                    lGcz.push_back(k + j * grid.nx + l * grid.nx * grid.ny);
                  }
            for (size_t l = 0; l < grid.nz - 1; l++)
              for (size_t j = 0; j < grid.ny; j++)
                for (size_t k = 0; k < grid.nx; k++)
                  {
                    lGcz.push_back(k + j * grid.nx + l * grid.nx * grid.ny);
                  }
            for (size_t l = 0; l < grid.nz - 1; l++)
              for (size_t j = 0; j < grid.ny; j++)
                for (size_t k = 0; k < grid.nx; k++)
                  {
                    jGcz.push_back(k + j * grid.nx + l * grid.nx * grid.ny);
                  }
            for (size_t l = 1; l < grid.nz; l++)
              for (size_t j = 0; j < grid.ny; j++)
                for (size_t k = 0; k < grid.nx; k++)
                  {
                    jGcz.push_back(k + j * grid.nx + l * grid.nx * grid.ny);
                  }
            for (size_t l = 1; l < grid.nz; l++)
              for (size_t j = 0; j < grid.ny; j++)
                for (size_t k = 0; k < grid.nx; k++)
                  {
                    kGcz.push_back(
                        ((cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[l - 1])
                            / (cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[l]
                                + cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[l - 1]))
                            * Guz[k + j * grid.nx + (l - 1) * grid.nx * grid.ny]);
                  }
            for (size_t l = 1; l < grid.nz; l++)
              for (size_t j = 0; j < grid.ny; j++)
                for (size_t k = 0; k < grid.nx; k++)
                  {
                    kGcz.push_back(
                        ((cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[l])
                            / (cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[l]
                                + cellwidth_x[k] * cellwidth_y[j] * cellwidth_z[l - 1]))
                            * Guz[k + j * grid.nx + (l - 1) * grid.nx * grid.ny]);
                  }
            //generate sparse matrix eGc,
            for (size_t l = 0; l < 2 * nax; l++)
              {
                eGcV.push_back(eGcT(lGcx[l], jGcx[l], kGcx[l]));
              }
            for (size_t l = 0; l < 2 * nay; l++)
              {
                eGcV.push_back(eGcT(lGcy[l] + nax, jGcy[l], kGcy[l]));
              }
            for (size_t l = 0; l < 2 * naz; l++)
              {
                eGcV.push_back(eGcT(lGcz[l] + nax + nay, jGcz[l], kGcz[l]));
              }
            eGc.setFromTriplets(eGcV.begin(), eGcV.end());

            //Generate a vector lm() related with data misfit and transpose of forward operator A
            //select related receivers of i source and calculate projection matrix from weighted data difference vector,
            size_t count = 0;
            double GeometryFactor;
            for (size_t Qi = 0; Qi < wdwmisfit.size(); Qi++)
              {
                if (i == geo.sno[Qi])
                  {
                    Q1 = Linint(cellcenterxpos, cellcenterypos, cellcenterzpos,
                        geo.rx1[Qi], geo.ry1[Qi], geo.rz1[Qi]);
                    Q2 = Linint(cellcenterxpos, cellcenterypos, cellcenterzpos,
                        geo.rx2[Qi], geo.ry2[Qi], geo.rz2[Qi]);
                    for (size_t j = 0; j < grid.nz; j++)
                      for (size_t k = 0; k < grid.ny; k++)
                        for (size_t l = 0; l < grid.nx; l++)
                          {
                            Q.push_back(
                                Q1[l + k * grid.nx + j * grid.nx * grid.ny]
                                    - Q2[l + k * grid.nx + j * grid.nx * grid.ny]);
                          }

                    double R1S1 = std::sqrt((geo.rx1[Qi]-geo.PosSx[i])*(geo.rx1[Qi]-geo.PosSx[i])+(geo.ry1[Qi]-geo.PosSy[i])*(geo.ry1[Qi]-geo.PosSy[i]));
                    double R2S1 = std::sqrt((geo.rx2[Qi]-geo.PosSx[i])*(geo.rx2[Qi]-geo.PosSx[i])+(geo.ry2[Qi]-geo.PosSy[i])*(geo.ry2[Qi]-geo.PosSy[i]));
                    double R1S2 = std::sqrt((geo.rx1[Qi]-geo.NegSx[i])*(geo.rx1[Qi]-geo.NegSx[i])+(geo.ry1[Qi]-geo.NegSy[i])*(geo.ry1[Qi]-geo.NegSy[i]));
                    double R2S2 = std::sqrt((geo.rx2[Qi]-geo.NegSx[i])*(geo.rx2[Qi]-geo.NegSx[i])+(geo.ry2[Qi]-geo.NegSy[i])*(geo.ry2[Qi]-geo.NegSy[i]));
                    GeometryFactor = 2*3.1415926/(1/R1S1-1/R2S1+1/R2S2-1/R1S2);

                    if (GeometryFactor==singularityvalue||GeometryFactor==-singularityvalue)
                    {
                    	wdwmisfitforproj.push_back(wdwmisfit(Qi));
                    }
                    else
                    {
                        wdwmisfitforproj.push_back(wdwmisfit(Qi)*GeometryFactor);
                        //wdwmisfitforproj.push_back(wdwmisfit(Qi));
                    }

                    count++;

                  }
              }
            Eigen::MatrixXd eQ(count, np);
            Eigen::VectorXd ewdwmisfitforproj(count);
            for (size_t Qc = 0; Qc < count; Qc++)
              for (size_t h = 0; h < np; h++)
                {
                  eQ(Qc, h) = Q[Qc * np + h];
                }
            for (size_t l = 0; l < count; l++)
              {
                ewdwmisfitforproj(l) = wdwmisfitforproj[l];
              }
            Qtd = (eQ.transpose()) * ewdwmisfitforproj;
            //calculate vector lm() related with projection matrix of weighted misfit and the transpose of forward operator A
            for (size_t l = 0; l < np; l++)
              {
                permutationeQtd(l) = -Qtd(order[l]);
              }
            permutationelm = pFORBCGSTtrans.solve(permutationeQtd);
            for (size_t j = 0; j < np; j++)
              {
                elm(j) = permutationelm(inv_order[j]);
              }
            //calculate gradient for each source term
            eGui = (Dm.transpose())
                * ((eGc.transpose()) * (((S * S).transpose()) * ((D.transpose()) * elm)));

            omp_set_lock(&lck);
            for (size_t l = 0; l < np; l++)
              {
                eGuii(l, i) = eGui(l);
              }
            omp_unset_lock(&lck);


          }
        jif3D::rvec gradient(np);
        for (size_t m = 0; m < np; m++)
          {
            double sumeGui = 0.0;
            for (size_t n = 0; n < geo.nsource; n++)
              {
                sumeGui = sumeGui + eGuii(m, n);
              }
            gradient[m] = sumeGui;
          }

        //return gradient;

        jif3D::rvec gradientfinal(np);
        for (size_t i = 0; i < grid.nz; ++i)
          for (size_t j = 0; j < grid.ny; ++j)
            for (size_t k = 0; k < grid.nx; ++k)
              {
                const size_t index1 = i * grid.nx * grid.ny + j * grid.nx + k;
                const size_t index2 = k * grid.nz * grid.ny + j * grid.nz + i;
                gradientfinal[index2] = -2*gradient[index1]/grid.rho[index1];
              }

        return gradientfinal;

      }

  }

