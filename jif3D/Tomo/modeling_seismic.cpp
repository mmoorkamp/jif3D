//============================================================================
// Name        : modeling_seismic.h
// Author      : Jul 6, 2009
// Version     :
// Copyright   : 2009, mmoorkamp (modification) and B. Heincke
//============================================================================

#ifdef HAVEHPX
#include <hpx/hpx.hpp>
#endif

#include "../Tomo/modeling_seismic.h"
#include "../Tomo/ReadWriteTomographyData.h"
#include "../Tomo/PodvinTime3D.h"
#include "../Global/FatalException.h"
#include "../Global/convert.h"
#include <cmath>
#include <vector>
#include <iostream>



jif3D::RayResult ForwardModShotArray(const std::vector<jif3D::GEOMETRY> &geo,
    const jif3D::GRID_STRUCT &grid)
  {
    jif3D::RayResult Result;
    for (size_t i = 0; i < geo.size(); ++i)
      {
        jif3D::RayResult CurrResult = ForwardModShot(geo[i], grid);
        std::copy(CurrResult.raypath.begin(), CurrResult.raypath.end(),
            back_inserter(Result.raypath));
        std::copy(CurrResult.tcalc.begin(), CurrResult.tcalc.end(),
            back_inserter(Result.tcalc));
      }
    return Result;
  }

#ifdef HAVEHPX

HPX_PLAIN_ACTION(ForwardModShotArray, ForwardModShotArray_action)
#endif

namespace jif3D
  {

    /*! Structure to organize the cell parameters during back tracing the rays*/
    struct CELL_STRUCT
      {
      int xno; /*!< position number in x-direction of the cell in the grid*/
      int yno; /*!< position number in y-direction of the cell in the grid*/
      int zno; /*!< position number in z-direction of the cell in the grid*/
      int dirx_i; /*!< The ray runs in negative x-direction into the cell = 1; the ray runs in positive x-direction into the cell = 1; else=0 */
      int diry_i; /*!< The ray runs in negative y-direction into the cell = 1; the ray runs in positive y-direction into the cell = 1; else=0 */
      int dirz_i; /*!< The ray runs in negative z-direction into the cell = 1; the ray runs in positive z-direction into the cell = 1; else=0 */
      double xpos; /*!< The position of the starting point of the ray in x-direction (normalized by the position in the grid)*/
      double ypos; /*!< The position of the starting point of the ray in y-direction (normalized by the position in the grid)*/
      double zpos; /*!< The position of the starting point of the ray in z-direction (normalized by the position in the grid)*/

      };

// Calculate rays from the observed travel time field
    int RayCalc(float *tt, int nx, int ny, int nz, float Xs, float Ys, float Zs,
        float *Xr, float *Yr, float *Zr, int nrec, std::vector<RP_STRUCT> &rp);
    jif3D::rvec TimeGrad(int x, int y, int z, float *tt, int ny, int nz);
    CELL_STRUCT RayBackTrace(double gradx, double grady, double gradz, CELL_STRUCT cell,
        float *tt, int ny, int nz);
    int ResortRays(std::vector<RP_STRUCT> &raypath, const DATA_STRUCT &data,
        const GRID_STRUCT &grid);


    /*-------------------------------------------------------------*/
    /*Performing the forward modeling(using the Podvin&Lecomte eikonal solver) and calculating conventional rays */
    /*Parameter:	 geo  := Geometry structure  */
    /*               grid  := Grid structure*/
    /*              *data  := Pointer on data structure (the calculated traveltimes will be (re-)determined in this routine)*/
    /*				*raypaths := Vector of the raypath structures (number of elements = number of shot-receiver combinations)*/
    /*              time_start := Time when inv3d starts*/

    /*REMARK: The Podvin&Lecomte forward modeling calculates the traveltimes at the grid cell nodes; The number of grid cell nodes is in each*/
    /*		direction one higher than the number of grid cell centers.*/
    /*		Therefore the grid have to be readjusted*/

    int ForwardModRay(const GEOMETRY &geo, const GRID_STRUCT &grid, DATA_STRUCT &data,
        std::vector<RP_STRUCT> &raypath)
      {

        /*******************************************************************************************/
        /*******************************************************************************************/
        /*--Conventional-rays--*/

        std::fill(data.tcalc.begin(), data.tcalc.end(), -1.0);

        if (data.ndata_seis == 0)
          {
            data.ndata_seis_act = 0;

            printf(
                "!!WARNING!! NO seismic traveltime data exists although the seismic forward modelling is activated\n\n");
            printf("Seismic forward modeling is therefore finished\n\n");
            printf("-------------------------------\n\n");
            return (1);
          }

        data.lshots.resize(geo.nshot);
        /*Start the loop over all shots*/
        std::vector<size_t> uniqueshots(data.ndata_seis);
        std::copy(data.sno.begin(), data.sno.end(), uniqueshots.begin());
        std::sort(uniqueshots.begin(), uniqueshots.end());
        uniqueshots.erase(std::unique(uniqueshots.begin(), uniqueshots.end()),
            uniqueshots.end());
        // assert(uniqueshots.size() == geo.nshot);
        //we need a signed version for the openmp loop
        const int nshot = uniqueshots.size();
#ifdef HAVEOPENMP
#pragma omp parallel default(shared)
          {

#pragma omp for
            for (int i = 0; i < nshot; i++)
              {
                size_t count = 0;
                std::vector<size_t> nact_rec; /*active receiver-numbers for the used shot*/
                std::vector<size_t> nact_datapos; /*Position of the active receivers in the data structure*/
                GEOMETRY geo_tmp;
                geo_tmp.x.push_back(geo.x[uniqueshots[i] - 1]);
                geo_tmp.y.push_back(geo.y[uniqueshots[i] - 1]);
                geo_tmp.z.push_back(geo.z[uniqueshots[i] - 1]);
                geo_tmp.nshot = 1;
                for (size_t j = 0; j < data.ndata_seis; j++)
                  {
                    if (uniqueshots[i] == data.sno[j])
                      {
                        nact_rec.push_back(data.rno[j]);
                        nact_datapos.push_back(j);
                        geo_tmp.x.push_back(geo.x[data.rno[j] -1 ]);
                        geo_tmp.y.push_back(geo.y[data.rno[j] -1]);
                        geo_tmp.z.push_back(geo.z[data.rno[j] -1]);
                        count++;
                      }
                  }
                geo_tmp.nrec = count;

                RayResult Rays = ForwardModShot(geo_tmp, grid);

                for (size_t j = 0; j < count; j++)
                  {
                    data.tcalc[nact_datapos[j]] = Rays.tcalc[j];
                    raypath[nact_datapos[j]].nray = Rays.raypath[j].nray;

                    raypath[nact_datapos[j]].len = Rays.raypath[j].len;
                    raypath[nact_datapos[j]].ele = Rays.raypath[j].ele;
                    raypath[nact_datapos[j]].x = Rays.raypath[j].x;
                    raypath[nact_datapos[j]].y = Rays.raypath[j].y;
                    raypath[nact_datapos[j]].z = Rays.raypath[j].z;

                  }
              }

            /*End of the loop over all shots*/
          }
        //end of parallel section
#endif

#ifdef HAVEHPX
        using hpx::future;
        using hpx::async;
        using hpx::wait_all;
        std::vector<hpx::id_type> localities = hpx::find_all_localities();
        const size_t nthreads = hpx::get_num_worker_threads();
        const size_t nlocs = localities.size();
        size_t nchunks = nthreads * nlocs;
        std::vector<hpx::future<RayResult>> ShotResult;
        ShotResult.reserve(nchunks);
        ForwardModShotArray_action ForwardModShotArray;
        const size_t shotsperchunk = nshot / nchunks +1;

        for (size_t c = 0; c < nchunks; ++c)
          {
            size_t count = 0;

            size_t startindex = c * shotsperchunk;
            size_t endindex = std::min(size_t(nshot), (c+1) * shotsperchunk);
            std::vector<GEOMETRY> geo_tmp(endindex - startindex);


            for (size_t i = startindex; i < endindex; i++)
              {
                size_t index = i - startindex;
                geo_tmp[index].nshot = 1;
                geo_tmp[index].x.push_back(geo.x[uniqueshots[i] - 1]);
                geo_tmp[index].y.push_back(geo.y[uniqueshots[i] - 1]);
                geo_tmp[index].z.push_back(geo.z[uniqueshots[i] - 1]);

                for (size_t j = 0; j < data.ndata_seis; j++)
                  {
                    if (uniqueshots[i] == data.sno[j])
                      {
                        geo_tmp[index].x.push_back(geo.x[data.rno[j] -1 ]);
                        geo_tmp[index].y.push_back(geo.y[data.rno[j] -1]);
                        geo_tmp[index].z.push_back(geo.z[data.rno[j] -1]);
                        count++;
                      }
                  }
                geo_tmp[index].nrec = count;
              }

            hpx::id_type const locality_id = localities.at(c % localities.size());
            ShotResult.push_back(async(ForwardModShotArray, locality_id, geo_tmp, grid));
          }
        wait_all(ShotResult);

        for (size_t c = 0; c < nchunks; ++c)
          {

            size_t startindex = c * shotsperchunk;
            size_t endindex = std::min(size_t(nshot), (c+1) * shotsperchunk);
            RayResult Rays = ShotResult[c].get();
            size_t dataindex = 0;
            for (size_t i = startindex; i < endindex; i++)
              {
                size_t count = 0;
                std::vector<size_t> nact_rec; /*active receiver-numbers for the used shot*/
                std::vector<size_t> nact_datapos; /*Position of the active receivers in the data structure*/
                for (size_t j = 0; j < data.ndata_seis; j++)
                  {
                    if (uniqueshots[i] == data.sno[j])
                      {
                        nact_rec.push_back(data.rno[j]);
                        nact_datapos.push_back(j);
                        count++;
                      }
                  }
                for (size_t j = 0; j < count; j++)
                  {
                    data.tcalc[nact_datapos[j]] = Rays.tcalc[dataindex + j];
                    raypath[nact_datapos[j]].nray = Rays.raypath[dataindex + j].nray;

                    raypath[nact_datapos[j]].len = Rays.raypath[dataindex + j].len;
                    raypath[nact_datapos[j]].ele = Rays.raypath[dataindex + j].ele;
                    raypath[nact_datapos[j]].x = Rays.raypath[dataindex + j].x;
                    raypath[nact_datapos[j]].y = Rays.raypath[dataindex + j].y;
                    raypath[nact_datapos[j]].z = Rays.raypath[dataindex + j].z;

                  }
                dataindex += count;
              }
          }
        /*End of the loop over all shots*/

#endif
        /*******************************************************************************************/
        /*******************************************************************************************/
        /*Check the modified structures*/
        for (size_t i = 0; i < data.ndata_seis; i++)
          {
            if (data.tcalc[i] == -1.0)
              {

                throw jif3D::FatalException(
                    "For the shot-receiver combination" + jif3D::stringify(i + 1)
                        + "no traveltime was calculated\n->Check the program\n", __FILE__,
                    __LINE__);
              }

            if (raypath[i].nray % 1 != 0)
              {
                throw jif3D::FatalException(
                    "For the shot-receiver combination" + jif3D::stringify(i + 1)
                        + "no raypath was calculated\n->Check the program\n", __FILE__,
                    __LINE__);
              }
          }

        /****************************************************************************/
        /*Resort rays (If the ray runs along the boundary between to cells, the cell with the higher velocity will be considered)*/
        ResortRays(raypath, data, grid);

        /*Determine the number of ACTIVE rays and the RMS-value*/
        data.ndata_seis_act = 0;
        for (size_t i = 0; i < data.ndata_seis; i++)
          {

            //printf("Number of rays: %d\n", raypath[i].nray);
            if (raypath[i].nray != 0)
              {
                data.ndata_seis_act++;
              }
          }

        //printf("Seismic forward modeling is finished: %d of %d shot-receiver combinations are active\n",
        //    data->ndata_seis_act, data->ndata_seis);
        return (1);
      }

    /*-------------------------------------------------------------*/
    /*Bilinear interpolation in 3-D to determine the velocity from receiver/shot positions to next grid point:  */
    /*Parameter:	x,y,z  := Position of the receiver/shots in "grid cells"  */
    /*              *grid  := Pointer on the grid structure*/
    /*              *data  := Data (usually velocities) from which the interpolated velocities are determined */
    /*Output:  Determined velocity from receiver/shot positions to next grid point */

#define dd(x,y,z)   data[nyz2*(x) + nz2*(y) + (z)]

    inline int lo(const float val)
      {
        return (int) floor((val));
      }

    inline int hi(const float val)
      {
        return (int) ceil((val));
      }
//! Interpolate from the solution at the grid points from the eikonal solver to actual station positions
    float interpolate(float x, float y, float z, const GRID_STRUCT &grid, float *data)
      {
        float u, v, w;
        int ok, nx2, ny2, nz2, nyz2;
        float ival;

        nx2 = grid.nx + 1;
        ny2 = grid.ny + 1;
        nz2 = grid.nz + 1;
        nyz2 = ny2 * nz2;

        /* Check, if point is in grid */
        ok = 1;
        if ((x > nx2) || (x < 0))
          ok = 0;
        if ((y > ny2) || (y < 0))
          ok = 0;
        if ((z > nz2) || (z < 0))
          ok = 0;
        if (!ok)
          {
            std::string error = "Interpolation point is out of the grid! x: "
                + stringify(x) + " y: " + stringify(y) + " z: " + stringify(z) + " nx2: "
                + stringify(nx2) + " ny2: " + stringify(ny2) + " nz2: " + stringify(nz2)
                + "\n";
            throw jif3D::FatalException(error);
          }

        /* Get interpolation distances */
        u = x - (float) floor(x);
        v = y - (float) floor(y);
        w = z - (float) floor(z);

        /* And now interpolate */
        ival = (1 - u) * (1 - v) * (1 - w) * dd(lo(x), lo(y), lo(z))
            + (u) * (1 - v) * (1 - w) * dd(hi(x), lo(y), lo(z))
            + (u) * (v) * (1 - w) * dd(hi(x), hi(y), lo(z))
            + (1 - u) * (v) * (1 - w) * dd(lo(x), hi(y), lo(z)) +

            (1 - u) * (1 - v) * (w) * dd(lo(x), lo(y), hi(z))
            + (u) * (1 - v) * (w) * dd(hi(x), lo(y), hi(z))
            + (u) * (v) * (w) * dd(hi(x), hi(y), hi(z))
            + (1 - u) * (v) * (w) * dd(lo(x), hi(y), hi(z));

        return (ival);
      }

#undef dd

    /*-------------------------------------------------------------*/
    /*Calculation of the ray paths by the means of maximum traveltime descent*/
    /*Parameter:	*tt       := Traveltimes at the grid cell edges  */
    /*              nx,ny,nz  := Size of the grid*/
    /*              Xs,Ys,Zs  := Normalized shot positions*/
    /*              Xr,Yr,Zr  := Normalized receiver positions */
    /*              nrec      := number of receivers per shot */
    /*              *rp       := vector of the ray path structure*/

#define cell_index(x,y,z) ray_cell_index[nyz1*(x) + (nz1)*(y) + (z)]

    int RayCalc(float *tt, int nx, int ny, int nz, float Xs, float Ys, float Zs,
        float *Xr, float *Yr, float *Zr, int nrec, std::vector<RP_STRUCT> &rp)
      {
        int a, b, c;
        int i, count;
        int nx1, ny1, nz1;
        long nyz1;
        std::vector<int> ray_cell_index; /*if 0=no ray in the cell; 1= ray path found in the cell*/

        jif3D::rvec gradient; /*Components of the gradient*/
        CELL_STRUCT next_cell, cell;

        nx1 = nx - 1; /* nx:Number of nodes; nx1= number of cells in x-direction*/
        ny1 = ny - 1;
        nz1 = nz - 1;
        nyz1 = ny1 * nz1;
        const long max_nr_of_ray_seg = 2 * (nx1 + ny1 + nz1); /*max. number of ray-segments*/
        for (i = 0; i < nrec; i++)
          {

            ray_cell_index.resize((nx1) * (ny1) * (nz1));

            /*Set "boundary-index" to 1 for all grid cells:*/
            for (a = 0; a < nx1; a++)
              for (b = 0; b < ny1; b++)
                for (c = 0; c < nz1; c++)
                  {
                    cell_index(a,b,c)= 0;
                  }

                /****************************************************/
                /*Determine the rays starting at the receiver position:*/

                /*Check, if shot and receiver are in the same cell*/
            if (floor((double) Xs) == floor((double) Xr[i])
                && floor((double) Ys) == floor((double) Yr[i])
                && floor((double) Zs) == floor((double) Zr[i]))
              {
                rp[i].len.resize(1);
                rp[i].x.resize(2);
                rp[i].y.resize(2);
                rp[i].z.resize(2);
                rp[i].ele.resize(1);

                rp[i].len[0] = sqrt(
                    (double) (Xs - Xr[i]) * (Xs - Xr[i]) + (Ys - Yr[i]) * (Ys - Yr[i])
                        + (Zs - Zr[i]) * (Zs - Zr[i])); /*Ray segment length*/
                rp[i].x[0] = (double) Xr[i];
                rp[i].y[0] = (double) Yr[i];
                rp[i].z[0] = (double) Zr[i];
                rp[i].x[1] = (double) Xs;
                rp[i].y[1] = (double) Ys;
                rp[i].z[1] = (double) Zs;

                cell.xno = (int) floor((double) Xr[i]);
                cell.yno = (int) floor((double) Yr[i]);
                cell.zno = (int) floor((double) Zr[i]);

                rp[i].ele[0] = nyz1 * cell.xno + nz1 * cell.yno + cell.zno; /*Determine the cell number, where the ray intersects*/
                rp[i].nray = 1; /*Number of segments of the ray*/

                goto fertig;
              }

            /*otherwise ...*/
            /*Set the first cell structure of the considered ray:*/
            cell.xno = (int) floor((double) Xr[i]);
            cell.yno = (int) floor((double) Yr[i]);
            cell.zno = (int) floor((double) Zr[i]);
            cell.xpos = Xr[i] - cell.xno;
            cell.ypos = Yr[i] - cell.yno;
            cell.zpos = Zr[i] - cell.zno;
            cell.dirx_i = 0;
            cell.diry_i = 0;
            cell.dirz_i = 0;

            /*Calculate the traveltime gradient*/
            gradient = TimeGrad(cell.xno, cell.yno, cell.zno, tt, ny, nz);

            /*Calculate the ray segment through the first grid cell*/
            next_cell = RayBackTrace(gradient(0), gradient(1), gradient(2), cell, tt, ny,
                nz);

            rp[i].len.resize(1);
            rp[i].x.resize(1);
            rp[i].y.resize(1);
            rp[i].z.resize(1);
            rp[i].ele.resize(1);

            rp[i].len[0] = sqrt(
                (next_cell.xpos + next_cell.xno - cell.xpos - cell.xno)
                    * (next_cell.xpos + next_cell.xno - cell.xpos - cell.xno)
                    + (next_cell.ypos + next_cell.yno - cell.ypos - cell.yno)
                        * (next_cell.ypos + next_cell.yno - cell.ypos - cell.yno)
                    + (next_cell.zpos + next_cell.zno - cell.zpos - cell.zno)
                        * (next_cell.zpos + next_cell.zno - cell.zpos - cell.zno)); /*Ray segment length*/
            rp[i].x[0] = (double) Xr[i];
            rp[i].y[0] = (double) Yr[i];
            rp[i].z[0] = (double) Zr[i];
            rp[i].ele[0] = nyz1 * cell.xno + nz1 * cell.yno + cell.zno; /*Determine the position number of the cell, which the ray intersects*/

            cell_index(cell.xno,cell.yno,cell.zno)= 1;

            /*Check, if the ray leave the cell*/
            if (next_cell.xno == 0 || next_cell.xno == nx1 || next_cell.yno == 0
                || next_cell.yno == ny1 || next_cell.zno == 0 || next_cell.zno == nz1)
              {
                //printf(
                //    "The ray from the shot-receiver\ncombination %d leaves the model\n\n",
                //    rp[i].n + 1);
                //rp[i].nray = 0;
                goto fertig;
              }

            count = 1;
            /******************************************************/
            /*Calculate iteratively the ray segments through the grid*/

            while (fabs(Xs - next_cell.xpos - next_cell.xno) > 1
                || fabs(Ys - next_cell.ypos - next_cell.yno) > 1
                || fabs(Zs - next_cell.zpos - next_cell.zno) > 1)
              {
                cell = next_cell;

                /*Calculate the traveltime gradient*/
                gradient = TimeGrad(cell.xno, cell.yno, cell.zno, tt, ny, nz);

                /*Calculate the ray segment through the corresponding grid cell*/
                next_cell = RayBackTrace(gradient[0], gradient[1], gradient[2], cell, tt,
                    ny, nz);

                rp[i].len.push_back(
                    sqrt(
                        (next_cell.xpos + next_cell.xno - cell.xpos - cell.xno)
                            * (next_cell.xpos + next_cell.xno - cell.xpos - cell.xno)
                            + (next_cell.ypos + next_cell.yno - cell.ypos - cell.yno)
                                * (next_cell.ypos + next_cell.yno - cell.ypos - cell.yno)
                            + (next_cell.zpos + next_cell.zno - cell.zpos - cell.zno)
                                * (next_cell.zpos + next_cell.zno - cell.zpos - cell.zno))); /*Ray segment length*/
                rp[i].x.push_back((double) cell.xpos + cell.xno);
                rp[i].y.push_back((double) cell.ypos + cell.yno);
                rp[i].z.push_back((double) cell.zpos + cell.zno);
                rp[i].ele.push_back(nyz1 * cell.xno + nz1 * cell.yno + cell.zno); /*Determine the position number of the cell, which the ray intersects*/

                cell_index(cell.xno,cell.yno,cell.zno)= 1;

                /*Check, if the ray leave the cell*/
                if (next_cell.xno == 0 || next_cell.xno == nx1 || next_cell.yno == 0
                    || next_cell.yno == ny1 || next_cell.zno == 0 || next_cell.zno == nz1)
                  {
                    //printf(
                    //    "The ray from the shot-receiver\ncombination %d leaves the model\n\n",
                    //    rp[i].n + 1);
                    //rp[i].nray = 0;
                    goto fertig;
                  }

                /*The ray will not be traced back, if the number of ray segments become too large; the ray is "fallen" probably in a local minima*/
                if (count >= max_nr_of_ray_seg)
                  {
                    //printf(
                    //    "The discretized traveltime field of the ray from the shot-receiver\ncombination %d had a probably local minima\n\n",
                    //   rp[i].n + 1);
                    rp[i].nray = 0;
                    goto fertig;
                  }

                count++;
              }

            rp[i].x.push_back((double) next_cell.xpos + next_cell.xno);
            rp[i].y.push_back((double) next_cell.ypos + next_cell.yno);
            rp[i].z.push_back((double) next_cell.zpos + next_cell.zno);

            rp[i].len.push_back(
                sqrt(
                    (Xs - next_cell.xpos - next_cell.xno)
                        * (Xs - next_cell.xpos - next_cell.xno)
                        + (Ys - next_cell.ypos - next_cell.yno)
                            * (Ys - next_cell.ypos - next_cell.yno)
                        + (Zs - next_cell.zpos - next_cell.zno)
                            * (Zs - next_cell.zpos - next_cell.zno))); /*Ray segment length*/
            rp[i].x.push_back((double) Xs);
            rp[i].y.push_back((double) Ys);
            rp[i].z.push_back((double) Zs);
            rp[i].ele.push_back(
                nyz1 * (int) floor((double) Xs) + nz1 * (int) floor((double) Ys)
                    + (int) floor((double) Zs)); /*Determine the position number of the cell, which the ray intersects*/
            rp[i].nray = count + 1; /*Number of the segments of the ray*/

            fertig: ;

          }

        //printf("All rays for the shot are calculated\n");
        //printf("----------------\n\n");

        return (1);
      }

#undef cell_index
    /*-------------------------------------------------------------*/
    /*Calculation the components of the time gradient*/
    /*Parameter:	*tt			:= Traveltimes at the grid cell edges  */
    /*				x,y,z	:= Position of the cell in the grid*/
    /*				ny,nz	:= Size of the grid (number of grid cell EDGES)*/

    /*	Output: Components of the gradient */

#define Traveltimes(a,b,c) tt[nyz*(a) + nz*(b) + (c)]

    jif3D::rvec TimeGrad(int x, int y, int z, float *tt, int ny, int nz)
      {
        int nyz;
        jif3D::rvec grad(3);

        nyz = ny * nz;

        grad[0] = (-Traveltimes(x + 1, y, z) - Traveltimes(x + 1, y, z + 1)
            - Traveltimes(x + 1, y + 1, z) - Traveltimes(x + 1, y + 1, z + 1)
            + Traveltimes(x, y, z) + Traveltimes(x, y, z + 1) + Traveltimes(x, y + 1, z)
            + Traveltimes(x, y + 1, z + 1)) / 4; /*x-component*/

        grad[1] = (-Traveltimes(x, y + 1, z) - Traveltimes(x, y + 1, z + 1)
            - Traveltimes(x + 1, y + 1, z) - Traveltimes(x + 1, y + 1, z + 1)
            + Traveltimes(x, y, z) + Traveltimes(x, y, z + 1) + Traveltimes(x + 1, y, z)
            + Traveltimes(x + 1, y, z + 1)) / 4; /*y-component*/

        grad[2] = (-Traveltimes(x, y, z + 1) - Traveltimes(x, y + 1, z + 1)
            - Traveltimes(x + 1, y, z + 1) - Traveltimes(x + 1, y + 1, z + 1)
            + Traveltimes(x, y, z) + Traveltimes(x, y + 1, z) + Traveltimes(x + 1, y, z)
            + Traveltimes(x + 1, y + 1, z)) / 4; /*z-component*/

        return (grad);
      }

#undef Traveltimes

    /*-------------------------------------------------------------*/
    /*Determine the ray segments through a cell in 3-D*/
    /*Parameter:	gradx,grady,gradz	:= (negative) Traveltime-Gradient */
    /*				cell				:= Cell structure of the actual cell*/
    /*				*tt					:= traveltimes */
    /*              ny,nz			:= grid size */

    /*	Output: Cell structure of the next cell reached by the ray */

    CELL_STRUCT RayBackTrace(double gradx, double grady, double gradz, CELL_STRUCT cell,
        float *tt, int ny, int nz)
      {
        double eps = 0.01; /*Stabilize the program;*/
        double tmp_xpos, tmp_ypos, tmp_zpos;
        jif3D::rvec gradient1, gradient2;
        int diff;
        CELL_STRUCT next_cell;

        /*Safety settings to avoid instabilties (move the points sligthy away from the grid boundary)*/
        if (cell.dirx_i == 0 && cell.xpos == 0)
          cell.xpos = cell.xpos + eps;
        if (cell.dirx_i == 0 && cell.xpos == 1)
          cell.xpos = cell.xpos - eps;
        if (cell.diry_i == 0 && cell.ypos == 0)
          cell.ypos = cell.ypos + eps;
        if (cell.diry_i == 0 && cell.ypos == 1)
          cell.ypos = cell.ypos - eps;
        if (cell.dirz_i == 0 && cell.zpos == 0)
          cell.zpos = cell.zpos + eps;
        if (cell.dirz_i == 0 && cell.zpos == 1)
          cell.zpos = cell.zpos - eps;

        /*Ray enter the cell from the upper/lower x-plane*/
        /*CASE I: The ray run into the cell*/
        if (cell.dirx_i == 0 || (cell.dirx_i == 1 && gradx <= 0)
            || (cell.dirx_i == -1 && gradx >= 0))
          {
            if (gradx > 0)
              {

                tmp_ypos = cell.ypos - ((cell.xpos - 1) * (grady / gradx));
                tmp_zpos = cell.zpos - ((cell.xpos - 1) * (gradz / gradx));

                if (1 >= tmp_ypos && tmp_ypos >= 0 && 1 >= tmp_zpos && tmp_zpos >= 0)
                  {
                    next_cell.xpos = 0;
                    next_cell.ypos = tmp_ypos;
                    next_cell.zpos = tmp_zpos;
                    next_cell.dirx_i = -1;
                    next_cell.diry_i = 0;
                    next_cell.dirz_i = 0;
                    next_cell.xno = cell.xno + 1;
                    next_cell.yno = cell.yno;
                    next_cell.zno = cell.zno;

                    goto bestimmt;

                  }
              }

            if (gradx < 0)
              {
                tmp_ypos = cell.ypos - (cell.xpos * (grady / gradx));
                tmp_zpos = cell.zpos - (cell.xpos * (gradz / gradx));

                if (1 >= tmp_ypos && tmp_ypos >= 0 && 1 >= tmp_zpos && tmp_zpos >= 0)
                  {
                    next_cell.xpos = 1;
                    next_cell.ypos = tmp_ypos;
                    next_cell.zpos = tmp_zpos;
                    next_cell.dirx_i = 1;
                    next_cell.diry_i = 0;
                    next_cell.dirz_i = 0;
                    next_cell.xno = cell.xno - 1;
                    next_cell.yno = cell.yno;
                    next_cell.zno = cell.zno;

                    goto bestimmt;
                  }

              }
          }

        /*Ray enter the cell from the positive or negative y-plane*/
        /*CASE I: The ray run into the cell*/
        if (cell.diry_i == 0 || (cell.diry_i == 1 && grady <= 0)
            || (cell.diry_i == -1 && grady >= 0))
          {
            if (grady > 0)
              {
                tmp_xpos = cell.xpos - ((cell.ypos - 1) * (gradx / grady));
                tmp_zpos = cell.zpos - ((cell.ypos - 1) * (gradz / grady));

                if (1 >= tmp_xpos && tmp_xpos >= 0 && 1 >= tmp_zpos && tmp_zpos >= 0)
                  {
                    next_cell.xpos = tmp_xpos;
                    next_cell.ypos = 0;
                    next_cell.zpos = tmp_zpos;
                    next_cell.dirx_i = 0;
                    next_cell.diry_i = -1;
                    next_cell.dirz_i = 0;
                    next_cell.xno = cell.xno;
                    next_cell.yno = cell.yno + 1;
                    next_cell.zno = cell.zno;

                    goto bestimmt;
                  }
              }

            if (grady < 0)
              {
                tmp_xpos = cell.xpos - (cell.ypos * (gradx / grady));
                tmp_zpos = cell.zpos - (cell.ypos * (gradz / grady));

                if (1 >= tmp_xpos && tmp_xpos >= 0 && 1 >= tmp_zpos && tmp_zpos >= 0)
                  {
                    next_cell.xpos = tmp_xpos;
                    next_cell.ypos = 1;
                    next_cell.zpos = tmp_zpos;
                    next_cell.dirx_i = 0;
                    next_cell.diry_i = 1;
                    next_cell.dirz_i = 0;
                    next_cell.xno = cell.xno;
                    next_cell.yno = cell.yno - 1;
                    next_cell.zno = cell.zno;

                    goto bestimmt;
                  }

              }

          }

        /*Ray enter the cell from the positive or negative z-plane*/
        /*CASE I:The ray run into the cell*/
        if (cell.dirz_i == 0 || (cell.dirz_i == 1 && gradz <= 0)
            || (cell.dirz_i == -1 && gradz >= 0))
          {
            if (gradz > 0)
              {
                tmp_xpos = cell.xpos - ((cell.zpos - 1) * (gradx / gradz));
                tmp_ypos = cell.ypos - ((cell.zpos - 1) * (grady / gradz));

                if (1 >= tmp_xpos && tmp_xpos >= 0 && 1 >= tmp_ypos && tmp_ypos >= 0)
                  {
                    next_cell.xpos = tmp_xpos;
                    next_cell.ypos = tmp_ypos;
                    next_cell.zpos = 0;
                    next_cell.dirx_i = 0;
                    next_cell.diry_i = 0;
                    next_cell.dirz_i = -1;
                    next_cell.xno = cell.xno;
                    next_cell.yno = cell.yno;
                    next_cell.zno = cell.zno + 1;

                    goto bestimmt;
                  }
              }

            if (gradz < 0)
              {
                tmp_xpos = cell.xpos - (cell.zpos * (gradx / gradz));
                tmp_ypos = cell.ypos - (cell.zpos * (grady / gradz));

                if (1 >= tmp_xpos && tmp_xpos >= 0 && 1 >= tmp_ypos && tmp_ypos >= 0)
                  {
                    next_cell.xpos = tmp_xpos;
                    next_cell.ypos = tmp_ypos;
                    next_cell.zpos = 1;
                    next_cell.dirx_i = 0;
                    next_cell.diry_i = 0;
                    next_cell.dirz_i = 1;
                    next_cell.xno = cell.xno;
                    next_cell.yno = cell.yno;
                    next_cell.zno = cell.zno - 1;

                    goto bestimmt;
                  }

              }
          }

        /*Ray enter the cell from the positive or negative x-plane*/
        /*CASE II: The ray run along the surface of the cell*/
        if (((cell.dirx_i == 1 && gradx > 0) || (cell.dirx_i == -1 && gradx < 0)))
          {
            if (grady < 0)
              {
                tmp_zpos = cell.zpos - cell.ypos * (gradz / grady);

                if (1 >= tmp_zpos && tmp_zpos >= 0)
                  {

                    /*Comparison of the gradients(in the neighboring cells) to determine the next cell:*/
                    /*Calculate the traveltime gradient*/
                    gradient1 = TimeGrad(cell.xno + cell.dirx_i, cell.yno - 1, cell.zno,
                        tt, ny, nz);
                    gradient2 = TimeGrad(cell.xno, cell.yno - 1, cell.zno, tt, ny, nz);
                    if (gradient1[0] + gradient2[0] >= 0)
                      diff = cell.dirx_i;
                    else
                      diff = -cell.dirx_i;

                    if ((gradient1[0] + gradient2[0]) >= 0)
                      next_cell.xpos = 0;
                    else
                      next_cell.xpos = 1;
                    next_cell.ypos = 1;
                    next_cell.zpos = tmp_zpos;
                    next_cell.dirx_i = 0;
                    next_cell.diry_i = 1;
                    next_cell.dirz_i = 0;
                    if (diff > 0)
                      next_cell.xno = cell.xno + cell.dirx_i;
                    else
                      next_cell.xno = cell.xno;
                    next_cell.yno = cell.yno - 1;
                    next_cell.zno = cell.zno;

                    goto bestimmt;
                  }
              }

            if (grady > 0)
              {
                tmp_zpos = cell.zpos - (cell.ypos - 1) * (gradz / grady);

                if (1 >= tmp_zpos && tmp_zpos >= 0)
                  {

                    /*Comparison of the gradients(in the neighboring cells) to determine the next cell:*/
                    /*Calculate the traveltime gradient*/
                    gradient1 = TimeGrad(cell.xno + cell.dirx_i, cell.yno + 1, cell.zno,
                        tt, ny, nz);
                    gradient2 = TimeGrad(cell.xno, cell.yno + 1, cell.zno, tt, ny, nz);
                    if (gradient1[0] + gradient2[0] >= 0)
                      diff = cell.dirx_i;
                    else
                      diff = -cell.dirx_i;

                    if (gradient1[0] + gradient2[0] >= 0)
                      next_cell.xpos = 0;
                    else
                      next_cell.xpos = 1;
                    next_cell.ypos = 0;
                    next_cell.zpos = tmp_zpos;
                    next_cell.dirx_i = 0;
                    next_cell.diry_i = -1;
                    next_cell.dirz_i = 0;
                    if (diff > 0)
                      next_cell.xno = cell.xno + cell.dirx_i;
                    else
                      next_cell.xno = cell.xno;
                    next_cell.yno = cell.yno + 1;
                    next_cell.zno = cell.zno;

                    goto bestimmt;
                  }
              }

            if (gradz < 0)
              {
                tmp_ypos = cell.ypos - cell.zpos * (grady / gradz);

                if (1 >= tmp_ypos && tmp_ypos >= 0)
                  {

                    /*Comparison of the gradients(in the neighboring cells) to determine the next cell:*/
                    /*Calculate the traveltime gradient*/
                    gradient1 = TimeGrad(cell.xno + cell.dirx_i, cell.yno, cell.zno - 1,
                        tt, ny, nz);
                    gradient2 = TimeGrad(cell.xno, cell.yno, cell.zno - 1, tt, ny, nz);
                    if (gradient1[0] + gradient2[0] >= 0)
                      diff = cell.dirx_i;
                    else
                      diff = -cell.dirx_i;

                    if (gradient1[0] + gradient2[0] >= 0)
                      next_cell.xpos = 0;
                    else
                      next_cell.xpos = 1;
                    next_cell.ypos = tmp_ypos;
                    next_cell.zpos = 1;
                    next_cell.dirx_i = 0;
                    next_cell.diry_i = 0;
                    next_cell.dirz_i = 1;

                    if (diff > 0)
                      next_cell.xno = cell.xno + cell.dirx_i;
                    else
                      next_cell.xno = cell.xno;
                    next_cell.yno = cell.yno;
                    next_cell.zno = cell.zno - 1;

                    goto bestimmt;
                  }
              }

            if (gradz > 0)
              {
                tmp_ypos = cell.ypos - (cell.zpos - 1) * (grady / gradz);

                if (1 >= tmp_ypos && tmp_ypos >= 0)
                  {
                    /*Comparison of the gradients(in the neighboring cells) to determine the next cell:*/
                    /*Calculate the traveltime gradient*/
                    gradient1 = TimeGrad(cell.xno + cell.dirx_i, cell.yno, cell.zno + 1,
                        tt, ny, nz);
                    gradient2 = TimeGrad(cell.xno, cell.yno, cell.zno + 1, tt, ny, nz);
                    if (gradient1[0] + gradient2[0] >= 0)
                      diff = cell.dirx_i;
                    else
                      diff = -cell.dirx_i;

                    if (gradient1[0] + gradient2[0] >= 0)
                      next_cell.xpos = 0;
                    else
                      next_cell.xpos = 1;
                    next_cell.ypos = tmp_ypos;
                    next_cell.zpos = 0;
                    next_cell.dirx_i = 0;
                    next_cell.diry_i = 0;
                    next_cell.dirz_i = -1;
                    if (diff > 0)
                      next_cell.xno = cell.xno + cell.dirx_i;
                    else
                      next_cell.xno = cell.xno;
                    next_cell.yno = cell.yno;
                    next_cell.zno = cell.zno + 1;

                    goto bestimmt;
                  }
              }

          }

        /*Ray enter the cell from the positive or negative y-plane*/
        /*CASE II: The ray run along the surface of the cell*/
        if (((cell.diry_i == 1 && grady > 0) || (cell.diry_i == -1 && grady < 0)))
          {
            if (gradx < 0)
              {
                tmp_zpos = cell.zpos - cell.xpos * (gradz / gradx);

                if (1 >= tmp_zpos && tmp_zpos >= 0)
                  {
                    /*Comparison of the gradients(in the neighboring cells) to determine the next cell:*/
                    /*Calculate the traveltime gradient*/
                    gradient1 = TimeGrad(cell.xno - 1, cell.yno + cell.diry_i, cell.zno,
                        tt, ny, nz);
                    gradient2 = TimeGrad(cell.xno - 1, cell.yno, cell.zno, tt, ny, nz);
                    if (gradient1[1] + gradient2[1] >= 0)
                      diff = cell.diry_i;
                    else
                      diff = -cell.diry_i;

                    next_cell.xpos = 1;
                    if (gradient1[1] + gradient2[1] >= 0)
                      next_cell.ypos = 0;
                    else
                      next_cell.ypos = 1;
                    next_cell.zpos = tmp_zpos;
                    next_cell.dirx_i = 1;
                    next_cell.diry_i = 0;
                    next_cell.dirz_i = 0;
                    next_cell.xno = cell.xno - 1;
                    if (diff > 0)
                      next_cell.yno = cell.yno + cell.diry_i;
                    else
                      next_cell.yno = cell.yno;
                    next_cell.zno = cell.zno;
                    goto bestimmt;
                  }
              }

            if (gradx > 0)
              {
                tmp_zpos = cell.zpos - (cell.xpos - 1) * (gradz / gradx);

                if (1 >= tmp_zpos && tmp_zpos >= 0)
                  {
                    /*Comparison of the gradients(in the neighboring cells) to determine the next cell:*/
                    /*Calculate the traveltime gradient*/
                    gradient1 = TimeGrad(cell.xno + 1, cell.yno + cell.diry_i, cell.zno,
                        tt, ny, nz);
                    gradient2 = TimeGrad(cell.xno + 1, cell.yno, cell.zno, tt, ny, nz);
                    if (gradient1[1] + gradient2[1] >= 0)
                      diff = cell.diry_i;
                    else
                      diff = -cell.diry_i;

                    next_cell.xpos = 0;
                    if (gradient1[1] + gradient2[1] >= 0)
                      next_cell.ypos = 0;
                    else
                      next_cell.ypos = 1;
                    next_cell.zpos = tmp_zpos;
                    next_cell.dirx_i = -1;
                    next_cell.diry_i = 0;
                    next_cell.dirz_i = 0;
                    next_cell.xno = cell.xno + 1;
                    if (diff > 0)
                      next_cell.yno = cell.yno + cell.diry_i;
                    else
                      next_cell.yno = cell.yno;
                    next_cell.zno = cell.zno;

                    goto bestimmt;
                  }
              }
            if (gradz < 0)
              {
                tmp_xpos = cell.xpos - cell.zpos * (gradx / gradz);

                if (1 >= tmp_xpos && tmp_xpos >= 0)
                  {
                    /*Comparison of the gradients(in the neighboring cells) to determine the next cell:*/
                    /*Calculate the traveltime gradient*/
                    gradient1 = TimeGrad(cell.xno, cell.yno + cell.diry_i, cell.zno - 1,
                        tt, ny, nz);
                    gradient2 = TimeGrad(cell.xno, cell.yno, cell.zno - 1, tt, ny, nz);
                    if (gradient1[1] + gradient2[1] >= 0)
                      diff = cell.diry_i;
                    else
                      diff = -cell.diry_i;

                    next_cell.xpos = tmp_xpos;
                    if (gradient1[1] + gradient2[1] >= 0)
                      next_cell.ypos = 0;
                    else
                      next_cell.ypos = 1;
                    next_cell.zpos = 1;
                    next_cell.dirx_i = 0;
                    next_cell.diry_i = 0;
                    next_cell.dirz_i = 1;

                    next_cell.xno = cell.xno;
                    if (diff > 0)
                      next_cell.yno = cell.yno + cell.diry_i;
                    else
                      next_cell.yno = cell.yno;
                    next_cell.zno = cell.zno - 1;

                    goto bestimmt;
                  }
              }

            if (gradz > 0)
              {
                tmp_xpos = cell.xpos - (cell.zpos - 1) * (gradx / gradz);

                if (1 >= tmp_xpos && tmp_xpos >= 0)
                  {

                    /*Comparison of the gradients(in the neighboring cells) to determine the next cell:*/
                    /*Calculate the traveltime gradient*/
                    gradient1 = TimeGrad(cell.xno, cell.yno + cell.diry_i, cell.zno + 1,
                        tt, ny, nz);
                    gradient2 = TimeGrad(cell.xno, cell.yno, cell.zno + 1, tt, ny, nz);
                    if (gradient1[1] + gradient2[1] >= 0)
                      diff = cell.diry_i;
                    else
                      diff = -cell.diry_i;

                    next_cell.xpos = tmp_xpos;
                    if (gradient1[1] + gradient2[1] >= 0)
                      next_cell.ypos = 0;
                    else
                      next_cell.ypos = 1;
                    next_cell.zpos = 0;
                    next_cell.dirx_i = 0;
                    next_cell.diry_i = 0;
                    next_cell.dirz_i = -1;
                    next_cell.xno = cell.xno;
                    if (diff > 0)
                      next_cell.yno = cell.yno + cell.diry_i;
                    else
                      next_cell.yno = cell.yno;
                    next_cell.zno = cell.zno + 1;

                    goto bestimmt;
                  }
              }

          }

        /*Ray enter the cell from the positive or negative z-plane*/
        /*CASE II: The ray run along the surface of the cell*/
        if (((cell.dirz_i == 1 && gradz > 0) || (cell.dirz_i == -1 && gradz < 0)))
          {
            if (gradx < 0)
              {
                tmp_ypos = cell.ypos - cell.xpos * (grady / gradx);

                if (1 >= tmp_ypos && tmp_ypos >= 0)
                  {
                    /*Comparison of the gradients(in the neighboring cells) to determine the next cell:*/
                    /*Calculate the traveltime gradient*/
                    gradient1 = TimeGrad(cell.xno - 1, cell.yno, cell.zno + cell.dirz_i,
                        tt, ny, nz);
                    gradient2 = TimeGrad(cell.xno - 1, cell.yno, cell.zno, tt, ny, nz);
                    if (gradient1[2] + gradient2[2] >= 0)
                      diff = cell.dirz_i;
                    else
                      diff = -cell.dirz_i;

                    next_cell.xpos = 1;
                    next_cell.ypos = tmp_ypos;
                    if (gradient1[2] + gradient2[2] >= 0)
                      next_cell.zpos = 0;
                    else
                      next_cell.zpos = 1;
                    next_cell.dirx_i = 1;
                    next_cell.diry_i = 0;
                    next_cell.dirz_i = 0;

                    next_cell.xno = cell.xno - 1;
                    next_cell.yno = cell.yno;
                    if (diff > 0)
                      next_cell.zno = cell.zno + cell.dirz_i;
                    else
                      next_cell.zno = cell.zno;

                    goto bestimmt;
                  }
              }

            if (gradx > 0)
              {
                tmp_ypos = cell.ypos - (cell.xpos - 1) * (grady / gradx);

                if (1 >= tmp_ypos && tmp_ypos >= 0)
                  {

                    /*Comparison of the gradients(in the neighboring cells) to determine the next cell:*/
                    /*Calculate the traveltime gradient*/
                    gradient1 = TimeGrad(cell.xno + 1, cell.yno, cell.zno + cell.dirz_i,
                        tt, ny, nz);
                    gradient2 = TimeGrad(cell.xno + 1, cell.yno, cell.zno, tt, ny, nz);
                    if (gradient1[2] + gradient2[2] >= 0)
                      diff = cell.dirz_i;
                    else
                      diff = -cell.dirz_i;

                    next_cell.xpos = 0;
                    next_cell.ypos = tmp_ypos;
                    if (gradient1[2] + gradient2[2] >= 0)
                      next_cell.zpos = 0;
                    else
                      next_cell.zpos = 1;
                    next_cell.dirx_i = -1;
                    next_cell.diry_i = 0;
                    next_cell.dirz_i = 0;
                    next_cell.xno = cell.xno + 1;
                    next_cell.yno = cell.yno;
                    if (diff > 0)
                      next_cell.zno = cell.zno + cell.dirz_i;
                    else
                      next_cell.zno = cell.zno;

                    goto bestimmt;
                  }
              }
            if (grady < 0)
              {
                tmp_xpos = cell.xpos - cell.ypos * (gradx / grady);

                if (1 >= tmp_xpos && tmp_xpos >= 0)
                  {
                    /*Comparison of the gradients(in the neighboring cells) to determine the next cell:*/
                    /*Calculate the traveltime gradient*/
                    gradient1 = TimeGrad(cell.xno, cell.yno - 1, cell.zno + cell.dirz_i,
                        tt, ny, nz);
                    gradient2 = TimeGrad(cell.xno, cell.yno - 1, cell.zno, tt, ny, nz);
                    if (gradient1[2] + gradient2[2] >= 0)
                      diff = cell.dirz_i;
                    else
                      diff = -cell.dirz_i;

                    next_cell.xpos = tmp_xpos;
                    next_cell.ypos = 1;
                    if (gradient1[2] + gradient2[2] >= 0)
                      next_cell.zpos = 0;
                    else
                      next_cell.zpos = 1;
                    next_cell.dirx_i = 0;
                    next_cell.diry_i = 1;
                    next_cell.dirz_i = 0;
                    next_cell.xno = cell.xno;
                    next_cell.yno = cell.yno - 1;
                    if (diff > 0)
                      next_cell.zno = cell.zno + cell.dirz_i;
                    else
                      next_cell.zno = cell.zno;

                    goto bestimmt;
                  }
              }

            if (grady > 0)
              {
                tmp_xpos = cell.xpos - (cell.ypos - 1) * (gradx / grady);

                if (1 >= tmp_xpos && tmp_xpos >= 0)
                  {

                    /*Comparison of the gradients(in the neighboring cells) to determine the next cell:*/
                    /*Calculate the traveltime gradient*/
                    gradient1 = TimeGrad(cell.xno, cell.yno + 1, cell.zno + cell.dirz_i,
                        tt, ny, nz);
                    gradient2 = TimeGrad(cell.xno, cell.yno + 1, cell.zno, tt, ny, nz);
                    if (gradient1[2] + gradient2[2] >= 0)
                      {
                        diff = cell.dirz_i;
                      }
                    else
                      {
                        diff = -cell.dirz_i;
                      }

                    next_cell.xpos = tmp_xpos;
                    next_cell.ypos = 0;
                    if (gradient1[2] + gradient2[2] >= 0)
                      next_cell.zpos = 0;
                    else
                      next_cell.zpos = 1;
                    next_cell.dirx_i = 0;
                    next_cell.diry_i = -1;
                    next_cell.dirz_i = 0;

                    next_cell.xno = cell.xno;
                    next_cell.yno = cell.yno + 1;
                    if (diff > 0)
                      next_cell.zno = cell.zno + cell.dirz_i;
                    else
                      next_cell.zno = cell.zno;

                    goto bestimmt;
                  }
              }

          }

        /*Unexspected cases*/
        printf("Unexpected ray path in a grid cell\n");
        next_cell = cell;

        bestimmt: ;

        return (next_cell);
      }

    /*-------------------------------------------------------------*/
    /*Resort the ray segments; if a ray runs along the boundary of two cells this routine compares the velocity*/
    /*in the two cells and assigns the ray segment to the cell with the higher velocity*/
    /*Parameter:	*raypath			:= Raypath structures (Attention: x,y,z are normalized coordinates) */
    /*				data				:= Data structure*/
    /*				grid				:= Grid structure */

    int ResortRays(std::vector<RP_STRUCT> &raypath, const DATA_STRUCT &data,
        const GRID_STRUCT &grid)
      {
        long c, d, e;
        long ny, nz, nyz, ny1, nz1, nyz1;

        double eps = 0.01;

        ny = grid.ny;
        nz = grid.nz;
        nyz = ny * nz;

        ny1 = ny + 1;
        nz1 = nz + 1;
        nyz1 = ny1 * nz1;

        /*Loop over all rays*/
        for (size_t a = 0; a < data.ndata_seis; a++)
          /*Loop over all segments of a ray*/
          for (size_t b = 0; b < raypath[a].nray; b++)
            {
              /*Find the right cell*/
              c = (long) floor((double) (raypath[a].ele[b] / nyz));
              d = (long) floor((double) ((raypath[a].ele[b] - c * nyz) / ny));
              e = raypath[a].ele[b] - c * nyz - d * ny;
              if (raypath[a].ele[b] == (c * nyz + d * nz + e))
                {
                  if (raypath[a].x[b] - (double) (c) <= eps
                      && raypath[a].x[b + 1] - (double) c <= eps)
                    {
                      if (grid.slow[c * nyz1 + d * nz1 + e]
                          > grid.slow[(c - 1) * nyz1 + d * nz1 + e])
                        {
                          raypath[a].ele[b] = ((c - 1) * nyz + d * nz + e);
                        }
                    }

                  if (raypath[a].x[b] - (double) c >= 1.0 - eps
                      && raypath[a].x[b + 1] - (double) c >= 1.0 - eps)
                    {
                      if (grid.slow[c * nyz1 + d * nz1 + e]
                          > grid.slow[(c + 1) * nyz1 + d * nz1 + e])
                        {
                          raypath[a].ele[b] = ((c + 1) * nyz + d * nz + e);
                        }
                    }

                  if (raypath[a].y[b] - (double) d <= eps
                      && raypath[a].y[b + 1] - (double) d <= eps)
                    {
                      if (grid.slow[c * nyz1 + d * nz1 + e]
                          > grid.slow[c * nyz1 + (d - 1) * nz1 + e])
                        {
                          raypath[a].ele[b] = (c * nyz + (d - 1) * nz + e);
                        }
                    }

                  if (raypath[a].y[b] - (double) d >= 1.0 - eps
                      && raypath[a].y[b + 1] - (double) d >= 1.0 - eps)
                    {
                      if (grid.slow[c * nyz1 + d * nz1 + e]
                          > grid.slow[c * nyz1 + (d + 1) * nz1 + e])
                        {
                          raypath[a].ele[b] = (c * nyz + (d + 1) * nz + e);
                        }
                    }

                  if (raypath[a].z[b] - (double) e <= eps
                      && raypath[a].z[b + 1] - (double) e <= eps)
                    {
                      if (grid.slow[c * nyz1 + d * nz1 + e]
                          > grid.slow[c * nyz1 + d * nz1 + (e - 1)])
                        {
                          raypath[a].ele[b] = (c * nyz + d * nz + (e - 1));
                        }
                    }
                  if (raypath[a].z[b] - (double) e >= 1.0 - eps
                      && raypath[a].z[b + 1] - (double) e >= 1.0 - eps)
                    {
                      if (grid.slow[c * nyz1 + d * nz1 + e]
                          > grid.slow[c * nyz1 + d * nz1 + (e + 1)])
                        {
                          raypath[a].ele[b] = (c * nyz + d * nz + (e + 1));
                        }
                    }

                  goto slutt;
                }
              slutt: ;
            }

        return (1);
      }
  }

jif3D::RayResult ForwardModShot(const jif3D::GEOMETRY &geo,
    const jif3D::GRID_STRUCT &grid)
  {
    const float delta_num = (float) 0.001;
    const size_t nx3 = (grid.nx + 1);
    const size_t ny3 = (grid.ny + 1);
    const size_t nz3 = (grid.nz + 1);
    /*Allocate memory for the travel-times that will be calculated by the forward algorithm*/
    std::vector<float> tt(nx3 * ny3 * nz3, 0.0);
    size_t count; /*Number of active receivers for a shot*/
#ifdef HAVEHPX
//    std::cout << "Calculating shot " << i << " on node " << hpx::find_here() << " Thread: " << hpx::get_worker_thread_num() << "\n";
#endif
    float Xs, Ys, Zs;
    std::vector<float> Xr;
    std::vector<float> Yr;
    std::vector<float> Zr;
    jif3D::RayResult Result;

    Xs = ((geo.x[0]) / grid.h); /*normalized x-coordinate of the shot locations according to grid cell nodes*/
    Ys = ((geo.y[0]) / grid.h); /*normalized y-coordinate of the shot locations according to grid cell nodes*/
    Zs = ((geo.z[0]) / grid.h); /*normalized z-coordinate of the shot locations according to grid cell nodes*/

    /***************************************************************************************/
    /*Podvin&Lecomte forward algorithm*/
    /*tt is the calculated traveltime for each grid cell node*/

    std::vector<float> SlowBuffer(grid.slow);
    jif3D::PodvinTime3D().time_3d(&SlowBuffer[0], &tt[0], nx3, ny3, nz3, Xs, Ys, Zs,
        delta_num, 0);

    /***************************************************************************************/

    //jif3D::PlotTimeField("times.vtk", &tt[0], grid.h, nx3, ny3, nz3);
    /*Determine the receivers that are activate for the corresponding shot:*/
    count = geo.x.size() - 1;

    /*
     data.lshots[i] = count;*/
    Result.tcalc.resize(count);
    Xr.resize(count + 1);
    Yr.resize(count + 1);
    Zr.resize(count + 1);
    /***************************************************************************************/
    /*Determine the accurate traveltimes at the receiver-locations (by trilinear interpolation of the traveltimes at the grid cell edges)*/
    for (size_t j = 0; j < count; j++)
      {
        Xr[j] = ((geo.x[j + 1]) / grid.h); /*normalized x-coordinate of the receiver locations according to grid cell EDGES*/
        Yr[j] = ((geo.y[j + 1]) / grid.h); /*normalized y-coordinate of the receiver locations according to grid cell EDGES*/
        Zr[j] = ((geo.z[j + 1]) / grid.h); /*normalized z-coordinate of the receiver locations according to grid cell EDGES*/

        Result.tcalc[j] = interpolate(Xr[j], Yr[j], Zr[j], grid, &tt[0]);

      }

    /***************************************/

    /***************************************************************************************/
    /*Ray calculations (3-D version of the Aldridge&Oldenburg raypath-generation, 1993, Journal of seismic exploration,Vol.2,pages 257-274)*/

    /*Allocate temporary ray-structures for the specific shot*/
    Result.raypath.resize(count);

    /*Calculate the rays*/
    RayCalc(&tt[0], nx3, ny3, nz3, Xs, Ys, Zs, &Xr[0], &Yr[0], &Zr[0], count,
        Result.raypath);

    return Result;
  }

