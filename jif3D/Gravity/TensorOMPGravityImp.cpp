//============================================================================
// Name        : TensorOMPGravityImp.cpp
// Author      : Dec 10, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================

#ifdef HAVEHPX
#include <hpx/hpx_init.hpp>
#include <hpx/include/actions.hpp>
#include <hpx/include/util.hpp>
#include <hpx/include/lcos.hpp>
#include <hpx/parallel/algorithms/transform_reduce.hpp>
#include <hpx/runtime_distributed.hpp>
#include <hpx/parallel/algorithms/transform.hpp>
#include <hpx/execution.hpp>
#include <parallel/algorithm>
#include <hpx/iostream.hpp>

#endif
#include "TensorOMPGravityImp.h"
#include "BasicGravElements.h"
#include "GravityBackground.h"

#ifdef HAVEHPX
std::vector<double> GravityChunkten(std::size_t start, std::size_t end, double x_meas, double y_meas,
    double z_meas, const std::vector<double>& XCoord, const std::vector<double>& YCoord,
    const std::vector<double>& ZCoord, const std::vector<double>& XSizes,
    const std::vector<double>& YSizes, const std::vector<double>& ZSizes,
    const std::vector<double>& Densities);
// Define the corresponding hpx action so that we can use GravityChunk in parallel with hpx
HPX_PLAIN_ACTION(GravityChunkten, Gravityten_Action)

std::vector<double> GravityChunkten(std::size_t start, std::size_t end, double x_meas, double y_meas,
    double z_meas, const std::vector<double>& XCoord,
    const std::vector<double>& YCoord,
    const std::vector<double>& ZCoord,
    const std::vector<double>& XSizes,
    const std::vector<double>& YSizes,
    const std::vector<double>& ZSizes,
    const std::vector<double>& Densities)
{
    //const size_t ysize = YCoord.size();
    const std::size_t ysize = YSizes.size();
    //The number of prisms in z-direction
    //const size_t zsize = ZCoord.size();
    const std::size_t zsize = ZSizes.size();
    //we will generate end- start geometric factors that we will return
    std::vector<double> Result;
    // Set the result value to return
    double U0 = 0.0, U1 = 0.0, U2 = 0.0, U3 = 0.0, U4 = 0.0, U5 = 0.0, U6 = 0.0, U7 =
        0.0, U8 = 0.0;
    for (size_t i = start; i < end; ++i)
    {
        //calculate the indices the current element would correspond to 
        //in a grid with storage order z,y,x
        size_t zi = i % zsize;
        size_t xi = (i - zi) / zsize;
        size_t yi = xi % ysize;
        xi = (xi - yi) / ysize;

        //calculate the geometric factor based on the position of the measurement
        //sites and the parameters of the currrent grid cells
        jif3D::GravimetryMatrix c = jif3D::CalcTensorBoxTerm(x_meas, y_meas, z_meas, XCoord[xi], YCoord[yi],
            ZCoord[zi], XSizes[xi], YSizes[yi], ZSizes[zi]);
        U0 += c(0, 0) * Densities[i];
        U1 += c(0, 1) * Densities[i];
        U2 += c(0, 2) * Densities[i];
        U3 += c(1, 0) * Densities[i];
        U4 += c(1, 1) * Densities[i];
        U5 += c(1, 2) * Densities[i];
        U6 += c(2, 0) * Densities[i];
        U7 += c(2, 1) * Densities[i];
        U8 += c(2, 2) * Densities[i];
        //store result in vector
    }
    Result.push_back(U0);
    Result.push_back(U1);
    Result.push_back(U2);
    Result.push_back(U3);
    Result.push_back(U4);
    Result.push_back(U5);
    Result.push_back(U6);
    Result.push_back(U7);
    Result.push_back(U8);
    return Result;
}
#endif

namespace jif3D
  {

    TensorOMPGravityImp::TensorOMPGravityImp()
      {

      }

    TensorOMPGravityImp::~TensorOMPGravityImp()
      {

      }

    /*!  Calculate the contribution of a layered background to a tensor gravity measurement.
     * @param measindex The index of the measurement
     * @param xwidth The total width of the discretized model area in x-direction in m
     * @param ywidth The total width of the discretized model area in y-direction in m
     * @param zwidth The total width of the discretized model area in z-direction in m
     * @param Model The gravity model
     * @param Sensitivities The \f$ 9 \times m\f$ matrix of sensitivities for the current measurement
     * @return The gravitational tensor due to the background
     */
    rvec TensorOMPGravityImp::CalcBackground(const size_t measindex, const double xwidth,
        const double ywidth, const double zwidth, const ThreeDGravityModel &Model,
        const TensorGravityData &Data, rmat &Sensitivities)
      {
        return CalcTensorBackground(measindex, xwidth, ywidth, zwidth, Model, Data,
            Sensitivities);
      }
    /*! Calculate the FTG response of the gridded domain.
     * @param measindex The index of the measurement
     * @param Model The gravity model
     * @param Sensitivities The \f$ 9 \times m\f$ matrix of sensitivities for the current measurement
     * @return A 9 component vector with the FTG matrix components
     */
    rvec TensorOMPGravityImp::CalcGridded(const size_t measindex,
        const ThreeDGravityModel &Model, const TensorGravityData &Data,
        rmat &Sensitivities)
      {
        const size_t xsize = Model.GetDensities().shape()[0];
        const size_t ysize = Model.GetDensities().shape()[1];
        const size_t zsize = Model.GetDensities().shape()[2];
        const double x_meas = Data.GetMeasPosX()[measindex];
        const double y_meas = Data.GetMeasPosY()[measindex];
        const double z_meas = Data.GetMeasPosZ()[measindex];
        const int nmod = xsize * ysize * zsize;
        const bool storesens = (Sensitivities.size1() >= ndatapermeas)
            && (Sensitivities.size2() >= size_t(nmod));
        GravimetryMatrix currvalue(3, 3);
        double U0 = 0.0, U1 = 0.0, U2 = 0.0, U3 = 0.0, U4 = 0.0, U5 = 0.0, U6 = 0.0, U7 =
            0.0, U8 = 0.0;
        //we cannot add up a user defined quantity in parallel
        //so break up the tensor into its component with different variables
        //and assign the results after the parallel loop
        //sum up the contributions of all prisms
#ifdef HAVEOPENMP
#pragma omp parallel default(shared) private(currvalue) reduction(+:U0,U1,U2,U3,U4,U5,U6,U7,U8)
          {
            //instead of nested loops over each dimension, we have one big
            //loop over all elements, this allows for a nearly infinite number
            //of parallel processors
#pragma omp for
            for (int offset = 0; offset < nmod; ++offset)
              {
                int xindex, yindex, zindex;
                //we still need the indices for each dimension
                //so we have to convert our loop variable
                Model.OffsetToIndex(offset, xindex, yindex, zindex);
                //currvalue contains only the geometric term
                currvalue = CalcTensorBoxTerm(x_meas, y_meas, z_meas,
                    Model.GetXCoordinates()[xindex], Model.GetYCoordinates()[yindex],
                    Model.GetZCoordinates()[zindex], Model.GetXCellSizes()[xindex],
                    Model.GetYCellSizes()[yindex], Model.GetZCellSizes()[zindex]);
                //to we have to multiply each element by the density
                const double Density = Model.GetDensities()[xindex][yindex][zindex];
                U0 += currvalue(0, 0) * Density;
                U1 += currvalue(0, 1) * Density;
                U2 += currvalue(0, 2) * Density;
                U3 += currvalue(1, 0) * Density;
                U4 += currvalue(1, 1) * Density;
                U5 += currvalue(1, 2) * Density;
                U6 += currvalue(2, 0) * Density;
                U7 += currvalue(2, 1) * Density;
                U8 += currvalue(2, 2) * Density;
                if (storesens)
                  {
                    for (size_t i = 0; i < ndatapermeas; ++i)
                      Sensitivities(i, offset) = currvalue.data()[i];
                  }
              }
          } //end of parallel region

#endif

#ifdef HAVEHPX
        using hpx::future;
        using hpx::async;
        std::vector<hpx::id_type> localities = hpx::find_all_localities();
        const size_t nthreads = hpx::get_num_worker_threads();
        const size_t nlocs = localities.size();

        //calculate how many chunks we have to divide the work into
        //assuming 1 chunk per thread and node, this might need adjustment
        size_t nchunks = nthreads * nlocs;

        // Here I change the way of dividing the work, initializing the ChunkResult so that it can work in the correct order
        std::vector<hpx::future<std::vector<double>>> ChunkResult;

        // Initial the density value used to 
        std::vector<double> Densities;
        for (int i = 0; i < nmod; i++) 
        {
            Densities.push_back(Model.GetDensities().data()[i]);
        }

        Gravityten_Action GravityChunktens;
        size_t cellsperchunk = 0;
        if (nmod % nchunks == 0)
        {
            cellsperchunk = nmod / nchunks;
    }
        else
        {
            cellsperchunk = nmod / nchunks + 1;
        }

        //create the individual chunks of work
        for (size_t c = 0; c < nchunks; ++c)
        {
            if (c * cellsperchunk >= size_t(nmod))
            {
                break;
            }
            //find the address of the first element of the current chunk
            size_t startindex = c * cellsperchunk;
            //the last element of the chunk, considering that we the last chunk might be a bit too big
            size_t endindex = std::min(size_t(nmod), (c + 1) * cellsperchunk);
            // simple round robin allocation to different localities, not sure this makes much sense
            hpx::id_type locality_id;

            locality_id = localities.at(c % localities.size());
            //create a hpx future
            ChunkResult.push_back(async(GravityChunktens, locality_id, startindex, endindex, x_meas, y_meas, z_meas,
                Model.GetXCoordinates(), Model.GetYCoordinates(), Model.GetZCoordinates(),
                Model.GetXCellSizes(), Model.GetYCellSizes(), Model.GetZCellSizes(), Densities));
        }

        auto add_function = [&U0, &U1, &U2, &U3, &U4, &U5, &U6, &U7, &U8](std::size_t i, hpx::future<std::vector<double>>&& f) {
            std::vector<double> Temp = f.get();
            U0 += Temp[0];
            U1 += Temp[1];
            U2 += Temp[2];
            U3 += Temp[3];
            U4 += Temp[4];
            U5 += Temp[5];
            U6 += Temp[6];
            U7 += Temp[7];
            U8 += Temp[8];
            // Add the result to final U0 to U8
        };

        // When a process is finished, copy it to the correct position
        hpx::wait_each(add_function, ChunkResult);

#endif
        rvec returnvalue(ndatapermeas);
        returnvalue(0) = U0;
        returnvalue(1) = U1;
        returnvalue(2) = U2;
        returnvalue(3) = U3;
        returnvalue(4) = U4;
        returnvalue(5) = U5;
        returnvalue(6) = U6;
        returnvalue(7) = U7;
        returnvalue(8) = U8;

        return returnvalue;
      }
  }
