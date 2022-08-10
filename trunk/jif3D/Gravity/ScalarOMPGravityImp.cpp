//============================================================================
// Name        : ScalarOMPGravityImp.cpp
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

#include <numeric>
#include <functional>

#include "ScalarOMPGravityImp.h"
#include "BasicGravElements.h"
#include "GravityBackground.h"

#ifdef HAVEHPX
/*!  Calculate the contribution of several prisms to the gravitational acceleration
 * at a single site.
 * @param start The index of the first cell in the grid to be used for the calculation
 * @param end The index of the last cell in the grid to be used for the calculation
 * @param x_meas The x-coordinate of the measurement site in m
 * @param y_meas The y-coordinate of the measurement site in m
 * @param z_meas The z-coordinate of the measurement site in m
 * @param XCoord The X-Coordinates of the grid cell boundaries in m
 * @param YCoord The Y-Coordinates of the grid cell boundaries in m
 * @param ZCoord The Z-Coordinates of the grid cell boundaries in m
 * @param XSizes The size of the grid cell in x-direction in m
 * @param YSizes The size of the grid cell in y-direction in m 
 * @param ZSizes The size of the grid cell in z-direction in m
 * @return The vector of geometric factors for each of the grid cells.
 */
std::vector<double> GravityChunk(size_t start, size_t end, double x_meas, double y_meas,
    double z_meas, const std::vector<double> &XCoord, const std::vector<double> &YCoord,
    const std::vector<double> &ZCoord, const std::vector<double> &XSizes,
    const std::vector<double> &YSizes, const std::vector<double> &ZSizes);
// Define the corresponding hpx action so that we can use GravityChunk in parallel with hpx
HPX_PLAIN_ACTION(GravityChunk, Gravity_Action)

std::vector<double> GravityChunk(size_t start, size_t end, double x_meas, double y_meas,
    double z_meas, const std::vector<double> &XCoord,
    const std::vector<double> &YCoord,
    const std::vector<double> &ZCoord,
    const std::vector<double> &XSizes,
    const std::vector<double> &YSizes,
    const std::vector<double> &ZSizes
    )
  {
    //The basis for the calculation is a rectilinear grid of cells
    //For these calculations we do not need the density for each of the grid cells
    //but only the geometry of the cells and the distance of the cell corners from
    //the measurement site
    // There is no interaction between the calculations for each prism, so we can divide
    // it up arbitrarily, the easiest way, seems to be to linearly address all cells
    // with a storage order z,y,x. This is compatible with the storage order for densities
    //that we need for calculation elsewhere
    //The number of prisms in y-direction

    //const size_t ysize = YCoord.size();
    const size_t ysize = YSizes.size();
    //The number of prisms in z-direction
    //const size_t zsize = ZCoord.size();
    const size_t zsize = ZSizes.size();
    //we will generate end- start geometric factors that we will return
    std::vector<double> Result;
    // Record the start and end position of the job allocated, which are used in the wait_each process
    Result.push_back(start);
    Result.push_back(end);

    //go through all the elements included in this chunk
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
        double c = jif3D::CalcGravBoxTerm(x_meas, y_meas, z_meas, XCoord[xi], YCoord[yi],
            ZCoord[zi], XSizes[xi], YSizes[yi], ZSizes[zi]);
        //store result in vector
        Result.push_back(c);
        }
    return Result;
  }
#endif

namespace jif3D
{
    ScalarOMPGravityImp::ScalarOMPGravityImp()
      {

      }

    ScalarOMPGravityImp::~ScalarOMPGravityImp()
      {

      }

    /*!  Calculate the contribution of a layered background to a scalar gravity measurement.
     * @param measindex the index of the measurement
     * @param xwidth The total width of the discretized model area in x-direction in m
     * @param ywidth The total width of the discretized model area in y-direction in m
     * @param zwidth The total width of the discretized model area in z-direction in m
     * @param Model The gravity model
     * @param Sensitivities If the matrix hold enough elements the sensitivities for the background are stored for a single measurement
     * @return The gravitational acceleration in m/s^2 due to the background
     */
    rvec ScalarOMPGravityImp::CalcBackground(const size_t measindex, const double xwidth,
        const double ywidth, const double zwidth, const ThreeDGravityModel &Model,
        const ScalarGravityData &Data, rmat &Sensitivities)
      {
        return CalcScalarBackground(measindex, xwidth, ywidth, zwidth, Model, Data,
            Sensitivities);
      }

    /*! Calculate the gravitational effect of the 3D model at a single measurement site.
     * @param measindex The index of the measurement
     * @param Model The gravity model
     * @param Sensitivities If the matrix hold enough elements the sensitivities for the background are stored for a single measurement
     * @return The gravitational acceleration in m/s^2 due to the model at this site
     */
    rvec ScalarOMPGravityImp::CalcGridded(const size_t measindex,
        const ThreeDGravityModel &Model, const ScalarGravityData &Data,
        rmat &Sensitivities)
      {
        //get the dimensions of the model
        const size_t xsize = Model.GetDensities().shape()[0];
        const size_t ysize = Model.GetDensities().shape()[1];
        const size_t zsize = Model.GetDensities().shape()[2];
        const double x_meas = Data.GetMeasPosX().at(measindex);
        const double y_meas = Data.GetMeasPosY().at(measindex);
        const double z_meas = Data.GetMeasPosZ().at(measindex);
        const int nmod = xsize * ysize * zsize;
        const bool storesens = (Sensitivities.size1() >= ndatapermeas)
            && (Sensitivities.size2() >= size_t(nmod));
        double returnvalue = 0.0;
#ifdef HAVEOPENMP
        double currvalue = 0.0;
        //sum up the contributions of all prisms in an openmp parallel loop
#pragma omp parallel default(shared) private(currvalue) reduction(+:returnvalue)
          {
            //instead of nested loops over each dimension, we have one big
            //loop over all elements, this allows for a nearly infinite number
            //of parallel processors
#pragma omp for
            for (int offset = 0; offset < nmod; ++offset)
              {

                //we store the current value for possible sensitivity calculations
                //currvalue contains the geometric term, i.e. the sensitivity
                int xindex, yindex, zindex;
                Model.OffsetToIndex(offset, xindex, yindex, zindex);
                currvalue = CalcGravBoxTerm(x_meas, y_meas, z_meas,
                    Model.GetXCoordinates()[xindex], Model.GetYCoordinates()[yindex], Model.GetZCoordinates()[zindex],
                    Model.GetXCellSizes()[xindex], Model.GetYCellSizes()[yindex], Model.GetZCellSizes()[zindex]);
                returnvalue += currvalue
                * Model.GetDensities()[xindex][yindex][zindex];
                if (storesens)
                  {
                    Sensitivities(0, offset) = currvalue;
                  }
              }
          } //end of parallel section 
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

        Gravity_Action GravityChunks;
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
            ChunkResult.push_back(async(GravityChunks, locality_id, startindex,endindex,x_meas,y_meas,z_meas,
                Model.GetXCoordinates(), Model.GetYCoordinates(), Model.GetZCoordinates(),
                Model.GetXCellSizes(), Model.GetYCellSizes(), Model.GetZCellSizes()));
          }

        // Create Sens before wait_all to save some time
        std::vector<double> Sens(nmod);
        std::vector<double> Densities;

        // Store the densities to calculate returnvalue in parallel
        for(int offset=0;offset<nmod;offset++)
        {
            Densities.push_back(Model.GetDensities().data()[offset]);
        }

        // When each future finished, call this function to record the value to Sens
        // Set two situations if we need to store the value calculated

        auto f = [&Sens, &Sensitivities, &storesens](std::size_t i, hpx::future<std::vector<double>>&& f) {
            std::vector<double> CurrSens = f.get();
            // Copy the result to Sens from the start position to the end position
            int copy_index = 2;
            if (storesens) 
            {
                for (int t = CurrSens[0]; t < CurrSens[1]; t++)
                {
                    Sens[t] = CurrSens[copy_index];
                    Sensitivities(t) = Sens[t];
                    copy_index++;
                }
            }
            else 
            {
                for (int t = CurrSens[0]; t < CurrSens[1]; t++)
                {
                    Sens[t] = CurrSens[copy_index];
                    copy_index++;
                }
            }
        };

        // When a process is finished, copy it to the correct position
        hpx::wait_each(f, ChunkResult);

        // Use parallel transform and resuce to get the return value
        returnvalue = hpx::transform_reduce(hpx::execution::parallel_policy(),
            Sens.begin(), Sens.end(),
            Densities.begin(),
            0.0L,
            std::plus<>(),
            [](double Sen_val, double Den_val) {
                return Sen_val * Den_val;
            });

#endif

        rvec returnvector(1);
        returnvector(0) = returnvalue;
        return returnvector;
      }
  }
