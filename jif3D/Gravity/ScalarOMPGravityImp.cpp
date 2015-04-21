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
#include <hpx/parallel/execution_policy.hpp>
#endif

#include <numeric>
#include <functional>

#include "ScalarOMPGravityImp.h"
#include "BasicGravElements.h"
#include "GravityBackground.h"

#ifdef HAVEHPX
std::vector<double> GravityChunk(size_t start, size_t end, double x_meas, double y_meas,
    double z_meas, const std::vector<double> &XCoord, const std::vector<double> &YCoord,
    const std::vector<double> &ZCoord, const std::vector<double> &XSizes,
    const std::vector<double> &YSizes, const std::vector<double> &ZSizes, size_t ysize,
    size_t zsize);

HPX_PLAIN_ACTION(GravityChunk, Gravity_Action)

std::vector<double> GravityChunk(size_t start, size_t end, double x_meas, double y_meas,
    double z_meas, const std::vector<double> &XCoord,
    const std::vector<double> &YCoord,
    const std::vector<double> &ZCoord,
    const std::vector<double> &XSizes,
    const std::vector<double> &YSizes,
    const std::vector<double> &ZSizes, size_t ysize, size_t zsize)
  {
    if (end > start)
      {
        std::vector<double> Result;
        Result.reserve(end - start);
        for (size_t i = start; i < end; ++i)
          {
            size_t zi = i % zsize;
            size_t xi = (i - zi) / zsize;
            size_t yi = xi % ysize;
            xi = (xi - yi) / ysize;
            double c = jif3D::CalcGravBoxTerm(x_meas, y_meas, z_meas, XCoord[xi], YCoord[yi],
                ZCoord[zi], XSizes[xi], YSizes[yi], ZSizes[zi]);
            Result.push_back(c);
          }
        return Result;
      }
    return std::vector<double>();
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
        rmat &Sensitivities)
      {
        return CalcScalarBackground(measindex, xwidth, ywidth, zwidth, Model,
            Sensitivities);
      }

    /*! Calculate the gravitational effect of the 3D model at a single measurement site.
     * @param measindex The index of the measurement
     * @param Model The gravity model
     * @param Sensitivities If the matrix hold enough elements the sensitivities for the background are stored for a single measurement
     * @return The gravitational acceleration in m/s^2 due to the model at this site
     */
    rvec ScalarOMPGravityImp::CalcGridded(const size_t measindex,
        const ThreeDGravityModel &Model, rmat &Sensitivities)
      {
        //get the dimensions of the model
        const size_t xsize = Model.GetDensities().shape()[0];
        const size_t ysize = Model.GetDensities().shape()[1];
        const size_t zsize = Model.GetDensities().shape()[2];
        const double x_meas = Model.GetMeasPosX()[measindex];
        const double y_meas = Model.GetMeasPosY()[measindex];
        const double z_meas = Model.GetMeasPosZ()[measindex];
        const int nmod = xsize * ysize * zsize;
        const bool storesens = (Sensitivities.size1() >= ndatapermeas)
            && (Sensitivities.size2() >= size_t(nmod));
        double returnvalue = 0.0;
        double currvalue = 0.0;
#ifdef HAVEOPENMP
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
                    XCoord[xindex], YCoord[yindex], ZCoord[zindex],
                    XSizes[xindex], YSizes[yindex], ZSizes[zindex]);
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

        using hpx::lcos::future;
        using hpx::async;
        using hpx::wait_all;
        std::vector<hpx::naming::id_type> localities = hpx::find_all_localities();
        const size_t nthreads = hpx::get_num_worker_threads();
        const size_t nlocs = localities.size();
        //BOOST_LOG_TRIVIAL(debug)<< "Running on: " << nlocs << " localities. With " << nthreads << " worker threads " << std::endl;
        size_t nchunks = nthreads * nlocs;
        std::vector<hpx::lcos::future<std::vector<double>>> ChunkResult;
        Gravity_Action GravityChunks;
        const size_t cellsperchunk = nmod / nchunks +1;
        std::vector<double> XC(XCoord.data(),XCoord.data() + XCoord.num_elements());
        std::vector<double> YC(YCoord.data(),YCoord.data() + YCoord.num_elements());
        std::vector<double> ZC(ZCoord.data(),ZCoord.data() + ZCoord.num_elements());
        std::vector<double> XS(XSizes.data(),XSizes.data() + XSizes.num_elements());
        std::vector<double> YS(YSizes.data(),YSizes.data() + YSizes.num_elements());
        std::vector<double> ZS(ZSizes.data(),ZSizes.data() + ZSizes.num_elements());
        for (size_t c = 0; c < nchunks; ++c)
          {
            size_t count = 0;

            size_t startindex = c * cellsperchunk;
            size_t endindex = std::min(size_t(nmod), (c+1) * cellsperchunk);

            hpx::naming::id_type const locality_id = localities.at(c % localities.size());
            ChunkResult.push_back(async(GravityChunks, locality_id, startindex,endindex,x_meas,y_meas,z_meas,
                    XC, YC, ZC,
                    XS, YS, ZS, ysize, zsize));
          }
        wait_all(ChunkResult);
        std::vector<double> Sens;
        Sens.reserve(nmod);
        for (size_t c = 0; c < nchunks; ++c)
          {
            std::vector<double> CurrSens = ChunkResult.at(c).get();
            std::copy(CurrSens.begin(),CurrSens.end(),back_inserter(Sens));
          }
        for (int offset = 0; offset < nmod; ++offset)
          {
            //we store the current value for possible sensitivity calculations
            //currvalue contains the geometric term, i.e. the sensitivity

            returnvalue += Sens[offset] * Model.GetDensities().data()[offset];
            if (storesens)
              {
                Sensitivities(0, offset) = currvalue;
              }
          }

#endif
        rvec returnvector(1);
        returnvector(0) = returnvalue;
        return returnvector;

      }
  }
