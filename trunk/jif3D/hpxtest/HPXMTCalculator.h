//============================================================================
// Name        : X3DMTCalculator.h
// Author      : Jul 7, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================

#ifndef HPXMTCALCULATOR_H_
#define HPXMTCALCULATOR_H_

#include <hpx/config.hpp>
//#include <hpx/include/util.hpp>
#include <limits>
#include <boost/filesystem.hpp>
#include <boost/serialization/serialization.hpp>
#include "../MT/X3DModel.h"
#include "../Global/VecMat.h"



namespace jif3D
  {
    /** \addtogroup mtmodelling Forward modelling of magnetotelluric data */
    /* @{ */
    //! Calculate magnetotelluric data from a 3D conductivity model using x3D by Avdeev et al.
    /*! This class uses x3d by Avdeev et al., see
     * Avdeev, D.B., A.V. Kuvshinov, O.V. Pankratov, and G. A. Newman, 1997, High-performance three-dimensional electromagnetic modelling using
     * modified Neumann series. Wide-band numerical solution and examples,  J. Geomagn. Geoelectr., 49, 1519-1539,
     * to calculate synthetic magnetotelluric data. We also implement the adjoint approach described in
     * Avdeev and Avdeeva, 3D magnetotelluric inversion using a limited-memory quasi-Newton optimization, Geophysics, 74(3), 2009
     * to calculate the gradient of the objective function.
     *
     * As we do not have access to the source code, we interface with the code by writing out files in the format required by x3d,
     * run the executable and then read in the results in ascii format. This is not the most efficient way of communication, but given
     * typical run times on the order of minutes, the relative overhead associated with this is low.
     *
     * There are a few restrictions on the conductivity model in order to ensure correct calculation of gradients. x3d
     * performs some "optimizations" that interfere with our approach for gradient calculation. Two adjacent layers
     * in the grid cannot have absolutely identical values, as x3d will then treat them as a single layer. Therefore at
     * least one grid cell must have a slightly different value from the layer above and below. Similarly a homogeneous
     * layer (with all grid cell values identical) cannot match the value for the background layer at this depth. In both
     * cases x3d will not yield field values for these areas in the .ema files that we need to calculate the gradient.
     * As the difference does not have to be very large, it is probably best to always have one grid cell in each layer
     * that differs from everything else by 0.1% or so.
     */
    class HPXMTCalculator
      {
    public:
      //! This type definition is necessary so that ThreeDModelObjective can correctly deduce the native type for a model object for this class
      typedef X3DModel ModelType;
    private:
      //! Remove all files created for running x3d
      void CleanUp();
      //! The directory to store all temporary files
      std::string TempDirName;
      rvec LQDerivativeFreq(const X3DModel &Model, const rvec &Misfit, size_t freqindex);

    public:
      //! Given a conductivity model, calculate a vector of impedances
      /*! For a conductivity model given by the input parameter Model, we calculate the synthetic magnetotelluric data. When compiled with
       * an appropriate compiler the calculation is run in parallel for each frequency. We return the synthetic data as a real valued vector.
       * The ordering is \f$Re(Z_xx),Im(Z_xx),Re(Z_xy),\ldots,Im(Z_yy)\f$ for the first frequency for all sites, then second frequency for all sites etc.
       *
       * @param Model The description of the conductivity model including sites locations and frequencies.
       * @return The synthetic MT data in the format described above.
       */
      rvec Calculate(const ModelType &Model, size_t minfreqindex = 0,
          size_t maxfreqindex = std::numeric_limits<size_t>::max());
      //! The constructor takes an optional argument to change the directory were temporary files are stored
      HPXMTCalculator(boost::filesystem::path TDir = boost::filesystem::current_path());
      virtual ~HPXMTCalculator();
      };
  /* @} */
  }

#endif /* HPXMTCALCULATOR_H_ */