//============================================================================
// Name        : X3DMTCalculator.h
// Author      : Jul 7, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================

#ifndef X3DMTCALCULATOR_H_
#define X3DMTCALCULATOR_H_

#include <limits>
#include <boost/filesystem.hpp>
#include <boost/serialization/serialization.hpp>
#include "MT3DCalculator.h"
#include "X3DModel.h"

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
    class X3DMTCalculator
      {
    public:
      //! This type definition is necessary so that ThreeDModelObjective can correctly deduce the native type for a model object for this class
      typedef X3DModel ModelType;
    private:
      //! The name of the executable for the x3d code
      std::string X3DName;
      //! The start of the names for files and directories created by this object
      std::string NameRoot;
      //! Remove all files created for running x3d
      void CleanUp();
      //! Make a unique string identifier for this object, basis for MakeUniqueName
      std::string ObjectID();
      //! Create a unique name for each object, calculation type and frequency so that we can write to different directories and execute in parallel
      std::string MakeUniqueName(X3DModel::ProblemType Type, const size_t FreqIndex);
      //! The directory to store all temporary files
      boost::filesystem::path TempDir;
      //! Do we want to perform distortion correction and calculate derivatives with respect to the distortion parameters
      bool WantDistCorr;
      //! The impedances from the last forward calculation without any distortion correction
      rvec RawImpedance;
      //! Calculate synthetic MT data for a single frequency
      rvec CalculateFrequency(const X3DModel &Model, const std::vector<double> &C,size_t freqindex);
      //! Calculate a least squares derivative for a single frequency
      rvec LQDerivativeFreq(const X3DModel &Model, const rvec &Misfit, const std::vector<double> &C, size_t freqindex);
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          // we do not need to serialize NameRoot, this is generated individually for each object
          if (Archive::is_saving::value)
            {
              std::string DirName(TempDir.string());
              ar & DirName;
              ar & WantDistCorr;
              ar & RawImpedance;
            }
          if (Archive::is_loading::value)
            {
              std::string DirName;
              ar & DirName;
              TempDir = DirName;
              ar & WantDistCorr;
              ar & RawImpedance;
            }

        }
    public:
      //! Given a conductivity model, calculate a vector of impedances
      /*! For a conductivity model given by the input parameter Model, we calculate the synthetic magnetotelluric data. When compiled with
       * an appropriate compiler the calculation is run in parallel for each frequency. We return the synthetic data as a real valued vector.
       * The ordering is \f$Re(Z_xx),Im(Z_xx),Re(Z_xy),\ldots,Im(Z_yy)\f$ for the first frequency for all sites, then second frequency for all sites etc.
       *
       * @param Model The description of the conductivity model including sites locations and frequencies.
       * @param minfreqindex The index of the first frequency for which to calculate the gradient
       * @param maxfreqindex The index one larger than the index of the last frequency for which to calculate the gradient (C++ loop convention)
       * @return The synthetic MT data in the format described above.
       */
      rvec Calculate(const ModelType &Model, size_t minfreqindex = 0,
          size_t maxfreqindex = std::numeric_limits<size_t>::max());
      //! Given a conductivity model and the misfit for each datum, calculate the derivative of the objective function with respect to the model parameters.
      /*! We use an adjoint approach to calculate the gradient of the objective functions with respect to the model parameters. As this approach requires
       * some of the fields from the forward calculation, the gradient will only be correct if the function Calculate of the same object has been called for
       * the same model beforehand. It is safe to calculate different models with separate objects between those calls.
       * @param Model The description of the conductivity model. Has to be the same as for the previous call to calculate.
       * @param Misfit The data misfit associated with the model.
       * @param minfreqindex The index of the first frequency for which to calculate the gradient
       * @param maxfreqindex The index one larger than the index of the last frequency for which to calculate the gradient (C++ loop convention)
       * @return The gradient of the objective function with respect to the model parameters for the given model. The storage ordering is identical to X3DModel.
       */
      rvec LQDerivative(const ModelType &Model, const rvec &Misfit, size_t minfreqindex =
          0, size_t maxfreqindex = std::numeric_limits<size_t>::max());
      //! The constructor takes optional arguments to change the directory were temporary files are stored and if we want to correct for distortion
      /*! When running calculations on a cluster in particular, it makes sense to store the files for x3D in a local
       * directory instead of the central filesystem. We can achieve this through setting TDir to an appropriate path.
       * If we set DC to true we calculate the derivative with respect to the distortion parameters.
       * @param TDir Directory to store the temporary files for x3D
       * @param x3d The name of the executable for the x3d code
       * @param DC Perform distortion correction and calculate the derivative with respect to the distortion parameters
       */
      X3DMTCalculator(boost::filesystem::path TDir = boost::filesystem::current_path(), std::string x3d = "x3d", bool DC = false);
      virtual ~X3DMTCalculator();
      };
  /* @} */
  }

#endif /* X3DMTCALCULATOR_H_ */
