/*
 * X3DTipperCalculator.h
 *
 *  Created on: 25 May 2018
 *      Author: mm489
 */

#ifndef MT_X3DTIPPERCALCULATOR_H_
#define MT_X3DTIPPERCALCULATOR_H_

#include "../Global/Serialization.h"
#include "../Global/VecMat.h"

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/filesystem.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

#include "../Global/Jif3DGlobal.h"
#include "../Global/convert.h"
#include "../Global/Jif3DPlatformHelper.h"
#include "ReadWriteX3D.h"
#include "TipperData.h"
#include "X3DModel.h"
#include "X3DFieldCalculator.h"

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
    class J3DEXPORT X3DTipperCalculator
      {
    public:
      //! This type definition is necessary so that ThreeDModelObjective can correctly deduce the native type for a model object for this class
      typedef X3DModel ModelType;
      typedef TipperData DataType;

    private:
      //! The type of green's function for forward calculation in x3d (stage 1)
      jif3D::GreenCalcType GreenType1;
      //! The type of green's function for forward calculation in x3d (stage 4)
      jif3D::GreenCalcType GreenType4;
      //! The name of the executable for the x3d code
      std::string X3DName;
      //! The start of the names for files and directories created by this object
      std::string NameRoot;
      //! Remove all files created for running x3d
      void CleanUp();
      //! The directory to store all temporary files
      boost::filesystem::path TempDir;
      //! Do we want to delete the temporary files when object is destroyed
      bool CleanFiles;
      //! Store the previous execution time for different frequencies to optimize speed
      std::vector<std::pair<size_t, size_t>> ForwardExecTime;
      std::vector<std::pair<size_t, size_t>> DerivExecTime;
      friend class access;
      boost::shared_ptr<jif3D::X3DFieldCalculator> FieldCalculator;
      //create a unique ID that we can use to name things and still
      //perform parallel calculations
      std::string ObjectID()
        {
          //a unique ID created on construction
          boost::uuids::uuid tag = boost::uuids::random_generator()();
          //make a unique filename for the sensitivity file created by this object
          //we use boost uuid to generate a unique identifier tag
          //and translate it to a string to generate the filename
          return "mt" + jif3D::stringify(jif3D::platform::get_process_id()) + "x"
              + jif3D::stringify(this) + "t" + jif3D::stringify(tag);
        }

    public:
      //! Provide serialization to be able to store objects and, more importantly for hpx parallelization
      template<class Archive>
      void save(Archive & ar, const unsigned int version) const
        {
          ar & GreenType1;
          ar & GreenType4;
          ar & X3DName;
          std::string DirName(TempDir.string());
          ar & DirName;
          ar & CleanFiles;
          ar & ForwardExecTime;
        }
      template<class Archive>
      void load(Archive & ar, const unsigned int version)
        {
          ar & GreenType1;
          ar & GreenType4;
          ar & X3DName;
          std::string DirName;
          ar & DirName;
          TempDir = DirName;
          ar & CleanFiles;
          ar & ForwardExecTime;
        }
#ifdef HAVEHPX
      HPX_SERIALIZATION_SPLIT_MEMBER()
#else
      BOOST_SERIALIZATION_SPLIT_MEMBER()
#endif

      //! Set type of green's function for forward calculation i X3D (stage 1)
      void SetGreenType1(jif3D::GreenCalcType G)
        {
          GreenType1 = G;
        }
      //! Set type of green's function for forward calculation i X3D (stage 4)
      void SetGreenType4(jif3D::GreenCalcType G)
        {
          GreenType4 = G;
        }
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
      rvec Calculate(const ModelType &Model, const TipperData &Data, size_t minfreqindex =
          0, size_t maxfreqindex = std::numeric_limits<size_t>::max());
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
      rvec LQDerivative(const ModelType &Model, const TipperData &Data,
          const rvec &Misfit, size_t minfreqindex = 0, size_t maxfreqindex =
              std::numeric_limits<size_t>::max());

      //! The constructor takes optional arguments to change the directory were temporary files are stored and if we want to correct for distortion
      /*! When running calculations on a cluster in particular, it makes sense to store the files for x3D in a local
       * directory instead of the central filesystem. We can achieve this through setting TDir to an appropriate path.
       * If we set DC to true we calculate the derivative with respect to the distortion parameters.
       * @param TDir Directory to store the temporary files for x3D
       * @param x3d The name of the executable for the x3d code
       * @param Clean Delete all temporary files when object is destroyed
       */
      X3DTipperCalculator(boost::filesystem::path TDir = boost::filesystem::current_path(),
          std::string x3d = "x3d", bool Clean = true , boost::shared_ptr<jif3D::X3DFieldCalculator> FC = boost::make_shared<jif3D::X3DFieldCalculator> ());
      virtual ~X3DTipperCalculator();
      };
  /* @} */

  }
/* namespace jif3D */

#endif /* MT_X3DTIPPERCALCULATOR_H_ */
