//============================================================================
// Name        : ThreeDGravityFactory.h
// Author      : Jun 18, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#ifndef THREEDGRAVITYFACTORY_H_
#define THREEDGRAVITYFACTORY_H_

#include "ScalarOMPGravityImp.h"
#include "TensorOMPGravityImp.h"
#include "ThreeDGravityModel.h"
#include "../Global/Jif3DGlobal.h"

#ifdef HAVEGPU
#include "ScalarCudaGravityImp.h"
#include "TensorCudaGravityImp.h"
#else
#include "../Global/FatalException.h"
#endif

namespace jif3D
  {

    //! Create a new calculator object that is ready to use
    /*! This class provides a simple way to create a calculator object
     * with an appropriate implementation object without major user
     * interaction. This is the standard way of creating new calculator
     * classes. We implement it as a template, because the type of
     * calculator object should be known at compile time. Possible
     * template parameters are MinMemGravityCalculator, FullSensitivityGravityCalculator and
     * WaveletCompressedGravityCalculator.
     *
     * Returning boost::shared_ptr objects eases memory management for the user, as they
     * automatically delete their argument when no more shared pointers reference it.
     */
    template<class CalculatorClass>
    class J3DEXPORT CreateGravityCalculator
      {
    public:
      //! Make a new calculator object to calculate scalar gravity data
      static boost::shared_ptr<CalculatorClass> MakeScalar(bool wantcuda = false);
      //! Make a new calculator object to calculate FTG data
      static boost::shared_ptr<CalculatorClass> MakeTensor(bool wantcuda = false);
      };

    /*! Creates a new shared pointer to a calculator object for scalar gravity.We specify
     * the way the sensitivities are handled through the calculator class in the template parameter. \see ThreeDGravityCalculator
     * @param wantcuda Do we want to use NVidia CUDA for the calculation, might be ignored if no appropriate device found
     * @return A shared pointer to a calculator object
     */
    template<class CalculatorClass>
    boost::shared_ptr<CalculatorClass> CreateGravityCalculator<CalculatorClass>::MakeScalar(
        bool wantcuda)
      {
        boost::shared_ptr<ThreeDGravMagImplementation<ThreeDGravityModel> > Imp;
        if (wantcuda)
          {
#ifdef HAVEGPU
            Imp = boost::shared_ptr<ThreeDGravMagImplementation<ThreeDGravityModel> >(
                new ScalarCudaGravityImp);
#else
            throw jif3D::FatalException("Code has been compiled without GPU support ! ", __FILE__, __LINE__);
#endif
          }
        else
          {
            Imp = boost::shared_ptr<ThreeDGravMagImplementation<ThreeDGravityModel> >(
                new ScalarOMPGravityImp);
          }
        return boost::shared_ptr<CalculatorClass>(new CalculatorClass(Imp));
      }
    /*! Creates a new shared pointer to a calculator object for tensorial gravity. We specify
     * the way the sensitivities are handled through the calculator class in the template parameter. \see ThreeDGravityCalculator
     * @param wantcuda Do we want to use NVidia CUDA for the calculation, might be ignored if no appropriate device found
     * @return A shared pointer to a calculator object
     */
    template<class CalculatorClass>
    boost::shared_ptr<CalculatorClass> CreateGravityCalculator<CalculatorClass>::MakeTensor(
        bool wantcuda)
      {
        boost::shared_ptr<ThreeDGravMagImplementation<ThreeDGravityModel> > Imp;
        if (wantcuda)
          {
#ifdef HAVEGPU
            Imp = boost::shared_ptr<ThreeDGravMagImplementation<ThreeDGravityModel> >(
                new TensorCudaGravityImp);
#else
            throw jif3D::FatalException("Code has been compiled without GPU support ! ", __FILE__, __LINE__);
#endif
          }
        else
          {
            Imp = boost::shared_ptr<ThreeDGravMagImplementation<ThreeDGravityModel> >(new TensorOMPGravityImp);
          }
        return boost::shared_ptr<CalculatorClass>(new CalculatorClass(Imp));
      }

  }

#endif /* THREEDGRAVITYFACTORY_H_ */
