//============================================================================
// Name        : ThreeDGravityCalculator.h
// Author      : Dec 10, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#ifndef THREEDGRAVITYCALCULATOR_H_
#define THREEDGRAVITYCALCULATOR_H_

#include "ThreeDGravityModel.h"
#include "ThreeDGravityImplementation.h"
#include "ScalarOMPGravityImp.h"
#include "ScalarCudaGravityImp.h"
#include "TensorOMPGravityImp.h"
#include <boost/shared_ptr.hpp>

namespace jiba
  {
    /** \addtogroup gravity Gravity forward modeling, display and inversion */
    /* @{ */
    //! The base class for all calculator objects, these handle how the sensitivities are stored and processed
    /*! The ThreeDGravityCalculator class is the base class that provides
     * the user interface to the gravity forward calculation. It uses a
     * ThreeDGravityImplementation object to do the actual forward calculation.
     * The derived classes determine how the sensitivity information
     * is handled, whether it is discarded to save memory, fully stored, or
     * compressed in some way. Also they decide in how far calculations are
     * cached using sensitivity information.
     *
     * Within the forward calculation for each measurement the implementation class
     * calls the member function HandleSensitivities() with the index for the measurement
     * to give the calculator class the possibility to store, process or discard the
     * current row of the sensitivity matrix. This avoids having to store the complete matrix
     * in cases where some compression is applied.
     */
    class ThreeDGravityCalculator
      {
    private:
      // We need a structure to hold the sensitivities for
      // the current measurement that can be passed to
      // the implementation object. The derived classes
      // decide how to handle this information
      rmat CurrentSensitivities;
    protected:
      // The shared pointer to the implementation object
      // that does the actual calculation
      boost::shared_ptr<ThreeDGravityImplementation> Imp;
      // check the information in the model is consistent
      void CheckModelConsistency(const ThreeDGravityModel &Model);
    public:
      //! Calculate the forward response of the given model, has to be implemented in a derived class
      virtual rvec Calculate(const ThreeDGravityModel &Model) = 0;
      //! Read and write access to the sensitivity information for the current measurement, only intended for implementation classes
      rmat &SetCurrentSensitivities()
        {
          return CurrentSensitivities;
        }
      //! Process the sensitivity information for the current measurement, called from the implementation class
      virtual void HandleSensitivities(const size_t measindex) = 0;
      //! This class is useless without an implementation object so we have to pass one to the constructor
      ThreeDGravityCalculator(boost::shared_ptr<ThreeDGravityImplementation> TheImp);
      virtual ~ThreeDGravityCalculator();
      };
    //! Create a new calculator object that is ready to use
    /*! This class provides a simple way to create a calculator object
     * with an appropriate implementation object without major user
     * interaction. This is the standard way of creating new calculator
     * classes. We implement it as a template, because the type of
     * calculator object should be known at compile time.
     */
    template<class CalculatorClass>
    class CreateGravityCalculator
      {
    public:
      typedef boost::shared_ptr<CalculatorClass> sp_CalculatorClass;
      //! Make a new calculator object to calculate scalar gravity data
      static boost::shared_ptr<CalculatorClass> MakeScalar(bool wantcuda =
          false);
      //! Make a new calculator object to calculate FTG data
      static boost::shared_ptr<CalculatorClass> MakeTensor();
      };


    /*! Creates a new shared pointer to a calculator object for scalar gravity of the type specified in
     * the template parameter.
     * @param wantcuda Do we want to use NVidia CUDA for the calculation, might be ignored if no appropriate device found
     * @return A shared pointer to a calculator object
     */
    template<class CalculatorClass>
    boost::shared_ptr<CalculatorClass> CreateGravityCalculator<CalculatorClass>::MakeScalar(
        bool wantcuda)
      {
        boost::shared_ptr<ThreeDGravityImplementation> Imp;
        if (wantcuda)
          {
            Imp = boost::shared_ptr<ThreeDGravityImplementation>(
                new ScalarCudaGravityImp);
          }
        else
          {
            Imp = boost::shared_ptr<ThreeDGravityImplementation>(
                new ScalarOMPGravityImp);
          }
        return boost::shared_ptr<CalculatorClass>(new CalculatorClass(Imp));
      }
    /*! Creates a new shared pointer to a calculator object for tensorial gravity of the type specified in
     * the template parameter.
     * @return A shared pointer to a calculator object
     */
    template<class CalculatorClass>
    boost::shared_ptr<CalculatorClass> CreateGravityCalculator<CalculatorClass>::MakeTensor()
      {
        boost::shared_ptr<ThreeDGravityImplementation> Imp(
            new TensorOMPGravityImp);
        return boost::shared_ptr<CalculatorClass>(new CalculatorClass(Imp));
      }

  /* @} */
  }

#endif /* THREEDGRAVITYCALCULATOR_H_ */
