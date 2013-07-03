//============================================================================
// Name        : RegularizationFunction.h
// Author      : 1 Mar 2013
// Version     : 
// Copyright   : 2013, mm489
//============================================================================

#ifndef REGULARIZATIONFUNCTION_H_
#define REGULARIZATIONFUNCTION_H_

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include "../Inversion/ObjectiveFunction.h"

namespace jif3D
  {
    //! This class provides a common denominator for the Matrix based regularization functions and other regularization functions such as MinimumSupport
    /*! The main purpose of this class is to provide a common base
     * for matrix based regularization and MinimumSupport based regularization. This
     * makes setting up the appropriate regularization easier as we can access
     * common functionality through a point to the base class.
     */
    class RegularizationFunction: public jif3D::ObjectiveFunction
      {
    private:
      //! The storage for a reference model
      jif3D::rvec Reference;
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & boost::serialization::base_object<ObjectiveFunction>(*this);
          ar & Reference;
        }
    public:
      //! The clone function provides a virtual constructor
      virtual RegularizationFunction *clone() const = 0;
      //! Set the reference model for the roughness calculation, this is optional
      void SetReferenceModel(const jif3D::rvec &Model)
        {
          Reference = Model;
        }
      //! Read only access to the reference model
      const jif3D::rvec &GetReferenceModel() const
        {
          return Reference;
        }
      RegularizationFunction()
        {
          // as this class only provides a common interface and minimal
          //storage, we declare constructor and destructor here
          //and omit the .cpp file
        }
      virtual ~RegularizationFunction()
        {

        }
      };

  } /* namespace jif3D */
#endif /* REGULARIZATIONFUNCTION_H_ */
