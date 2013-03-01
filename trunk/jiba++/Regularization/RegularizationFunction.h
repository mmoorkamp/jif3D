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
#include "ObjectiveFunction.h"

namespace jiba
  {

    class RegularizationFunction: public jiba::ObjectiveFunction
      {
    private:
      //! The storage for a reference model
      jiba::rvec Reference;
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
      void SetReferenceModel(const jiba::rvec &Model)
        {
          Reference = Model;
        }
      const jiba::rvec &GetReferenceModel() const
        {
          return Reference;
        }
      RegularizationFunction();
      virtual ~RegularizationFunction();
      };

  } /* namespace jiba */
#endif /* REGULARIZATIONFUNCTION_H_ */
