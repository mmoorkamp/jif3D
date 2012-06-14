//============================================================================
// Name        : MinDiffRegularization.h
// Author      : May 7, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#ifndef MINDIFFREGULARIZATION_H_
#define MINDIFFREGULARIZATION_H_

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include "MatOpRegularization.h"

namespace jiba
  {
    /** \addtogroup Regularization classes to regularize the inversion */
    /* @{ */
    //! Regularize by calculating the difference to a reference model
    /*! This class can be used to regularize the inversion by
     * looking for a model that is as close as possible to a reference model. If no reference
     * model is specified it is assumed to be zero and therefore the smallest model
     * is sought.
     */
    class MinDiffRegularization: public MatOpRegularization
      {
    private:
      void ConstructOperator(const jiba::ThreeDModelBase &ModelGeometry)
        {
          const size_t ngrid = ModelGeometry.GetData().num_elements();

          for (size_t i = 0; i < ngrid; ++i)
            {
              XOperatorMatrix(i, i) = 1.0;
            }
        }
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & boost::serialization::base_object<MatOpRegularization>(*this);
        }
    public:
      //! The clone function provides a virtual constructor
      virtual MinDiffRegularization *clone() const
        {
          return new MinDiffRegularization(*this);
        }

      MinDiffRegularization(const jiba::ThreeDModelBase &Geometry) :
          MatOpRegularization(Geometry)
        {
          ConstructOperator(Geometry);
        }
      virtual ~MinDiffRegularization()
        {
        }
      }
    ;
  /* @} */
  }

#endif /* MINDIFFREGULARIZATION_H_ */
