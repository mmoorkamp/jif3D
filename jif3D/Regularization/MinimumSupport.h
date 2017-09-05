//============================================================================
// Name        : MinimumSupport.h
// Author      : Jul 9, 2012
// Version     :
// Copyright   : 2012, mmoorkamp
//============================================================================

#ifndef MINIMUMSUPPORT_H_
#define MINIMUMSUPPORT_H_


#include "../Global/Serialization.h"
#include "../Global/Jif3DGlobal.h"
#include "RegularizationFunction.h"
#include "MatOpRegularization.h"
#include <boost/shared_ptr.hpp>

namespace jif3D
  {

    class J3DEXPORT MinimumSupport: public jif3D::RegularizationFunction
      {
    private:
      size_t ModelNumber;
      jif3D::rvec RegDiff;
      double beta;
      boost::shared_ptr<jif3D::MatOpRegularization> RegFunc;
      jif3D::ThreeDModelBase Geometry;
      friend class access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & base_object<RegularizationFunction>(*this);
          ar & ModelNumber;
          ar & RegDiff;
          ar & beta;
          ar & RegFunc;
          ar & Geometry;
        }
    public:
      //! The clone function provides a virtual constructor
      virtual MinimumSupport *clone() const override
        {
          return new MinimumSupport(*this);
        }
      virtual void ImplDataDifference(const jif3D::rvec &Model, jif3D::rvec &Diff) override;
      virtual jif3D::rvec ImplGradient(const jif3D::rvec &Model, const jif3D::rvec &Diff) override;
      MinimumSupport(boost::shared_ptr<jif3D::MatOpRegularization> RF, double b = 1.0);
      virtual ~MinimumSupport();
      };

  } /* namespace jif3D */
#endif /* MINIMUMSUPPORT_H_ */
