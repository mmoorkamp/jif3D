//============================================================================
// Name        : CurvatureRegularization.h
// Author      : Jan 11, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================


#ifndef CURVATUREREGULARIZATION_H_
#define CURVATUREREGULARIZATION_H_

#include "MatOpRegularization.h"

namespace jiba
  {

    class CurvatureRegularization: public jiba::MatOpRegularization
      {
    private:
      const double Eps;
      void ConstructOperator(const jiba::ThreeDModelBase &ModelGeometry);
    public:
      CurvatureRegularization(const jiba::ThreeDModelBase &Geometry,
          const double ModEps = 1e-8) :
        MatOpRegularization(Geometry), Eps(ModEps)
        {
          ConstructOperator(Geometry);
        }
      virtual ~CurvatureRegularization();
      };

  }

#endif /* CURVATUREREGULARIZATION_H_ */
