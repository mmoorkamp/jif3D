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
      void ConstructOperator(const jiba::ThreeDModelBase &ModelGeometry);
    public:
      CurvatureRegularization(const jiba::ThreeDModelBase &Geometry) :
        MatOpRegularization(Geometry)
        {
          ConstructOperator(Geometry);
        }
      virtual ~CurvatureRegularization();
      };

  }

#endif /* CURVATUREREGULARIZATION_H_ */
