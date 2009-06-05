//============================================================================
// Name        : LimitedMemoryQuasiNewton.h
// Author      : Apr 20, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef LIMITEDMEMORYQUASINEWTON_H_
#define LIMITEDMEMORYQUASINEWTON_H_

#include "NonLinearOptimization.h"
#include <boost/shared_ptr.hpp>
#include <vector>
#include "../Gravity/ThreeDGravityModel.h"

namespace jiba
  {

    class LimitedMemoryQuasiNewton: public jiba::NonLinearOptimization
      {
      private:
    	  jiba::ThreeDGravityModel SearchModel;
    	double mu;
        const size_t MaxPairs;
        std::vector<boost::shared_ptr<jiba::rvec> > SHistory;
        std::vector<boost::shared_ptr<jiba::rvec> > YHistory;
        virtual void StepImplementation(jiba::rvec &CurrentModel) ;
    public:
        void SetModelGeometry(const jiba::ThreeDGravityModel &Model)
          {
            SearchModel = Model;
          }
      LimitedMemoryQuasiNewton(boost::shared_ptr<jiba::ObjectiveFunction> ObjFunction,const size_t n);
      virtual ~LimitedMemoryQuasiNewton();
      };

  }

#endif /* LIMITEDMEMORYQUASINEWTON_H_ */
