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

namespace jiba
  {

    class LimitedMemoryQuasiNewton: public jiba::NonLinearOptimization
      {
      private:
    	  double mu;
        const size_t MaxPairs;
        std::vector<boost::shared_ptr<jiba::rvec> > SHistory;
        std::vector<boost::shared_ptr<jiba::rvec> > YHistory;
        virtual void StepImplementation(jiba::rvec &CurrentModel) ;
    public:
      LimitedMemoryQuasiNewton(boost::shared_ptr<jiba::ObjectiveFunction> ObjFunction,const size_t n);
      virtual ~LimitedMemoryQuasiNewton();
      };

  }

#endif /* LIMITEDMEMORYQUASINEWTON_H_ */
