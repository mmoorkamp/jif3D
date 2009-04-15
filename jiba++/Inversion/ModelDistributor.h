//============================================================================
// Name        : ModelDistributor.h
// Author      : Apr 15, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef MODELDISTRIBUTOR_H_
#define MODELDISTRIBUTOR_H_

#include <cassert>
#include <vector>
#include <boost/function.hpp>
#include "../Global/VecMat.h"

namespace jiba
  {

    class ModelDistributor
      {
    public:
      typedef boost::function1<jiba::rvec, const jiba::rvec&> tModelTransform;
      ModelDistributor()
        {
        }
      virtual ~ModelDistributor()
        {
        }
    private:
      std::vector<tModelTransform> Transformers;
    public:
      jiba::rvec operator()(const jiba::rvec &FullModel, const size_t TransIndex)
        {
          assert(TransIndex < Transformers.size());
          return Transformers.at(TransIndex)(FullModel);
        }
      void AddTransformer(tModelTransform Trans)
        {
          Transformers.push_back(Trans);
        }
      };

    class CopyTransform
      {
    public:
      CopyTransform()
        {
        }
      virtual ~CopyTransform()
        {
        }
      jiba::rvec operator()(const jiba::rvec &FullModel)
        {
          return FullModel;
        }
      };
  }

#endif /* MODELDISTRIBUTOR_H_ */
