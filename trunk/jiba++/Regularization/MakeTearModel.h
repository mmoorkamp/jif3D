//============================================================================
// Name        : MakeTearMode.h
// Author      : Jul 12, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================


#ifndef MAKETEARMODE_H_
#define MAKETEARMODE_H_

#include <cassert>
#include "../ModelBase/ThreeDModelBase.h"

namespace jiba
  {
    //! Create a dummy tear model with the same geometry as the geometry model, but with a value of 1.0 everywhere so tjat there is no tear
    inline void MakeTearModel(const jiba::ThreeDModelBase &Geometry,
        jiba::ThreeDModelBase &TearModel)
      {
        //create a model with a geometry that matches the inversion domain
        //we ignore the information about cell sizes etc.
        //in the end it only matters that the size and shape of the
        TearModel.SetData().resize(
            boost::extents[Geometry.GetData().shape()[0]][Geometry.GetData().shape()[1]][Geometry.GetData().shape()[2]]);
        //fill the object with 1 values, which means that there is no tear anywhere
        std::fill_n(TearModel.SetData().origin(),
            TearModel.GetData().num_elements(), 1.0);

      }
  }
#endif /* MAKETEARMODE_H_ */
