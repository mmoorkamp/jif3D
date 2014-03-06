//============================================================================
// Name        : ThreeDDCResistivityModel.h
// Author      : 6 Mar 2014
// Version     : 
// Copyright   : 2014, mm489
//============================================================================

#ifndef THREEDDCRESISTIVITYMODEL_H_
#define THREEDDCRESISTIVITYMODEL_H_

#include "../ModelBase/ThreeDModelBase.h"

namespace jif3D
  {
  /** \addtogroup dcresistivity DC Resistivity classes and functions */
  /* @{ */
    class ThreeDDCResistivityModel: public ThreeDModelBase
      {
    public:
      ThreeDDCResistivityModel();
      virtual ~ThreeDDCResistivityModel();
      };
    /* @} */
  } /* namespace jif3D */

#endif /* THREEDDCRESISTIVITYMODEL_H_ */
