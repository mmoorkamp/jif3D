//============================================================================
// Name        : EqualGeometry.h
// Author      : Mar 3, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================


#ifndef EQUALGEOMETRY_H_
#define EQUALGEOMETRY_H_

#include "ThreeDModelBase.h"

namespace jiba
  {
    //! Test if two models have the same grid geometry
    /*! This function can be used to test whether two 3D models have the same grid geometry, i.e. the same number of cells
     * in each direction with the same sizes.
     * @param Model1 The first model for the comparison
     * @param Model2 the second model for the comparison
     * @return True if the geometries are equal, false otherwise
     */
    inline bool EqualGridGeometry(const jiba::ThreeDModelBase &Model1,
        const jiba::ThreeDModelBase &Model2)
      {
        //we first make a cheap test and then the more expensive ones
        if (Model1.GetNModelElements() != Model2.GetNModelElements())
          return false;
        // as soon as one fails we return false
        //we check the sizes of the cells in all three directions
        if (!std::equal(Model1.GetXCellSizes().origin(),
            Model1.GetXCellSizes().origin()
                + Model1.GetXCellSizes().num_elements(),
            Model2.GetXCellSizes().origin()))
          return false;
        if (!std::equal(Model1.GetYCellSizes().origin(),
            Model1.GetYCellSizes().origin()
                + Model1.GetYCellSizes().num_elements(),
            Model2.GetYCellSizes().origin()))
          return false;
        if (!std::equal(Model1.GetZCellSizes().origin(),
            Model1.GetZCellSizes().origin()
                + Model1.GetZCellSizes().num_elements(),
            Model2.GetZCellSizes().origin()))
          return false;
        //we only get here if everything is equal
        return true;
      }
  }
#endif /* EQUALGEOMETRY_H_ */
