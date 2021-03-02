//============================================================================
// Name        : EqualGeometry.h
// Author      : Mar 3, 2010
// Version     :
// Copyright   : 2010, mmoorkamp
//============================================================================

#ifndef EQUALGEOMETRY_H_
#define EQUALGEOMETRY_H_

#include <iostream>
#include "../Global/Jif3DGlobal.h"
#include "../Global/NumUtil.h"
#include "../Global/FatalException.h"
#include "ThreeDModelBase.h"

namespace jif3D
  {
    /** \addtogroup modelbase Basic classes and routines for 3D models */
    /* @{ */
    //! Test if two models have the same grid geometry
    /*! This function can be used to test whether two 3D models have the same grid geometry, i.e. the same number of cells
     * in each direction with the same sizes. It also checks if the origins match
     * @param Model1 The first model for the comparison
     * @param Model2 the second model for the comparison
     * @return True if the geometries are equal, false otherwise
     */
    inline bool EqualGridGeometry(const jif3D::ThreeDModelBase &Model1,
        const jif3D::ThreeDModelBase &Model2)
      {
        //we first make a cheap test and then the more expensive ones
        if (Model1.GetNModelElements() != Model2.GetNModelElements())
          {
            std::cerr << "Models have different sizes: " << Model1.GetNModelElements()
                << " " << Model2.GetNModelElements() << std::endl;
            return false;
          }
        //compare the origin of the models to a precision of 1cm
        const double tolerance = 0.01;
        if (!roughlyEqual<double, double>(tolerance)(Model1.GetXOrigin(),
            Model2.GetXOrigin()))
          {
            std::cerr << "Models have different x-origins: " << Model1.GetXOrigin() << " "
                << Model2.GetXOrigin() << std::endl;
            return false;
          }

        if (!roughlyEqual<double, double>(tolerance)(Model1.GetYOrigin(),
            Model2.GetYOrigin()))
          {
            std::cerr << "Models have different y-origins: " << Model1.GetYOrigin() << " "
                << Model2.GetYOrigin() << std::endl;
            return false;
          }

        if (!roughlyEqual<double, double>(tolerance)(Model1.GetZOrigin(),
            Model2.GetZOrigin()))
          {
            std::cerr << "Models have different z-origins: " << Model1.GetZOrigin() << " "
                << Model2.GetZOrigin() << std::endl;
            return false;
          }
        // as soon as one fails we return false
        //we check the sizes of the cells in all three directions
        const std::size_t nx = Model1.GetXCellSizes().size();
        const std::size_t ny = Model1.GetYCellSizes().size();
        const std::size_t nz = Model1.GetZCellSizes().size();
        //check if we have the same number of cells in each direction
        if (nx != Model2.GetXCellSizes().size())
          {
            std::cerr << "Different number of cells in x-direction " << nx << " "
                << Model2.GetXCellSizes().size() << std::endl;
            return false;
          }

        if (ny != Model2.GetYCellSizes().size())
          {
            std::cerr << "Different number of cells in y-direction " << ny << " "
                << Model2.GetYCellSizes().size() << std::endl;
            return false;
          }

        if (nz != Model2.GetZCellSizes().size())
          {
            std::cerr << "Different number of cells in z-direction " << nz << " "
                << Model2.GetZCellSizes().size() << std::endl;
            return false;
          }
        //now check all the cell sizes, given the issues with floating point
        //comparisons we cannot simply do a straight comparison, so we do an approximate one
        //for our models a match within 1cm should be more than sufficient

        for (std::size_t i = 0; i < nx; ++i)
          {
            if (!roughlyEqual<double, double>(tolerance)(Model1.GetXCellSizes()[i],
                Model2.GetXCellSizes()[i]))
              {
                std::cerr << "Cell sizes in x-direction do not match "
                    << Model1.GetXCellSizes()[i] << " " << Model2.GetXCellSizes()[i]
                    << std::endl;
                return false;
              }
          }
        for (std::size_t i = 0; i < ny; ++i)
          {
            if (!roughlyEqual<double, double>(tolerance)(Model1.GetYCellSizes()[i],
                Model2.GetYCellSizes()[i]))
              {
                std::cerr << "Cell sizes in y-direction do not match "
                    << Model1.GetYCellSizes()[i] << " " << Model2.GetYCellSizes()[i]
                    << std::endl;
                return false;
              }
          }
        for (std::size_t i = 0; i < nz; ++i)
          {
            if (!roughlyEqual<double, double>(tolerance)(Model1.GetZCellSizes()[i],
                Model2.GetZCellSizes()[i]))
              {
                std::cerr << "Cell sizes in z-direction do not match "
                    << Model1.GetZCellSizes()[i] << " " << Model2.GetZCellSizes()[i]
                    << std::endl;
                return false;
              }
          }
        //we only get here if everything is equal
        return true;
      }
  /* @} */
  }
#endif /* EQUALGEOMETRY_H_ */
