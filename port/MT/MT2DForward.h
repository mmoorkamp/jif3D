//============================================================================
// Name        : MT2DForward.h
// Author      : Sep 24, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#ifndef MT2DFORWARD_H_
#define MT2DFORWARD_H_

#include <boost/multi_array.hpp>

#include "../Global/Jif3DGlobal.h"

namespace jif3D
  {
    /** \addtogroup mtmodelling Forward modelling of magnetotelluric data */
    /* @{ */

    //! Calculate the magnetoteluric response for a 2D model
    /*! This class encapsulates calls to the 2D MT Fortran code
     * by P. Tarits. It only contains some very basic features
     * as it is not actively used at the moment.
     */
    class J3DEXPORT MT2DForward
      {
    public:
      typedef boost::multi_array<double, 1> t2DModelDim;
      typedef boost::multi_array<double, 2> t2DModelData;
      typedef boost::multi_array<double, 3> t2DModelFields;
    private:
      t2DModelDim XSizes;
      t2DModelDim ZSizes;
      t2DModelData Resistivities;
      t2DModelFields Hx_real;
      t2DModelFields Hx_imag;
      t2DModelFields Hy_real;
      t2DModelFields Hy_imag;
      t2DModelFields Ey_real;
      t2DModelFields Ey_imag;
      t2DModelFields Ex_real;
      t2DModelFields Ex_imag;
      t2DModelFields Ez_real;
      t2DModelFields Ez_imag;
    public:
      t2DModelDim GetXSizes() const
        {
          return XSizes;
        }

      t2DModelDim &SetXSizes()
        {
          return this->XSizes;
        }

      t2DModelDim GetZSizes() const
        {
          return ZSizes;
        }

      t2DModelDim &SetZSizes()
        {
          return this->ZSizes;
        }

      t2DModelData GetResistivities() const
        {
          return Resistivities;
        }

      t2DModelData &SetResistivities()
        {
          return this->Resistivities;
        }

      t2DModelFields GetHx_real() const
        {
          return Hx_real;
        }

      t2DModelFields GetHx_imag() const
        {
          return Hx_imag;
        }

      t2DModelFields GetHy_real() const
        {
          return Hy_real;
        }

      t2DModelFields GetHy_imag() const
        {
          return Hy_imag;
        }

      t2DModelFields GetEy_real() const
        {
          return Ey_real;
        }

      t2DModelFields GetEy_imag() const
        {
          return Ey_imag;
        }

      t2DModelFields GetEx_real() const
        {
          return Ex_real;
        }

      t2DModelFields GetEx_imag() const
        {
          return Ex_imag;
        }
      //! Calculate the response for the E-polarization
      void CalcEpol(const std::vector<double> &Periods);
      //! Calculate the response of the B-polarization
      void CalcBpol(const std::vector<double> &Periods);
      MT2DForward();
      virtual ~MT2DForward();
      };
  /* @} */
  }

#endif /* MT2DFORWARD_H_ */
