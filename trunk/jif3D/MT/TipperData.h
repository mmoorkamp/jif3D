/*
 * TipperData.h
 *
 *  Created on: May 22, 2019
 *      Author: max
 */

#ifndef MT_TIPPERDATA_H_
#define MT_TIPPERDATA_H_

#include "../Global/Serialization.h"
#include "EMData.h"
#include "X3DModel.h"

namespace jif3D
  {

    class TipperData: public EMData
      {
    private:
      std::vector<int> HxIndices;
      std::vector<int> HyIndices;
      std::vector<int> HzIndices;
    public:
      typedef X3DModel ModelType;
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & base_object<EMData>(*this);
          ar & HxIndices;
          ar & HyIndices;
          ar & HzIndices;
        }

      const std::vector<int> &GetHxIndices() const
        {
          return HxIndices;
        }
      const std::vector<int> &GetHyIndices() const
        {
          return HyIndices;
        }
      const std::vector<int> &GetHzIndices() const
        {
          return HzIndices;
        }
      //! Set the indices for magnetic field measurements
      void SetFieldIndices(std::vector<int> &value_Hx, std::vector<int> &value_Hy,
          const std::vector<int> &value_Hz)
        {

          if (HxIndices.size() != HyIndices.size())
            {
              throw jif3D::FatalException(
                  "Number of Hy magnetic field indices is not the same as the Hx field indices!",
                  __FILE__, __LINE__);
            }

          if (HxIndices.size() != HzIndices.size())
            {
              throw jif3D::FatalException(
                  "Number of Hz magnetic field indices is not the same as the Hx field indices!",
                  __FILE__, __LINE__);
            }

          HxIndices = value_Hx;
          HyIndices = value_Hy;
          HzIndices = value_Hz;
        }
      void CompleteObject();
      //! Write all information to a netcdf file
      virtual void ReadNetCDF(const std::string &filename) override;
      virtual void WriteNetCDF(const std::string &filename) override;
      void ReadModEM(const std::string &filename);
      void WriteModEM(const std::string &filename);
      void PlotMeasurementConfiguration(const std::string &filename);
      TipperData();
      virtual ~TipperData();
      };

  } /* namespace jif3D */

#endif /* MT_TIPPERDATA_H_ */
