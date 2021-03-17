/*
 * MTData.h
 *
 *  Created on: May 21, 2019
 *      Author: max
 */

#ifndef MT_MTDATA_H_
#define MT_MTDATA_H_

#include "../Global/Serialization.h"
#include "EMData.h"
#include "X3DModel.h"
#include <vector>

namespace jif3D
  {

    class MTData: public EMData
      {
    private:
      std::vector<double> Distortion;
      std::vector<int> ExIndices;
      std::vector<int> EyIndices;
      std::vector<int> HxIndices;
      std::vector<int> HyIndices;
    public:
      typedef X3DModel ModelType;
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & base_object<EMData>(*this);
          ar & Distortion;
          ar & ExIndices;
          ar & EyIndices;
          ar & HxIndices;
          ar & HyIndices;
        }
      const std::vector<double> &GetDistortion() const
        {
          return Distortion;
        }
      void SetDistortion(const std::vector<double> &D)
        {
          Distortion = D;
        }
      void SetDefaultDistortion();
      const std::vector<int> &GetExIndices() const
        {
          return ExIndices;

        }
      const std::vector<int> &GetEyIndices() const
        {
          return EyIndices;
        }
      const std::vector<int> &GetHxIndices() const
        {
          return HxIndices;
        }
      const std::vector<int> &GetHyIndices() const
        {
          return HyIndices;
        }
      void CompleteObject();
      //! Set the indices for electric and magnetic field measurements
      void SetFieldIndices(const std::vector<int> &value_Ex,
          const std::vector<int> &value_Ey, std::vector<int> &value_Hx,
          std::vector<int> &value_Hy)
        {

          if (HxIndices.size() != ExIndices.size())
            {
              throw jif3D::FatalException(
                  "Number of Hx magnetic field indices is not the same as the electric field indices!",
                  __FILE__, __LINE__);
            }

          if (HyIndices.size() != ExIndices.size())
            {
              throw jif3D::FatalException(
                  "Number of Hy magnetic field indices is not the same as the electric field indices!",
                  __FILE__, __LINE__);
            }

          if (EyIndices.size() != ExIndices.size())
            {
              throw jif3D::FatalException(
                  "Number of Ey field indices is not the same as the Ex field indices!",
                  __FILE__, __LINE__);
            }

          ExIndices = value_Ex;
          EyIndices = value_Ey;
          HxIndices = value_Hx;
          HyIndices = value_Hy;
        }
      MTData();
      virtual ~MTData();
      //! Write all information to a netcdf file
      virtual void ReadNetCDF(const std::string &filename) override;
      virtual void WriteNetCDF(const std::string &filename) const override;
      void ReadModEM(const std::string &filename);
      void WriteModEM(const std::string &filename) const;
      void PlotMeasurementConfiguration(const std::string &filename) const;
      };

  } /* namespace jif3D */

#endif /* MT_MTDATA_H_ */
