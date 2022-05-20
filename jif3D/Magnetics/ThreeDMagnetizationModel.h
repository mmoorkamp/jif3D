/*
 * ThreeDMagnetizationModel.h
 *
 *  Created on: Apr 26, 2022
 *      Author: max
 */

#ifndef MAGNETICS_THREEDMAGNETIZATIONMODEL_H_
#define MAGNETICS_THREEDMAGNETIZATIONMODEL_H_

#include "../Global/Serialization.h"
#include "../Global/Jif3DGlobal.h"
#include "../Global/VecMat.h"
#include "../Global/FatalException.h"
#include "../Global/NetCDFPortHelper.h"

#include "../ModelBase/ThreeDModelBase.h"
#include "../ModelBase/GridOnlyModelCache.h"

namespace jif3D
  {

    class ThreeDMagnetizationModel: public ThreeDModelBase
      {
    public:
      //! The type of the background thickness and susceptibility vector, this is a std::vector because we want to easily append elements
      typedef std::vector<double> tBackgroundVec;
      typedef jif3D::GridOnlyModelCache<ThreeDMagnetizationModel> ModelCacheType;
    private:
      //! The susceptibilities of the background layers
      tBackgroundVec bg_susceptibilities;
      //! The thicknesses of the background layers
      tBackgroundVec bg_thicknesses;
      ThreeDModelBase::t3DModelData Magnetization_Y;
      ThreeDModelBase::t3DModelData Magnetization_Z;
    public:
      ThreeDMagnetizationModel();
      virtual ~ThreeDMagnetizationModel();
      jif3D::rvec GetModelParameters() const
        {
          const size_t ngrid = GetData().num_elements();
          jif3D::rvec parms(3 * ngrid);
          std::copy(GetData().origin(), GetData().origin() + ngrid,
              parms.begin());
          std::copy(Magnetization_Y.origin(), Magnetization_Y.origin() + ngrid,
              parms.begin() + ngrid);
          std::copy(Magnetization_Z.origin(), Magnetization_Z.origin() + ngrid,
                        parms.begin() + 2 * ngrid);
          return parms;
        }
      //! The implementation class for magnetic and gravity data needs to know the number of model parameters
      size_t GetNModelParm() const
        {
          return 3 * GetData().num_elements();
        }
      void SetCellCoords(const std::vector<double> &nx, const std::vector<double> &ny,
          const std::vector<double> &nz)
        {
          SetXCoordinates(nx);
          SetYCoordinates(ny);
          SetZCoordinates(nz);
          ThreeDModelBase::SetData().resize(
              boost::extents[nx.size() - 1][ny.size() - 1][nz.size() - 1]);
          Magnetization_Y.resize(
              boost::extents[nx.size() - 1][ny.size() - 1][nz.size() - 1]);
          Magnetization_Y.resize(
              boost::extents[nx.size() - 1][ny.size() - 1][nz.size() - 1]);
        }
      virtual void SetMeshSize(const size_t nx, const size_t ny, const size_t nz) override
      {
          ThreeDModelBase::SetMeshSize(nx, ny, nz);
          Magnetization_Y.resize(boost::extents[nx][ny][nz]);
          Magnetization_Z.resize(boost::extents[nx][ny][nz]);
      }
      t3DModelData& SetMagnetization_X()
        {
          return SetData();
        }
      t3DModelData& SetMagnetization_Y()
        {
          return Magnetization_Y;
        }
      t3DModelData& SetMagnetization_Z()
        {
          return Magnetization_Z;
        }
      const t3DModelData& GetMagnetization_X() const
        {
          return GetData();
        }
      const t3DModelData& GetMagnetization_Y() const
        {
          return Magnetization_Y;
        }
      const t3DModelData& GetMagnetization_Z() const
        {
          return Magnetization_Z;
        }
      //! Set the sizes of the grid cells in x-direction in m
      void SetXCellSizes(ThreeDModelBase::t3DModelDim &XCS)
        {
          ThreeDModelBase::SetXCellSizes(XCS);
        }
      //! Set the sizes of the grid cells in y-direction in m
      void SetYCellSizes(ThreeDModelBase::t3DModelDim &YCS)
        {
          ThreeDModelBase::SetYCellSizes(YCS);
        }
      //! Set the sizes of the grid cells in z-direction in m
      void SetZCellSizes(ThreeDModelBase::t3DModelDim &ZCS)
        {
          ThreeDModelBase::SetZCellSizes(ZCS);
        }
      virtual void WriteVTK(const std::string filename) const;
      virtual void WriteNetCDF(const std::string &filename) const override;
      virtual void ReadNetCDF(const std::string &model_file) override;

      };

  } /* namespace jif3D */

#endif /* MAGNETICS_THREEDMAGNETIZATIONMODEL_H_ */
