/*
 * SurfaceWavesModels.h
 *
 *  Created on: 19 Sep 2019
 *      Author: bweise
 */

#ifndef SURFACEWAVES_SURFACEWAVEMODEL_H_
#define SURFACEWAVES_SURFACEWAVEMODEL_H_

#include "../ModelBase/ThreeDModelBase.h"
#include "../ModelBase/VTKTools.h"
#include "SurfaceWaveData.h"

namespace jif3D
  {
    class SWExtraParameterSetter;

    class J3DEXPORT SurfaceWaveModel: public jif3D::ThreeDModelBase
      {
    public:
      typedef SWExtraParameterSetter ExtraParameterSetter;
      SurfaceWaveModel();
      virtual void ReadNetCDF(const std::string &model_file) override;
      t3DModelData GetVp() const
        {
          return DataVp;
        }
      t3DModelData GetDens() const
        {
          return DataDens;
        }
      void SetCellCoords(const std::vector<double> &nx, const std::vector<double> &ny,
          const std::vector<double> &nz)
        {
          ThreeDModelBase::t3DModelDim XCS(nx.size());
          ThreeDModelBase::t3DModelDim YCS(ny.size());
          ThreeDModelBase::t3DModelDim ZCS(nz.size());
          std::copy(nx.begin(), nx.end(), XCS.begin());
          std::copy(ny.begin(), ny.end(), YCS.begin());
          std::copy(nz.begin(), nz.end(), ZCS.begin());
          SetXCoordinates(XCS);
          SetYCoordinates(YCS);
          SetZCoordinates(ZCS);
          ThreeDModelBase::SetData().resize(
              boost::extents[nx.size()-1][ny.size()-1][nz.size()-1]);
          DataVp.resize(
              boost::extents[nx.size()-1][ny.size()-1][nz.size()-1]);
          DataDens.resize(
              boost::extents[nx.size()-1][ny.size()-1][nz.size()-1]);

        }
      t3DModelData& SetVp()
        {
          return DataVp;
        }
      t3DModelData& SetDens()
        {
          return DataDens;
        }
      virtual void WriteVTK(const std::string filename) const
        {
          ThreeDModelBase::WriteVTK(filename + ".vs.vtk", "Vs");
          Write3DModelToVTK(filename + ".vp.vtk", "Vp", GetXCoordinates(),
              GetYCoordinates(), GetZCoordinates(), GetVp());
          Write3DModelToVTK(filename + ".dens.vtk", "Density", GetXCoordinates(),
              GetYCoordinates(), GetZCoordinates(), GetDens());
        }
      virtual void WriteNetCDF(const std::string &filename) const override;
      SurfaceWaveModel& operator=(const ThreeDModelBase &source);
      SurfaceWaveModel& operator=(const SurfaceWaveModel &source);

    private:
      t3DModelData DataVp, DataDens;
      };

    class SWExtraParameterSetter
      {
    public:
      void operator()(SurfaceWaveModel &Model, SurfaceWaveData &Data,
          const std::vector<double> &vel)
        {
          throw jif3D::FatalException(
              "Tomography class does not support extra parameters !", __FILE__, __LINE__);
        }
      };

  }

#endif /* SURFACEWAVES_SURFACEWAVESMODEL_H_ */
