/*
 * SurfaceWavesModels.h
 *
 *  Created on: 19 Sep 2019
 *      Author: bweise
 */

#ifndef SURFACEWAVES_SURFACEWAVESMODEL_H_
#define SURFACEWAVES_SURFACEWAVESMODEL_H_

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
      virtual void WriteVTK(const std::string filename) const
        {
          ThreeDModelBase::WriteVTK(filename + ".vs", "Vs");
          Write3DModelToVTK(filename + ".vp", "Vp", GetXCoordinates(), GetYCoordinates(),
              GetZCoordinates(), GetVp());
          Write3DModelToVTK(filename + ".dens", "Density", GetXCoordinates(),
              GetYCoordinates(), GetZCoordinates(), GetVp());
        }
      virtual void WriteNetCDF(const std::string &filename);
      SurfaceWaveModel& operator=(const ThreeDModelBase& source);
      SurfaceWaveModel& operator=(const SurfaceWaveModel& source);

    private:
      t3DModelData DataVp, DataDens;
      };

    class SWExtraParameterSetter
      {
    public:
      void operator()(SurfaceWaveModel &Model, SurfaceWaveData &Data, const std::vector<double> &vel)
        {
          throw jif3D::FatalException(
              "Tomography class does not support extra parameters !", __FILE__, __LINE__);
        }
      };

  }

#endif /* SURFACEWAVES_SURFACEWAVESMODEL_H_ */
