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
      const t3DModelData &GetVp() const
        {
          return DataVp;
        }
      const t3DModelData &GetDens() const
        {
          return DataDens;
        }
      void SetCellCoords(const std::vector<double> &nx, const std::vector<double> &ny,
          const std::vector<double> &nz)
        {
          SetXCoordinates(nx);
          SetYCoordinates(ny);
          SetZCoordinates(nz);
          ThreeDModelBase::SetData().resize(
              boost::extents[nx.size() - 1][ny.size() - 1][nz.size() - 1]);
          DataVp.resize(boost::extents[nx.size() - 1][ny.size() - 1][nz.size() - 1]);
          DataDens.resize(boost::extents[nx.size() - 1][ny.size() - 1][nz.size() - 1]);
          BGDens.resize(boost::extents[nx.size() - 1][ny.size() - 1][nz.size() - 1]);
        }
      void SetVp(const t3DModelData &vp)
        {
          DataVp.resize(boost::extents[vp.shape()[0]][vp.shape()[1]][vp.shape()[2]]);
          DataVp = vp;

        }
      void SetDens(const t3DModelData& dens)
        {
          DataDens.resize(
              boost::extents[dens.shape()[0]][dens.shape()[1]][dens.shape()[2]]);
          DataDens = dens;
          BGDens.resize(
              boost::extents[dens.shape()[0]][dens.shape()[1]][dens.shape()[2]]);
          BGDens = dens;
        }
      void SetDensAnomaly(const t3DModelData& densan)
        {
          if (densan.num_elements() != BGDens.num_elements())
            throw jif3D::FatalException(
                "Not the right amount of elements to set background density", __FILE__,
                __LINE__);
          std::transform(densan.origin(), densan.origin() + densan.num_elements(),
              BGDens.origin(), DataDens.origin(), std::plus<double>());
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
      t3DModelData DataVp, DataDens, BGDens;
      };

    class SWExtraParameterSetter
      {
    public:
      void operator()(SurfaceWaveModel &Model, SurfaceWaveData &Data,
          const std::vector<double> &denspvel)
        {
          const size_t ncells = Model.GetNModelElements();
          if (denspvel.size() != 2 * ncells)
            {
              throw jif3D::FatalException(
                  "Number of P-velocities and densities does not match grid size",
                  __FILE__,
                  __LINE__);
            }

          SurfaceWaveModel::t3DModelData values(Model.GetVp());
          std::copy(denspvel.begin(), denspvel.begin() + ncells, values.origin());
          Model.SetVp(values);
          std::copy(denspvel.begin() + ncells, denspvel.begin() + 2 * ncells,
              values.origin());
          Model.SetDensAnomaly(values);

        }
      };

  }

#endif /* SURFACEWAVES_SURFACEWAVESMODEL_H_ */
