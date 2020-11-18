/*
 * ThreeDDCResistivityModel.h
 *
 *  Created on: Nov 13, 2020
 *      Author: zhanjie
 */

#ifndef THREEDDCRESISTIVITYMODEL_H_
#define THREEDDCRESISTIVITYMODEL_H_

#include "../Global/Serialization.h"
#include "../Global/FatalException.h"
#include "../Global/Jif3DGlobal.h"
#include "../ModelBase/ThreeDModelBase.h"

namespace jif3D
  {
    /** \addtogroup DC Resistivity classes and functions
     *
     * This module contains the functionality for DC resistivity forward modelling and data handling.
     * As for the other methods the 3D Models are stored in a class derived from ThreeDModelBase. The model parameter for these models is resistivity in ohm.m.
     * The DC resistivity model class stores all the model information necessary to calculate DCresistivity response, i.e. cell size and coordinate, resistivity value.
     * As this class is derived from ThreeDModelBase the overall handling is similar to the other model classes. There is a notable exceptions.
     * The model has to have variable cell sizes in z direction, which is thinner for shallow layers but thicker for deeper layers.
     *
     * For the source and receiver defination of measurement, we use the other class 'DCResistivityData' to set the source and receiver electrode parameters.
     * Moreover, there are often two sources and two receivers for each datum, so we store the Positive and Negative source at the same time. And the other receiver
     * electrode need to be stored, too. We store the source position only once, but store all receiver positions and source index for each datum.
     */
    /* @{ */
    class DCExtraParameterSetter;

    class DCResistivityData;

    class J3DEXPORT ThreeDDCResistivityModel: public jif3D::ThreeDModelBase
      {
    public:
      typedef DCExtraParameterSetter ExtraParameterSetter;

    public:
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & base_object<ThreeDModelBase>(*this);
        }
      //! The DC Resistivity model for DCResistivityCalculator has the same cell size for all cells in each horizontal directions so we just have one function
      //! to set it
      /*! This function sets both the size of all cells as well as the number of cells in the horizontal (x and y) directions
       * @param XSize The size of each cell in x-direction (North) in m
       * @param YSize The size of each cell in y-direction (East) in m
       * @param nx The number of cells in x-direction (North)
       * @param ny The number of cells in y-direction (East)
       */
      void SetHorizontalCellSize(const double XSize, const double YSize, const size_t nx,
          const size_t ny)
        {
          ThreeDModelBase::t3DModelDim XCS(nx,XSize);
          ThreeDModelBase::t3DModelDim YCS(ny,YSize);
          ThreeDModelBase::SetXCellSizes(XCS);
          ThreeDModelBase::SetYCellSizes(YCS);
        }
      //! The depth direction, cells can all have different sizes so we allow direct access to the CellSize structure
      void SetZCellSizes(ThreeDModelBase::t3DModelDim &ZCS)
        {
    	  ThreeDModelBase::SetZCellSizes(ZCS);
        }

      //! The seismic model has the same cell size for all cells in all directions so we just have one function to set it
      /*! This function sets both the size of all cells as well as the number of cells in x,y and z-direction
       * @param Size The size of each cell in all directions in m
       * @param nx The number of cells in x-direction (North)
       * @param ny The number of cells in y-direction (East)
       * @param nz The number of cells in z-direction
       */
      void SetCellSize(const double Size, const size_t nx, const size_t ny,
          const size_t nz)
        {
          ThreeDModelBase::t3DModelDim XCS(nx,Size);
          ThreeDModelBase::t3DModelDim YCS(ny,Size);
          ThreeDModelBase::t3DModelDim ZCS(nz,Size);
          SetXCellSizes(XCS);
          SetYCellSizes(YCS);
          SetZCellSizes(ZCS);
          ThreeDModelBase::SetData().resize(boost::extents[nx][ny][nz]);
        }


      //! return read only access to the stored resistivity values
      const t3DModelData &GetResistivities() const
        {
          return ThreeDModelBase::GetData();
        }
      //! return a reference to stored resistivities
      t3DModelData &SetResistivities()
        {
          return ThreeDModelBase::SetData();
        }

      //! Write the resistivity model in VTK format, at the moment the best format for plotting, this only writes the resistivity and not the source and receiver positions
      void WriteVTK(const std::string filename) const
        {
          ThreeDModelBase::WriteVTK(filename, "Resistivity");
        }
      //! Write the resistivity model and all associated information in a netcdf file
      virtual void WriteNetCDF(const std::string &filename) const override;
      //! Read the resistivity model and all associated information from a netcdf file
      virtual void ReadNetCDF(const std::string &filename) override;
      ThreeDDCResistivityModel();
      ThreeDDCResistivityModel& operator=(const ThreeDModelBase& source);
      virtual ~ThreeDDCResistivityModel();
      };

    class DCExtraParameterSetter
      {
    public:
      void operator()(ThreeDDCResistivityModel &Model, DCResistivityData &Data, const std::vector<double> &Dens)
        {
          throw jif3D::FatalException(
              "Tomography class does not support extra parameters !", __FILE__, __LINE__);
        }
      };
  /* @} */
  }

#endif /* THREEDDCRESISTIVITYMODEL_H_ */
