//============================================================================
// Name        : ThreeDMTModel.h
// Author      : Apr 8, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#ifndef THREEDMTMODEL_H_
#define THREEDMTMODEL_H_

#include "../Global/Serialization.h"
#include "../Global/Jif3DGlobal.h"
#include "../ModelBase/ThreeDModelBase.h"

namespace jif3D
  {
    /** \addtogroup mtmodelling Forward modelling of magnetotelluric data */
    /* @{ */

    //we use these names for storing the model in NetCDF and vtk files
    static const std::string ConductivityName = "Conductivity";
    static const std::string ConductivityUnit = "S/m";

    //! A helper class for the template ThreeDModelObjective that lets us set distortion values as extra inversion parameters
    class MTDistortionSetter;
    //! This class stores all information associated with 3D magnetotelluric models
    /*! This class extends ThreeDModelBase to store the calculation frequencies
     * which are required for any MT forward calculation. It also provides named
     * access to the model through SetConductivities and writting file for vtk.
     * All other functionality and properties that are specific to a certain
     * forward code have to be implemented in a derived class.
     */
    class J3DEXPORT ThreeDMTModel : public jif3D::ThreeDModelBase
      {
    private:
      //! The four real elements of the distortion matrix C for each Station
      std::vector<double> DistortionParameters;
      //! The calculation frequencies in Hz
      std::vector<double> Frequencies;
    public:
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & base_object<ThreeDModelBase>(*this);
          ar & DistortionParameters;
          ar & Frequencies;
        }
      typedef MTDistortionSetter ExtraParameterSetter;
      const std::vector<double> &GetDistortionParameters() const
        {
          return DistortionParameters;
        }
      void SetDistortionParameters(const std::vector<double> &Dist)
        {
          DistortionParameters = Dist;
        }
      //! Get the vector of calculation frequencies in Hz, read only
      const std::vector<double> &GetFrequencies() const
        {
          return Frequencies;
        }
      //! Set the vector of calculation frequencies in Hz
      std::vector<double> &SetFrequencies()
        {
          return Frequencies;
        }
      //! return read only access to the stored conductivity values
      const t3DModelData &GetConductivities() const
        {
          return ThreeDModelBase::GetData();
        }
      //! return a reference to stored conductivities
      t3DModelData &SetConductivities()
        {
          return ThreeDModelBase::SetData();
        }
      //! Copy the source and receiver positions and the indices for the source-receiver combinations from the source model, also copies frequencies for calculation
      void CopyMeasurementConfigurations(const ThreeDMTModel &Source)
        {
          ClearMeasurementPoints();
          const size_t nmeas = Source.GetMeasPosX().size();
          for (size_t i = 0; i < nmeas; ++i)
            {
              AddMeasurementPoint(Source.GetMeasPosX()[i], Source.GetMeasPosY()[i],
                  Source.GetMeasPosZ()[i]);
            }
          SetFrequencies() = Source.GetFrequencies();
          SetDistortionParameters(Source.GetDistortionParameters());
        }
      //! Write the MT model in VTK format, at the moment the best format for plotting, this only writes the conductivities
      void WriteVTK(const std::string filename) const
        {
          ThreeDModelBase::WriteVTK(filename, "Conductivity");
        }
      //! Read an ascii file as written by ModEM
      void ReadModEM(const std::string filename);
      //! Write an ascii file as read by ModEM
      void WriteModEM(const std::string filename);
      //! We have a copy operator for other MT models
      ThreeDMTModel& operator=(const ThreeDMTModel& source);
      //! Other models will be copied by the copy operator for the base class
      ThreeDMTModel& operator=(const ThreeDModelBase& source);
      ThreeDMTModel();
      //! As we defined a copy operator, we also need a copy constructor
      ThreeDMTModel(const ThreeDMTModel &source);
      virtual ~ThreeDMTModel();
      };

    //! A helper class for the template ThreeDModelObjective that lets us set distortion values as extra inversion parameters
    class J3DEXPORT MTDistortionSetter
      {
    public:
      void operator()(ThreeDMTModel &Model, const std::vector<double> &Dist)
        {
          Model.SetDistortionParameters(Dist);
        }
      };
  /* @} */
  }

#endif /* THREEDMTMODEL_H_ */
