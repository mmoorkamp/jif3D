//============================================================================
// Name        : ThreeDMTModel.h
// Author      : Apr 8, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef THREEDMTMODEL_H_
#define THREEDMTMODEL_H_

#include "../ModelBase/ThreeDModelBase.h"

namespace jiba
  {
    /** \addtogroup mtmodelling Forward modelling of magnetotelluric data */
    /* @{ */
    //! This class stores all information associated with 3D magnetotelluric models
    /*! This class extends ThreeDModelBase to store the calculation frequencies
     * which are required for any MT forward calculation. It also provides named
     * access to the model through SetConductivities and writting file for vtk.
     * All other functionality and properties that are specific to a certain
     * forward code have to be implemented in a derived class.
     */
    class ThreeDMTModel: public jiba::ThreeDModelBase
      {
    private:
      //! The calculation frequencies in Hz
      std::vector<double> Frequencies;
    public:
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
      //! Write the MT model in VTK format, at the moment the best format for plotting, this only writes the conductivities
      void WriteVTK(const std::string filename) const
        {
          ThreeDModelBase::WriteVTK(filename, "Conductivity");
        }
      //! We have a copy operator for other MT models
      ThreeDMTModel& operator= (const ThreeDMTModel& source);
      //! Other models will be copied by the copy operator for the base class
      ThreeDMTModel& operator= (const ThreeDModelBase& source);
      ThreeDMTModel();
      //! As we defined a copy operator, we also need a copy constructor
      ThreeDMTModel(const ThreeDMTModel &source);
      virtual ~ThreeDMTModel();
      };
  /* @} */
  }

#endif /* THREEDMTMODEL_H_ */
