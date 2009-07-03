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
    //! This class stores all information associated with 3D magnetotelluric models, so far it has only rudimentary functionality and will be extended
    class ThreeDMTModel: public jiba::ThreeDModelBase
      {
    private:
      std::vector<double> Frequencies;
    public:
      const std::vector<double> &GetFrequencies()
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
      ThreeDMTModel();
      virtual ~ThreeDMTModel();
      };

  }

#endif /* THREEDMTMODEL_H_ */
