/*
 * EMData.h
 *
 *  Created on: May 22, 2019
 *      Author: max
 */

#ifndef MT_EMDATA_H_
#define MT_EMDATA_H_

#include "../Global/Serialization.h"
#include "../DataBase/GeneralData.h"
#include "X3DModel.h"

namespace jif3D
  {

    class EMData: public GeneralData
      {
    private:
      std::vector<double> Frequencies;
      std::vector<double> RotAngles;
      std::vector<std::string> Names;
    public:
      typedef jif3D::X3DModel ModelType;
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & base_object<GeneralData>(*this);
          ar & Frequencies;
          ar & RotAngles;
        }
      void SetFrequencies(const std::vector<double> &Freq)
        {
          Frequencies = Freq;
        }
      const std::vector<double> &GetFrequencies() const
        {
          return Frequencies;
        }
      const std::vector<double> &GetRotAngles() const
        {
          return RotAngles;
        }
      void SetRotAngles(const std::vector<double> &RA)
        {
          RotAngles = RA;
        }
      const std::vector<std::string> &GetNames() const
        {
          return Names;
        }
      void SetNames(const std::vector<std::string> &s)
        {
          Names = s;
        }

      EMData();
      virtual ~EMData();
      };

  } /* namespace jif3D */

#endif /* MT_EMDATA_H_ */
