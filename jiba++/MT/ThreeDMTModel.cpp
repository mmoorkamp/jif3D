//============================================================================
// Name        : ThreeDMTModel.cpp
// Author      : Apr 8, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#include "ThreeDMTModel.h"

namespace jiba
  {

    ThreeDMTModel::ThreeDMTModel() :
      Frequencies()
      {

      }

    ThreeDMTModel::~ThreeDMTModel()
      {

      }

    ThreeDMTModel::ThreeDMTModel(const ThreeDMTModel &source) :
      ThreeDModelBase(source), Frequencies(source.Frequencies)
      {

      }

    ThreeDMTModel& ThreeDMTModel::operator=(const ThreeDMTModel& source)
      {
        if (&source != this)
          {
            //apart from copying the contents of the base class
            //we have to copy the vector of frequencies which
            //is the only additional data in this class
            ThreeDModelBase::operator=(source);
            Frequencies.resize(source.Frequencies.size());
            std::copy(source.Frequencies.begin(), source.Frequencies.end(),
                Frequencies.begin());
          }
        return *this;
      }

    ThreeDMTModel& ThreeDMTModel::operator=(const ThreeDModelBase& source)
      {
        if (&source != this)
          {
            ThreeDModelBase::operator=(source);
          }
        return *this;
      }
  }
