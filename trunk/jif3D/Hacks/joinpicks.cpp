//============================================================================
// Name        : joinpicks.cpp
// Author      : Nov 22, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================

#include "../Global/FileUtil.h"
#include <iostream>
#include "../Tomo/ThreeDSeismicModel.h"
#include "../Tomo/ReadWriteTomographyData.h"

int main()
  {
    std::string filename1 = jif3D::AskFilename("First pick file: ");
    std::string filename2 = jif3D::AskFilename("Second pick file: ");

    jif3D::ThreeDSeismicModel PickModel1, PickModel2;
    jif3D::rvec Times1, Times2;
    jif3D::ReadTraveltimes(filename1, Times1, PickModel1);
    jif3D::ReadTraveltimes(filename2, Times2, PickModel2);

    const size_t nsources1 = PickModel1.GetSourcePosX().size();
    const size_t nmeas1 = PickModel1.GetMeasPosX().size();
    const size_t ntimes1 = Times1.size();

    const size_t nsources2 = PickModel2.GetSourcePosX().size();
    const size_t nmeas2 = PickModel2.GetMeasPosX().size();
    const size_t ntimes2 = Times2.size();

    for (size_t i = 0; i < nsources2; ++i)
      {
        PickModel1.AddSource(PickModel2.GetSourcePosX()[i],
            PickModel2.GetSourcePosY()[i], PickModel2.GetSourcePosZ()[i]);
      }

    for (size_t i = 0; i < nmeas2; ++i)
      {
        PickModel1.AddMeasurementPoint(PickModel2.GetMeasPosX()[i],
            PickModel2.GetMeasPosY()[i], PickModel2.GetMeasPosZ()[i]);
      }

    for (size_t i = 0; i < ntimes2; ++i)
      {
        PickModel1.AddMeasurementConfiguration(PickModel2.GetSourceIndices()[i]
            + nsources1, PickModel2.GetReceiverIndices()[i] + nmeas1);
      }

    jif3D::rvec JointTimes(ntimes1 + ntimes2);
    std::copy(Times1.begin(), Times1.end(), JointTimes.begin());
    std::copy(Times2.begin(), Times2.end(), JointTimes.begin() + ntimes1);

    std::string outfilename = jif3D::AskFilename("Outfile name: ",false);
    jif3D::SaveTraveltimes(outfilename, JointTimes, PickModel1);
  }
