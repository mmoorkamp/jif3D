//============================================================================
// Name        : joinpicks.cpp
// Author      : Nov 22, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================

#include "../Global/FileUtil.h"
#include <iostream>
#include "../Tomo/TomographyData.h"
#include "../Tomo/ReadWriteTomographyData.h"

int main()
  {
    std::string filename1 = jif3D::AskFilename("First pick file: ");
    std::string filename2 = jif3D::AskFilename("Second pick file: ");

    jif3D::TomographyData PickData1, PickData2;
    jif3D::rvec Times1, Times2, Error1, Error2;
    PickData1.ReadNetCDF(filename1);
    PickData2.ReadNetCDF(filename2);

    const size_t nsources1 = PickData1.GetSourcePosX().size();
    const size_t nmeas1 = PickData1.GetMeasPosX().size();
    const size_t ntimes1 = PickData1.GetData().size();

    const size_t nsources2 = PickData2.GetSourcePosX().size();
    const size_t nmeas2 = PickData2.GetMeasPosX().size();
    const size_t ntimes2 = PickData2.GetData().size();;

    for (size_t i = 0; i < nsources2; ++i)
      {
        PickData1.AddSource(PickData2.GetSourcePosX()[i],
            PickData2.GetSourcePosY()[i], PickData2.GetSourcePosZ()[i]);
      }

    for (size_t i = 0; i < nmeas2; ++i)
      {
        PickData1.AddMeasurementPoint(PickData2.GetMeasPosX()[i],
            PickData2.GetMeasPosY()[i], PickData2.GetMeasPosZ()[i]);
      }

    for (size_t i = 0; i < ntimes2; ++i)
      {
        PickData1.AddMeasurementConfiguration(PickData2.GetSourceIndices()[i]
            + nsources1, PickData2.GetReceiverIndices()[i] + nmeas1);
      }

    std::vector<double> JointTimes(ntimes1 + ntimes2);
    std::copy(PickData1.GetData().begin(),PickData1.GetData().end(), JointTimes.begin());
    std::copy(PickData2.GetData().begin(), PickData2.GetData().end(), JointTimes.begin() + ntimes1);

    std::vector<double> JointError(ntimes1 + ntimes2);
    std::copy(PickData1.GetErrors().begin(), PickData1.GetErrors().end(), JointError.begin());
    std::copy(PickData2.GetErrors().begin(), PickData1.GetErrors().end(), JointError.begin() + ntimes1);

    PickData1.SetDataAndErrors(JointTimes,JointError);

    std::string outfilename = jif3D::AskFilename("Outfile name: ",false);
    PickData1.WriteNetCDF(outfilename);
    PickData1.WriteSourcePoints(outfilename+".sor.vtk");
    PickData1.WriteMeasurementPoints(outfilename+".rec.vtk");
  }
