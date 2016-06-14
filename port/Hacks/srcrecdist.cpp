//============================================================================
// Name        : srcrecdist.cpp
// Author      : Oct 11, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================

#include "../Global/VecMat.h"
#include "../Global/FileUtil.h"
#include "../Global/NetCDFTools.h"
#include "../ModelBase/NetCDFModelTools.h"
#include "../ModelBase/VTKTools.h"

int main()
  {

    static const std::string SourceNumberName = "SourceNumber";
    static const std::string ReceiverNumberName = "ReceiverNumber";
    static const std::string SourcePosXName = "SourcePosX";
    static const std::string SourcePosYName = "SourcePosY";
    static const std::string SourcePosZName = "SourcePosZ";
    static const std::string MeasIndexName = "MeasIndex";
    static const std::string TravelTimeName = "TravelTime";
    static const std::string SourceIndexName = "SourceIndex";
    static const std::string ReceiverIndexName = "ReceiverIndex";
    static const std::string MeasPosXName = "MeasPosX";
    static const std::string MeasPosYName = "MeasPosY";
    static const std::string MeasPosZName = "MeasPosZ";

    std::string filename = jif3D::AskFilename("Filename: ");

    NcFile DataFile(filename.c_str(), NcFile::ReadOnly);
    std::vector<double> SourcePosX, SourcePosY, SourcePosZ;

    //read the positions of the sources
    jif3D::ReadVec(DataFile, SourcePosXName, SourcePosX);
    jif3D::ReadVec(DataFile, SourcePosYName, SourcePosY);
    jif3D::ReadVec(DataFile, SourcePosZName, SourcePosZ);


    //read the positions of the receivers
    std::vector<double> RecPosX, RecPosY, RecPosZ;
    jif3D::ReadVec(DataFile, MeasPosXName, RecPosX);
    jif3D::ReadVec(DataFile, MeasPosYName, RecPosY);
    jif3D::ReadVec(DataFile, MeasPosZName, RecPosZ);


    std::vector<int> SourceIndices, ReceiverIndices;
    //now read the indices for the source receiver combinations
    //for each measurement
    jif3D::ReadVec(DataFile, SourceIndexName, SourceIndices);
    jif3D::ReadVec(DataFile, ReceiverIndexName, ReceiverIndices);

    jif3D::rvec Data;
    jif3D::ReadVec(DataFile, TravelTimeName, Data);
    const size_t ntimes = Data.size();

    std::ofstream outfile((filename + ".diff.out").c_str());
    std::ofstream pseudofile((filename + ".xyz").c_str());

    double refx = 0.0;
    double refy = 0.0;
    std::cout << " X Reference coordinate: ";
    std::cin >> refx;
    std::cout << " Y Reference coordinate: ";
    std::cin >> refy;
    std::vector<double> PosX(ntimes), PosY(ntimes), PosZ(ntimes);
    for (size_t i = 0; i < ntimes; ++i)
      {
        double XDist = SourcePosX[SourceIndices[i]] - RecPosX[ReceiverIndices[i]];
        double YDist = SourcePosY[SourceIndices[i]] - RecPosY[ReceiverIndices[i]];
        double Dist = sqrt(XDist * XDist + YDist * YDist);
        outfile << Dist << " " << Data(i) << "\n";
        PosX[i] = SourcePosX[SourceIndices[i]] - XDist / 2.0 - refx;
        PosY[i] = SourcePosY[SourceIndices[i]] - YDist / 2.0 - refy;
        PosZ[i] = Dist;
        double coord = sqrt(PosX[i] * PosX[i] + PosY[i] * PosY[i]);
        pseudofile << coord << " " << Dist << " " << Data(i) << "\n";
      }
    jif3D::Write3DDataToVTK(filename + ".vtk", "Traveltimes", Data, PosX, PosY, PosZ);

  }
