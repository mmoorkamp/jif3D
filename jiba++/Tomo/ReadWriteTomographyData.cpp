//============================================================================
// Name        : ReadWriteTomographyData.cpp
// Author      : May 5, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#include<fstream>
#include "ReadWriteTomographyData.h"
#include "../Global/NumUtil.h"
#include "../ModelBase/NetCDFTools.h"

namespace jiba
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

    void SaveTraveltimes(const std::string &filename, const jiba::rvec &Data,
        const jiba::ThreeDSeismicModel &Model)
      {
        //make sure all vectors have consistent sizes
        const size_t ndata = Data.size();
        assert( ndata == Model.GetReceiverIndices().size());
        assert(ndata == Model.GetSourceIndices().size());

        //create a netcdf file
        const size_t nsourcepos = Model.GetSourcePosX().size();
        const size_t nrecpos = Model.GetMeasPosX().size();
        NcFile DataFile(filename.c_str(), NcFile::Replace);
        //we use the station number as a dimension
        NcDim *SourceNumDim = DataFile.add_dim(SourceNumberName.c_str(),
            nsourcepos);
        NcDim *RecNumDim =
            DataFile.add_dim(ReceiverNumberName.c_str(), nrecpos);

        //this is just an index over the measurement vector
        //and does not have any special meaning
        std::vector<int> SourcePosNumber, ReceiverPosNumber;
        std::generate_n(std::back_inserter(SourcePosNumber), nsourcepos,
            IntSequence(0));
        std::generate_n(std::back_inserter(ReceiverPosNumber), nrecpos,
            IntSequence(0));
        NcVar *SourceNumVar = DataFile.add_var(SourceNumberName.c_str(), ncInt,
            SourceNumDim);
        SourceNumVar->put(&SourcePosNumber[0], nsourcepos);
        //write out the measurement coordinates
        WriteVec(DataFile, SourcePosXName, Model.GetSourcePosX(), SourceNumDim,
            "m");
        WriteVec(DataFile, SourcePosYName, Model.GetSourcePosY(), SourceNumDim,
            "m");
        WriteVec(DataFile, SourcePosZName, Model.GetSourcePosZ(), SourceNumDim,
            "m");

        NcVar *RecNumVar = DataFile.add_var(ReceiverNumberName.c_str(), ncInt,
            RecNumDim);
        RecNumVar->put(&ReceiverPosNumber[0], nrecpos);
        WriteVec(DataFile, MeasPosXName, Model.GetMeasPosX(), RecNumDim, "m");
        WriteVec(DataFile, MeasPosYName, Model.GetMeasPosY(), RecNumDim, "m");
        WriteVec(DataFile, MeasPosZName, Model.GetMeasPosZ(), RecNumDim, "m");

        NcDim *MeasIndexDim = DataFile.add_dim(MeasIndexName.c_str(), ndata);
        std::vector<int> MeasIndex;
        std::generate_n(back_inserter(MeasIndex), ndata, IntSequence(0));
        NcVar *MeasIndexVar = DataFile.add_var(MeasIndexName.c_str(), ncInt,
            MeasIndexDim);
        MeasIndexVar->put(&MeasIndex[0], ndata);
        //Write the measurements
        NcVar *DataVar = DataFile.add_var(TravelTimeName.c_str(), ncDouble,
            MeasIndexDim);
        DataVar->add_att("units", "s");
        DataVar->add_att("_FillValue", -1.0);
        DataVar->put(&Data[0], MeasIndexDim->size());

        NcVar *SourceIndexVar = DataFile.add_var(SourceIndexName.c_str(),
            ncInt, MeasIndexDim);
        SourceIndexVar->put(&Model.GetSourceIndices()[0], MeasIndexDim->size());

        NcVar *RecIndexVar = DataFile.add_var(ReceiverIndexName.c_str(), ncInt,
            MeasIndexDim);
        RecIndexVar->put(&Model.GetReceiverIndices()[0], MeasIndexDim->size());

      }

    void ReadTraveltimes(const std::string &filename, jiba::rvec &Data,
        jiba::ThreeDSeismicModel &Model)
      {
        NcFile DataFile(filename.c_str(), NcFile::ReadOnly);
        jiba::ThreeDSeismicModel::tMeasPosVec PosX, PosY, PosZ;
        Model.ClearMeasurementPoints();
        Model.ClearSourcePos();
        Model.ClearMeasurementConfigurations();
        ReadVec(DataFile, SourcePosXName, SourceNumberName, PosX);
        ReadVec(DataFile, SourcePosYName, SourceNumberName, PosY);
        ReadVec(DataFile, SourcePosZName, SourceNumberName, PosZ);
        const size_t nsource = PosX.size();
        for (size_t i = 0; i < nsource; ++i)
          {
            Model.AddSource(PosX[i], PosY[i], PosZ[i]);
          }
        ReadVec(DataFile, MeasPosXName, ReceiverNumberName, PosX);
        ReadVec(DataFile, MeasPosYName, ReceiverNumberName, PosY);
        ReadVec(DataFile, MeasPosZName, ReceiverNumberName, PosZ);
        const size_t nmeas = PosX.size();
        for (size_t i = 0; i < nmeas; ++i)
          {
            Model.AddMeasurementPoint(PosX[i], PosY[i], PosZ[i]);
          }
        std::vector<int> SourceIndices, ReceiverIndices;
        ReadVec(DataFile, SourceIndexName, MeasIndexName, SourceIndices);
        ReadVec(DataFile, ReceiverIndexName, MeasIndexName, ReceiverIndices);
        const size_t nconf = SourceIndices.size();
        for (size_t i = 0; i < nconf; ++i)
          {
            Model.AddMeasurementConfiguration(SourceIndices[i],
                ReceiverIndices[i]);
          }
        ReadVec(DataFile, TravelTimeName,MeasIndexName, Data);
      }


    void PlotRaypath(const std::string &filename, jiba::RP_STRUCT *raypath,
        const size_t nmeas, const double gridspacing, const size_t nairlayers)
      {
        std::ofstream outfile(filename.c_str());
        outfile << "# vtk DataFile Version 2.0\n";
        outfile << "Raypaths\n";
        outfile << "ASCII\n";
        outfile << "DATASET POLYDATA\n";
        size_t npoints = 0;
        size_t act_rays = 0;
        for (size_t i = 0; i < nmeas; ++i)
          {
            if (raypath[i].nray > 0)
              {
                ++act_rays;
                npoints += raypath[i].nray + 1;
              }

          }
        outfile << "POINTS " << npoints << " double\n";

        for (size_t i = 0; i < nmeas; ++i)
          {
            if (raypath[i].nray > 0)
              {
                for (size_t j = 0; j < raypath[i].nray + 1; ++j)
                  {
                    outfile << raypath[i].x[j]  * gridspacing << " "
                        << raypath[i].y[j] * gridspacing << " "
                        << (raypath[i].z[j]- nairlayers) * gridspacing << "\n ";
                  }
              }
          }
        outfile << "\nLINES " << act_rays << " " << npoints + act_rays
            << std::endl;
        size_t index = 0;
        for (size_t i = 0; i < nmeas; ++i)
          {
            if (raypath[i].nray > 0)
              {
                outfile << raypath[i].nray + 1;
                for (size_t j = 0; j < raypath[i].nray + 1; ++j)
                  {
                    outfile << " " << index;
                    ++index;
                  }
                outfile << std::endl;
              }
          }
        outfile << "POINT_DATA " << npoints << std::endl;
        outfile << "SCALARS  Rays float\n";
        outfile << "LOOKUP_TABLE default\n";
        for (size_t i = 0; i < nmeas; ++i)
          {
            if (raypath[i].nray > 0)
              {
                std::fill_n(std::ostream_iterator<double>(outfile, "\n"),
                    raypath[i].nray + 1, double(i));
              }
          }
      }
  }
