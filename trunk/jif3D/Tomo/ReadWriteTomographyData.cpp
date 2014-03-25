//============================================================================
// Name        : ReadWriteTomographyData.cpp
// Author      : May 5, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#include<fstream>
#include "ReadWriteTomographyData.h"
#include "../Global/NumUtil.h"
#include "../Global/NetCDFTools.h"
#include "../ModelBase/NetCDFModelTools.h"

namespace jif3D
  {
    //we define these string variables to guarantee that we
    //use the exact same name for reading and writing the nectdf file
    static const std::string SourceNumberName = "SourceNumber";
    static const std::string ReceiverNumberName = "ReceiverNumber";
    static const std::string SourcePosXName = "SourcePosX";
    static const std::string SourcePosYName = "SourcePosY";
    static const std::string SourcePosZName = "SourcePosZ";
    static const std::string MeasIndexName = "MeasIndex";
    static const std::string TravelTimeName = "TravelTime";
    static const std::string TravelTimeErrorName = "dT";
    static const std::string SourceIndexName = "SourceIndex";
    static const std::string ReceiverIndexName = "ReceiverIndex";

    void SaveTraveltimes(const std::string &filename, const jif3D::rvec &Data,
        const jif3D::rvec &Error, const jif3D::ThreeDSeismicModel &Model)
      {
        //make sure all vectors have consistent sizes
        const size_t ndata = Data.size();
        assert(ndata == Model.GetReceiverIndices().size());
        assert(ndata == Model.GetSourceIndices().size());

        //create a netcdf file
        const size_t nsourcepos = Model.GetSourcePosX().size();
        const size_t nrecpos = Model.GetMeasPosX().size();
        NcFile DataFile(filename.c_str(), NcFile::Replace);
        //we use the station number as a dimension
        NcDim *SourceNumDim = DataFile.add_dim(SourceNumberName.c_str(), nsourcepos);
        NcDim *RecNumDim = DataFile.add_dim(ReceiverNumberName.c_str(), nrecpos);

        //this is just an index over the measurement vector
        //and does not have any special meaning
        std::vector<int> SourcePosNumber, ReceiverPosNumber;
        std::generate_n(std::back_inserter(SourcePosNumber), nsourcepos, IntSequence(0));
        std::generate_n(std::back_inserter(ReceiverPosNumber), nrecpos, IntSequence(0));
        NcVar *SourceNumVar = DataFile.add_var(SourceNumberName.c_str(), ncInt,
            SourceNumDim);
        SourceNumVar->put(&SourcePosNumber[0], nsourcepos);
        //write out the measurement coordinates
        WriteVec(DataFile, SourcePosXName, Model.GetSourcePosX(), SourceNumDim, "m");
        WriteVec(DataFile, SourcePosYName, Model.GetSourcePosY(), SourceNumDim, "m");
        WriteVec(DataFile, SourcePosZName, Model.GetSourcePosZ(), SourceNumDim, "m");
        //write out the positions of the receivers, i.e. measurement positions
        NcVar *RecNumVar = DataFile.add_var(ReceiverNumberName.c_str(), ncInt, RecNumDim);
        RecNumVar->put(&ReceiverPosNumber[0], nrecpos);
        WriteVec(DataFile, MeasPosXName, Model.GetMeasPosX(), RecNumDim, "m");
        WriteVec(DataFile, MeasPosYName, Model.GetMeasPosY(), RecNumDim, "m");
        WriteVec(DataFile, MeasPosZName, Model.GetMeasPosZ(), RecNumDim, "m");
        //generate an index for the receivers
        NcDim *MeasIndexDim = DataFile.add_dim(MeasIndexName.c_str(), ndata);
        std::vector<int> MeasIndex;
        std::generate_n(back_inserter(MeasIndex), ndata, IntSequence(0));
        NcVar *MeasIndexVar = DataFile.add_var(MeasIndexName.c_str(), ncInt,
            MeasIndexDim);
        MeasIndexVar->put(&MeasIndex[0], ndata);
        //Write the travel times and the error
        WriteVec(DataFile, TravelTimeName, Data, MeasIndexDim, "s");

        //write the index of the source for each measurement
        NcVar *SourceIndexVar = DataFile.add_var(SourceIndexName.c_str(), ncInt,
            MeasIndexDim);
        SourceIndexVar->put(&Model.GetSourceIndices()[0], MeasIndexDim->size());
        //write the index of the receiver for each measurement
        NcVar *RecIndexVar = DataFile.add_var(ReceiverIndexName.c_str(), ncInt,
            MeasIndexDim);
        RecIndexVar->put(&Model.GetReceiverIndices()[0], MeasIndexDim->size());

        WriteVec(DataFile, TravelTimeErrorName, Error, MeasIndexDim, "s");

      }

    void ReadTraveltimes(const std::string &filename, jif3D::rvec &Data,
        jif3D::rvec &Error, jif3D::ThreeDSeismicModel &Model)
      {
        NcFile DataFile(filename.c_str(), NcFile::ReadOnly);
        jif3D::ThreeDSeismicModel::tMeasPosVec PosX, PosY, PosZ;
        //delete any old values in the model object
        Model.ClearMeasurementPoints();
        Model.ClearSourcePos();
        Model.ClearMeasurementConfigurations();
        //read the positions of the sources
        ReadVec(DataFile, SourcePosXName, PosX);
        ReadVec(DataFile, SourcePosYName, PosY);
        ReadVec(DataFile, SourcePosZName, PosZ);
        const size_t nsource = PosX.size();
        //and add them to the model object
        for (size_t i = 0; i < nsource; ++i)
          {
            Model.AddSource(PosX[i], PosY[i], PosZ[i]);
          }
        //read the positions of the receivers
        ReadVec(DataFile, MeasPosXName, PosX);
        ReadVec(DataFile, MeasPosYName, PosY);
        ReadVec(DataFile, MeasPosZName, PosZ);
        const size_t nmeas = PosX.size();
        //and add them as measurement positions to the model object
        for (size_t i = 0; i < nmeas; ++i)
          {
            Model.AddMeasurementPoint(PosX[i], PosY[i], PosZ[i]);
          }
        std::vector<int> SourceIndices, ReceiverIndices;
        //now read the indices for the source receiver combinations
        //for each measurement
        ReadVec(DataFile, SourceIndexName, SourceIndices);
        ReadVec(DataFile, ReceiverIndexName, ReceiverIndices);
        const size_t nconf = SourceIndices.size();
        //and configure the model object for these combinations
        for (size_t i = 0; i < nconf; ++i)
          {
            Model.AddMeasurementConfiguration(SourceIndices[i], ReceiverIndices[i]);
          }
        //finally read in the traveltimes
        ReadVec(DataFile, TravelTimeName, Data);
        // it is possible that there is no error information in the file
        // so we don't want to crash the program if not but set it to zero
        NcError NetCDFError(NcError::silent_nonfatal);
        if (DataFile.get_var(TravelTimeErrorName.c_str()) != nullptr)
          {
            ReadVec(DataFile, TravelTimeErrorName, Error);
          }
        else
          {
            Error.resize(Data.size());
            Error.clear();
          }
      }

    void PlotRaypath(const std::string &filename, std::vector<RP_STRUCT> &raypath,
        const size_t nmeas, const double gridspacing, const size_t nairlayers,
        int minxindex, int minyindex)
      {
        std::ofstream outfile(filename.c_str());
        //write out the old style vtk header
        outfile << "# vtk DataFile Version 2.0\n";
        outfile << "Raypaths\n";
        outfile << "ASCII\n";
        outfile << "DATASET POLYDATA\n";
        size_t npoints = 0;
        size_t act_rays = 0;
        //count the number of points we need to plot the rays
        for (size_t i = 0; i < nmeas; ++i)
          {
            if (raypath[i].nray > 0)
              {
                ++act_rays;
                npoints += raypath[i].nray + 1;
              }

          }
        outfile << "POINTS " << npoints << " double\n";
        //write out the positions of each ray segment start and endpoint
        for (size_t i = 0; i < nmeas; ++i)
          {
            if (raypath[i].nray > 0)
              {
                for (size_t j = 0; j < raypath[i].nray + 1; ++j)
                  {
                    outfile << (raypath[i].x[j] + minxindex) * gridspacing << " "
                        << (raypath[i].y[j] + minyindex) * gridspacing << " "
                        << (raypath[i].z[j] - nairlayers) * gridspacing << "\n ";
                  }
              }
          }
        //now we connect the points by lines
        //by writing out the indices of the start and endpoint
        outfile << "\nLINES " << act_rays << " " << npoints + act_rays << std::endl;
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
        //finally we write out the receiver positions
        //as an independent set of points for plotting
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
