//============================================================================
// Name        : ReadWriteSparseMatrix.cpp
// Author      : 1 Jul 2013
// Version     :
// Copyright   : 2013, mm489
//============================================================================

#include <netcdf>
#include <vector>
#include "ReadWriteSparseMatrix.h"
#include "NetCDFTools.h"

using netCDF::NcFile;
using netCDF::NcDim;

namespace jif3D
  {
    const std::string CovDimName = "nvalues";
    const std::string Index1Name = "Index1";
    const std::string Index2Name = "Index2";

    void ReadSparseMatrixFromNetcdf(const std::string &filename, jif3D::comp_mat &CoVar, const std::string &Name)
      {
        NcFile DataFile(filename, NcFile::read);
        std::vector<int> Index1, Index2;
        std::vector<double> Values;
        ReadVec(DataFile, Index1Name, Index1);
        ReadVec(DataFile, Index2Name, Index2);
        ReadVec(DataFile, Name, Values);
        size_t nvals = Values.size();
        for (size_t i = 0; i < nvals; ++i)
          {
            CoVar(Index1[i], Index2[i]) = Values[i];
            CoVar(Index2[i], Index1[i]) = Values[i];
          }
      }

    void WriteSparseMatrixToNetcdf(const std::string &filename,
        const jif3D::comp_mat &CoVar, const std::string &Name)
      {
        NcFile DataFile(filename, NcFile::replace);

        typedef jif3D::comp_mat::const_iterator1 i1_t;
        typedef jif3D::comp_mat::const_iterator2 i2_t;

        std::vector<int> Index1, Index2;
        std::vector<double> Values;

        for (i1_t i1 = CoVar.begin1(); i1 != CoVar.end1(); ++i1)
          {
            for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2)
              {
                Index1.push_back(i2.index1());
                Index2.push_back(i2.index2());
                Values.push_back(*i2);
              }
          }
        const size_t nstats = Values.size();

        NcDim CovNumDim = DataFile.addDim(CovDimName.c_str(), nstats);
        WriteVec(DataFile, Index1Name, Index1, CovNumDim, "");
        WriteVec(DataFile, Index2Name, Index2, CovNumDim, "");
        WriteVec(DataFile, Name, Values, CovNumDim, "");
      }

  }
