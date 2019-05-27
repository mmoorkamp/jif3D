/*
 * X3DFieldCalculator.cpp
 *
 *  Created on: 26 Feb 2018
 *      Author: mm489
 */

#include "X3DFieldCalculator.h"
#include "MTUtils.h"
#include "../ModelBase/CellBoundaries.h"
#include "ReadWriteX3D.h"

#include <boost/filesystem/operations.hpp>
#include <boost/algorithm/string.hpp>

namespace fs = boost::filesystem;

namespace jif3D
  {

    void X3DFieldCalculator::CleanUp()
      {
        //under certain conditions we might not be able
        //to delete all files. We don't want the program to stop
        //because of this as we can always delete files afterwards
        //so we pass an error code object to remove_all and ignore the error
        boost::system::error_code ec;
        fs::directory_iterator end_itr; // default construction yields past-the-end
        //go through the directory and delete any file that starts with NameRoot
        for (fs::directory_iterator itr(TempDir); itr != end_itr; ++itr)
          {
            if (boost::algorithm::starts_with(itr->path().filename().string(), NameRoot))
              {
                fs::remove_all(itr->path().filename(), ec);
              }
          }
      }

    void X3DFieldCalculator::CalculateFields(const X3DModel &Model,
        const std::vector<double> &Frequencies, const std::vector<double> MeasPosZ,
        size_t freqindex)
      {

        const size_t nmodx = Model.GetData().shape()[0];
        const size_t nmody = Model.GetData().shape()[1];

        fs::path RootName = TempDir / MakeUniqueName(NameRoot, X3DModel::MT, freqindex);
        fs::path DirName = RootName.string() + dirext;
        ForwardDirName = DirName.string();
        std::vector<double> CurrFreq(1, Frequencies[freqindex]);
        std::vector<double> ShiftDepth;
        std::vector<size_t> MeasDepthIndices;
        //construct a vector of indices of unique station depths
        size_t nlevels = ConstructDepthIndices(MeasDepthIndices, ShiftDepth, Model,
            MeasPosZ);
        const size_t nval = (nmodx * nmody * nlevels);

        if (Ex1.size() == nval && OldModel == Model)
          {
            //we only check Ex1, because all fields have the same size
            // so if we have the right amount of values in the fields
            //and the model is the same as for the previous calculation,
            //we assume that we already have the correct results
            return;
          }
        else
          {
            Ex1.clear();
            Ex2.clear();
            Ey1.clear();
            Ey2.clear();
            Hx1.clear();
            Hx2.clear();
            Hy1.clear();
            Hy2.clear();
            Hz1.clear();
            Hz2.clear();
          }
        //writing out files causes problems in parallel
        // so we make sure it is done one at a time
#pragma omp critical(forward_write_files)
          {
            MakeRunFile(RootName.string(), DirName.string(), X3DName);
            WriteProjectFile(DirName.string(), CurrFreq, X3DModel::MT, resultfilename,
                modelfilename, GreenStage1, GreenStage4);
            Write3DModelForX3D((DirName / modelfilename).string(), Model.GetXCellSizes(),
                Model.GetYCellSizes(), Model.GetZCellSizes(), ShiftDepth,
                Model.GetConductivities(), Model.GetBackgroundConductivities(),
                Model.GetBackgroundThicknesses());
          }
        //run x3d in parallel
        RunX3D(RootName.string());
        //read in the electric and magnetic field at the observe sites
#pragma omp critical(forward_read_emo)
          {
            ReadEMO((DirName / emoAname).string(), Ex1, Ey1, Hx1, Hy1, Hz1);
            ReadEMO((DirName / emoBname).string(), Ex2, Ey2, Hx2, Hy2, Hz2);
          }

        CheckField(Ex1, nval);
        CheckField(Ex2, nval);
        CheckField(Ey1, nval);
        CheckField(Ey2, nval);
        CheckField(Hx1, nval);
        CheckField(Hx2, nval);
        CheckField(Hy1, nval);
        CheckField(Hy2, nval);
        CheckField(Hz1, nval);
        CheckField(Hz2, nval);
        OldModel = Model;
      }

    X3DFieldCalculator::X3DFieldCalculator(boost::filesystem::path TDir, std::string x3d,
        bool Clean, jif3D::GreenCalcType GS1, jif3D::GreenCalcType GS4) :
        TempDir(TDir), X3DName(x3d), CleanFiles(Clean), GreenStage1(GS1), GreenStage4(GS4)
      {
        NameRoot = ObjectID();
      }

    X3DFieldCalculator::~X3DFieldCalculator()
      {
        CleanUp();
      }

  } /* namespace jif3D */
