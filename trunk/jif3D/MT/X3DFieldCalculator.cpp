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
#include <chrono>
#include <omp.h>

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

    const std::vector<std::complex<double> > &X3DFieldCalculator::ReturnField(double Freq,
        const std::vector<std::vector<std::complex<double>>> &Field) const
      {
        size_t index = FrequencyMap.at(Freq);
        return Field.at(index);
      }

    void X3DFieldCalculator::CalculateFields(const X3DModel &Model,
        const std::vector<double> &Frequencies, const std::vector<double> MeasPosZ)
      {
        const size_t nfreq = Frequencies.size();
        for (size_t i = 0; i < nfreq; ++i)
          {
            size_t nmap = FrequencyMap.size();
            //insert is only successful when the frequency does not already exist
            auto HasInserted = FrequencyMap.insert(
                std::pair<double, int>(Frequencies.at(i), nmap));
            //insert returns a pair, the second value is a bool
            //which is true when successfull
            if (HasInserted.second)
              {
                HaveCurrentFields.push_back(false);
              }
          }
        size_t ncalcfreq = FrequencyMap.size();
        std::vector<double> CalcFreqs(ncalcfreq);
        for (auto FreqPair : FrequencyMap)
          {
            CalcFreqs.at(FreqPair.second) = FreqPair.first;
          }

        if (!(OldModel == Model))
          {
            HaveCurrentFields.resize(ncalcfreq);
            std::fill(HaveCurrentFields.begin(), HaveCurrentFields.end(), false);
          }
        if (Ex1.size() != ncalcfreq)
          {
            ForwardDirName.resize(ncalcfreq);
            Ex1.resize(ncalcfreq);
            Ex2.resize(ncalcfreq);
            Ey1.resize(ncalcfreq);
            Ey2.resize(ncalcfreq);
            Hx1.resize(ncalcfreq);
            Hx2.resize(ncalcfreq);
            Hy1.resize(ncalcfreq);
            Hy2.resize(ncalcfreq);
            Hz1.resize(ncalcfreq);
            Hz2.resize(ncalcfreq);
            Ex1_all.resize(ncalcfreq);
            Ex2_all.resize(ncalcfreq);
            Ey1_all.resize(ncalcfreq);
            Ey2_all.resize(ncalcfreq);
            Ez1_all.resize(ncalcfreq);
            Ez2_all.resize(ncalcfreq);
          }

        if (ForwardExecTime.empty() || ForwardExecTime.size() != ncalcfreq)
          {
            ForwardExecTime.clear();
            for (size_t i = 0; i < ncalcfreq; ++i)
              {
                ForwardExecTime.push_back(std::make_pair(0, i));
              }
          }

        std::vector<std::pair<size_t, size_t>> NewExecTime;
        omp_lock_t lck;
        omp_init_lock(&lck);

        //we parallelize by frequency, this is relatively simple
        //but we can use up to 20 processors for typical MT problems
        // as we do not have the source for x3d, this is our only possibility anyway
        //the const qualified variables above are predetermined to be shared by the openmp standard
#pragma omp parallel for default(shared) schedule(dynamic,1)
        for (size_t i = 0; i < ncalcfreq; ++i)
          {
            const size_t queueindex = (i % 2) == 0 ? i / 2 : ncalcfreq - 1 - i / 2;
            const size_t calcindex = ForwardExecTime.at(queueindex).second;

            if (!HaveCurrentFields.at(calcindex))
              {
                std::chrono::system_clock::time_point start =
                    std::chrono::system_clock::now();
                CalculateFields(Model, CalcFreqs, MeasPosZ, calcindex);
                HaveCurrentFields.at(calcindex) = true;
                std::chrono::system_clock::time_point end =
                    std::chrono::system_clock::now();
                size_t duration = std::chrono::duration_cast<std::chrono::seconds>(
                    end - start).count();

                omp_set_lock(&lck);
                NewExecTime.push_back(std::make_pair(duration, calcindex));
                omp_unset_lock(&lck);

              }
          }
        for (auto time : NewExecTime)
          {
            ForwardTimesFile << time.second << " " << time.first << " \n";
          }
        ForwardTimesFile << "\n" << std::endl;
        std::stable_sort(NewExecTime.begin(), NewExecTime.end());
        ForwardExecTime = NewExecTime;
        omp_destroy_lock(&lck);

        OldModel = Model;
      }

    void X3DFieldCalculator::CalculateFields(const X3DModel &Model,
        const std::vector<double> &Frequencies, const std::vector<double> MeasPosZ,
        size_t freqindex)
      {

        const size_t nmodx = Model.GetData().shape()[0];
        const size_t nmody = Model.GetData().shape()[1];
        const size_t nmodz = Model.GetData().shape()[2];

        fs::path RootName = TempDir / MakeUniqueName(NameRoot, X3DModel::MT, freqindex);
        fs::path DirName = RootName.string() + dirext;
        ForwardDirName.at(freqindex) = DirName.string();
        std::vector<double> CurrFreq(1, Frequencies[freqindex]);
        std::vector<double> ShiftDepth;
        std::vector<size_t> MeasDepthIndices;
        //construct a vector of indices of unique station depths
        size_t nlevels = ConstructDepthIndices(MeasDepthIndices, ShiftDepth, Model,
            MeasPosZ);
        const size_t nval = (nmodx * nmody * nlevels);
        //std::cout << "CurrFreq: " << CurrFreq.at(0) << " Freqindex: " << freqindex << std::endl;
        if (Ex1.at(freqindex).size() == nval && OldModel == Model)
          {
            //we only check Ex1, because all fields have the same size
            // so if we have the right amount of values in the fields
            //and the model is the same as for the previous calculation,
            //we assume that we already have the correct results
            return;
          }
        else
          {
            Ex1.at(freqindex).clear();
            Ex2.at(freqindex).clear();
            Ey1.at(freqindex).clear();
            Ey2.at(freqindex).clear();
            Hx1.at(freqindex).clear();
            Hx2.at(freqindex).clear();
            Hy1.at(freqindex).clear();
            Hy2.at(freqindex).clear();
            Hz1.at(freqindex).clear();
            Hz2.at(freqindex).clear();
            Ex1_all.at(freqindex).clear();
            Ex2_all.at(freqindex).clear();
            Ey1_all.at(freqindex).clear();
            Ey2_all.at(freqindex).clear();
            Ez1_all.at(freqindex).clear();
            Ez2_all.at(freqindex).clear();
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
            ReadEMO((DirName / emoAname).string(), Ex1.at(freqindex), Ey1.at(freqindex),
                Hx1.at(freqindex), Hy1.at(freqindex), Hz1.at(freqindex));
            ReadEMO((DirName / emoBname).string(), Ex2.at(freqindex), Ey2.at(freqindex),
                Hx2.at(freqindex), Hy2.at(freqindex), Hz2.at(freqindex));
          }

#pragma omp critical(gradient_readema)
          {
            ReadEMA((DirName / emaAname).string(), Ex1_all.at(freqindex),
                Ey1_all.at(freqindex), Ez1_all.at(freqindex), nmodx, nmody, nmodz);
            ReadEMA((DirName / emaBname).string(), Ex2_all.at(freqindex),
                Ey2_all.at(freqindex), Ez2_all.at(freqindex), nmodx, nmody, nmodz);
          }
        CheckField(Ex1.at(freqindex), nval);
        CheckField(Ex2.at(freqindex), nval);
        CheckField(Ey1.at(freqindex), nval);
        CheckField(Ey2.at(freqindex), nval);
        CheckField(Hx1.at(freqindex), nval);
        CheckField(Hx2.at(freqindex), nval);
        CheckField(Hy1.at(freqindex), nval);
        CheckField(Hy2.at(freqindex), nval);
        CheckField(Hz1.at(freqindex), nval);
        CheckField(Hz2.at(freqindex), nval);
      }

    X3DFieldCalculator::X3DFieldCalculator(boost::filesystem::path TDir, std::string x3d,
        bool Clean, jif3D::GreenCalcType GS1, jif3D::GreenCalcType GS4) :
        TempDir(TDir), X3DName(x3d), CleanFiles(Clean), GreenStage1(GS1), GreenStage4(GS4)
      {
        NameRoot = ObjectID();
        ForwardTimesFile.open("mtforward" + NameRoot + ".out");
      }

    X3DFieldCalculator::~X3DFieldCalculator()
      {
        ForwardTimesFile.close();
        CleanUp();
      }

  } /* namespace jif3D */
