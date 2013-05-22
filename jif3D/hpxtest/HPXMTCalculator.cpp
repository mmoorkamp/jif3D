//============================================================================
// Name        : X3DMTCalculator.cpp
// Author      : Jul 7, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================

#include <cassert>
#include <complex>
#include <fstream>
#include <iomanip>

#include <boost/multi_array.hpp>

#include "../Global/FatalException.h"
#include "../Global/convert.h"

#include "HPXMTCalculator.h"
#include "CalcFreq.h"


namespace fs = boost::filesystem;



namespace jif3D
  {

  //check that the .hnk file for x3d are in a certain directory
  bool CheckHNK(const fs::path &TargetDir)
    {
      return fs::exists(TargetDir / "ndec15.hnk")
          && fs::exists(TargetDir / "ndec20.hnk")
          && fs::exists(TargetDir / "ndec30.hnk")
          && fs::exists(TargetDir / "ndec40.hnk");
    }





    void HPXMTCalculator::CleanUp()
      {
        //std::string NameRoot(ObjectID());
        /*		fs::directory_iterator end_itr; // default construction yields past-the-end
         //go through the directory and delete any file that starts with NameRoot
         for (fs::directory_iterator itr(fs::current_path()); itr != end_itr;
         ++itr)
         {
         if (boost::algorithm::starts_with(itr->path().filename().string(),
         NameRoot))
         {
         fs::remove_all(itr->path().filename());
         }
         }*/
      }

    HPXMTCalculator::HPXMTCalculator(boost::filesystem::path TDir)
      {
        if (!fs::is_directory(TDir))
          throw FatalException("TDir is not a directory: " + TDir.string());
        TempDir = TDir;
        //we make sure that the .hnk files are there
        //this is a common problem and when we check for use later
        //we are inside an openmp thread and swallow all sensible error messages.
        if (!CheckHNK(fs::path()))
          {
            throw jif3D::FatalException("Cannot find .hnk files in current directory! ");
          }
      }

    HPXMTCalculator::~HPXMTCalculator()
      {
        CleanUp();
      }



    void CompareDepths(const std::vector<double> &BGDepths,
        const jif3D::ThreeDModelBase::t3DModelDim &ModelDepths)
      {
        size_t mindex = 0;
        for (size_t i = 0; i < BGDepths.size(); ++i)
          {
            while (mindex < ModelDepths.size() && ModelDepths[mindex] < BGDepths[i])
              {
                ++mindex;
              }
            if (mindex < ModelDepths.size() && ModelDepths[mindex] != BGDepths[i])
              {
                throw jif3D::FatalException(
                    "Depth to background layer: " + jif3D::stringify(BGDepths[i])
                        + " does not match grid cell depth: "
                        + jif3D::stringify(ModelDepths[mindex]));
              }
          }
      }

    rvec HPXMTCalculator::Calculate(const X3DModel &Model, size_t minfreqindex,
        size_t maxfreqindex)
      {
        using hpx::lcos::future;
        using hpx::async;
        using hpx::wait_all;

        assert(minfreqindex <= maxfreqindex);
        const size_t maxindex = std::min(maxfreqindex, Model.GetFrequencies().size());
        const size_t nfreq = maxindex - minfreqindex;
        const size_t nmeas = Model.GetMeasPosX().size();
        const size_t nmodx = Model.GetXCoordinates().size();
        const size_t nmody = Model.GetYCoordinates().size();

        jif3D::rvec result(nmeas * nfreq * 8);
        result.clear();
        //we make a call to the coordinate functions to make sure
        //that we have updated the coordinate information and cached it
        //only then the subsequent calls are thread safe
        Model.GetXCoordinates();
        Model.GetYCoordinates();
        Model.GetZCoordinates();
        std::vector<double> BGDepths(Model.GetBackgroundThicknesses().size(), 0.0);
        std::partial_sum(Model.GetBackgroundThicknesses().begin(),
            Model.GetBackgroundThicknesses().end(), BGDepths.begin());
        CompareDepths(BGDepths, Model.GetZCoordinates());

        hpx::naming::id_type const locality_id = hpx::find_here();
        std::vector<future<jif3D::rvec>> FreqResult;
        FreqResult.reserve(nfreq);
        CalculateFrequency_action FreqCalc;
        for (int i = minfreqindex; i < maxindex; ++i)
          {

            //rvec freqresult = CalculateFrequency(Model, i, TempDir);
            FreqResult.push_back(async(FreqCalc, locality_id, Model, i, TempDir.string()));

          }
        wait_all(FreqResult);
        std::cout << "Nfreq: " << FreqResult.size() << std::endl;
        for (int i = minfreqindex; i < maxindex; ++i)
          {
            size_t currindex = i - minfreqindex;
            size_t startindex = nmeas * currindex * 8;
            jif3D::rvec imp = FreqResult[currindex].get();
            std::cout << imp.size() << std::endl;
            std::copy(imp.begin(), imp.end(), result.begin() + startindex);
          }

        return result;

      }





  }

