//============================================================================
// Name        : X3DMTCalculator.cpp
// Author      : Jul 7, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================
#ifdef HAVEHPX
#include <hpx/config.hpp>
#include <hpx/include/lcos.hpp>
#endif
#include <cassert>
#include <complex>
#include <fstream>
#include <iomanip>
#include <map>
#ifdef HAVEOPENMP
#include <omp.h>
#endif

#include <boost/multi_array.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/algorithm/string.hpp>

#include "../Global/FatalException.h"
#include "../Global/convert.h"

#include "X3DMTCalculator.h"
#include "X3DFreqFunctions.h"
#include "ReadWriteX3D.h"
#include "MTEquations.h"
#include "InterpolateField.h"
#include "MTUtils.h"

namespace fs = boost::filesystem;

namespace jif3D
  {

    void X3DMTCalculator::CleanUp()
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

    X3DMTCalculator::X3DMTCalculator(boost::filesystem::path TDir, std::string x3d,
        bool DC) :
        X3DName(x3d), WantDistCorr(DC)
      {
        NameRoot = ObjectID();
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

    X3DMTCalculator::~X3DMTCalculator()
      {
        CleanUp();
      }

    rvec X3DMTCalculator::Calculate(const X3DModel &Model, size_t minfreqindex,
        size_t maxfreqindex)
      {
        //we define nfreq as int to make the compiler happy in the openmp loop
        assert(minfreqindex <= maxfreqindex);
        maxfreqindex = std::min(maxfreqindex, Model.GetFrequencies().size());
        const size_t nfreq = maxfreqindex - minfreqindex;
        const size_t nmeas = Model.GetMeasPosX().size();

        //result will hold the final impedance values with
        //applied distortion correction
        jif3D::rvec result(nmeas * nfreq * 8);
        bool FatalError = false;
        result.clear();
        //we store the undistorted impedance with the calculator object
        //as we need it later for the gradient calculation
        //and the distortion correction
        RawImpedance.resize(nmeas * nfreq * 8);
        RawImpedance.clear();
        //we make a call to the coordinate functions to make sure
        //that we have updated the coordinate information and cached it
        //only then the subsequent calls are thread safe
        Model.GetXCoordinates();
        Model.GetYCoordinates();
        Model.GetZCoordinates();

        std::vector<double> C(Model.GetDistortionParameters());
        if (C.size() != nmeas * 4)
          {
            C.resize(nmeas * 4);
            for (size_t i = 0; i < nmeas; ++i)
              {
                C[i * 4] = 1.0;
                C[i * 4 + 1] = 0.0;
                C[i * 4 + 2] = 0.0;
                C[i * 4 + 3] = 1.0;
              }
          }

        std::vector<double> BGDepths(Model.GetBackgroundThicknesses().size(), 0.0);
        std::partial_sum(Model.GetBackgroundThicknesses().begin(),
            Model.GetBackgroundThicknesses().end(), BGDepths.begin());
        CompareDepths(BGDepths, Model.GetZCoordinates());
#ifdef HAVEOPENMP
        omp_lock_t lck;
        omp_init_lock(&lck);

        //we parallelize by frequency, this is relatively simple
        //but we can use up to 20 processors for typical MT problems
        // as we do not have the source for x3d, this is our only possibility anyway
        //the const qualified variables above are predetermined to be shared by the openmp standard
#pragma omp parallel for shared(result) schedule(dynamic,1)
        for (int i = minfreqindex; i < maxfreqindex; ++i)
          {
            //the openmp standard specifies that we cannot leave a parallel construct
            //by throwing an exception, so we catch all exceptions and just
            //generate an error message
            try
              {
                ForwardInfo Info;
                Info.Model = Model;
                Info.C = C;
                Info.freqindex = i;
                Info.TempDirName = TempDir.string();
                Info.NameRoot = NameRoot;
                Info.X3DName = X3DName;
                ForwardResult freqresult = CalculateFrequency(Info);
                size_t startindex = nmeas * i * 8;

                omp_set_lock(&lck);

                std::copy(freqresult.DistImpedance.begin(), freqresult.DistImpedance.end(),
                    result.begin() + startindex);
                std::copy(freqresult.RawImpedance.begin(),freqresult.RawImpedance.end(),RawImpedance.begin()+startindex);

                omp_unset_lock(&lck);
              }
            catch (...)
              {
                FatalError = true;
                std::cerr << "Problem in MT forward calculation.";
              }
            //finished with one frequency
          }

        omp_destroy_lock(&lck);
#endif

#ifdef HAVEHPX
        using hpx::lcos::future;
        using hpx::async;
        using hpx::wait_all;
        std::vector<hpx::naming::id_type> localities =
        hpx::find_all_localities();
        std::vector<hpx::lcos::future<ForwardResult>> FreqResult;
        FreqResult.reserve(nfreq);
        CalculateFrequency_action FreqCalc;
        for (int i = minfreqindex; i < maxfreqindex; ++i)
          {
            ForwardInfo Info;
            Info.Model = Model;
            Info.C = C;
            Info.freqindex = i;
            Info.TempDirName = TempDir.string();
            Info.NameRoot = NameRoot;
            Info.X3DName = X3DName;

            hpx::naming::id_type const locality_id = localities.at(i % localities.size());
            //rvec freqresult = CalculateFrequency(Model, i, TempDir);
            FreqResult.push_back(async(FreqCalc, locality_id, Info));

          }
        wait_all(FreqResult);

        for (int i = minfreqindex; i < maxfreqindex; ++i)
          {
            size_t currindex = i - minfreqindex;
            size_t startindex = nmeas * currindex * 8;
            ForwardResult freqresult = FreqResult[currindex].get();

            std::copy(freqresult.DistImpedance.begin(), freqresult.DistImpedance.end(),
                result.begin() + startindex);
            std::copy(freqresult.RawImpedance.begin(),freqresult.RawImpedance.end(),RawImpedance.begin()+startindex);
          }

#endif
        //we cannot throw from within the openmp section so if there was an exception
        //inside the parallel region we set FatalErrror to true and throw a new exception here
        if (FatalError)
          throw jif3D::FatalException("Problem in MT forward calculation.");
        return result;

      }

    rvec X3DMTCalculator::LQDerivative(const X3DModel &Model, const rvec &Misfit,
        size_t minfreqindex, size_t maxfreqindex)
      {
        //we define nfreq as int to make the compiler happy in the openmp loop
        maxfreqindex = std::min(maxfreqindex, Model.GetFrequencies().size());
        const int nfreq = maxfreqindex - minfreqindex;
        //a few commonly used quantities for shorter notation
        const size_t nmodx = Model.GetConductivities().shape()[0];
        const size_t nmody = Model.GetConductivities().shape()[1];
        const size_t nmodz = Model.GetConductivities().shape()[2];
        //the number of measurement sites
        const size_t nmeas = Model.GetMeasPosX().size();
        const size_t nmod = nmodx * nmody * nmodz;
        assert(Misfit.size() == nmeas * nfreq * 8);
        jif3D::rvec Gradient(nmod, 0.0);
        if (WantDistCorr)
          {
            //if we want to adapt the distortion parameters, we put
            //the gradient with respect to the distortion parameters at the end
            Gradient.resize(nmod + 4 * nmeas, 0.0);
          }
        bool FatalError = false;

        //we make a call to the coordinate functions to make sure
        //that we have updated the coordinate information and cached it
        //only then the subsequent calls are thread safe
        Model.GetXCoordinates();
        Model.GetYCoordinates();
        Model.GetZCoordinates();
        //we read the distortion parameters from the model
        std::vector<double> C(Model.GetDistortionParameters());
        //if they have not been set, we use the identity matrix
        //for each station
        if (C.size() != nmeas * 4)
          {
            C.resize(nmeas * 4);
            for (size_t i = 0; i < nmeas; ++i)
              {
                C[i * 4] = 1.0;
                C[i * 4 + 1] = 0.0;
                C[i * 4 + 2] = 0.0;
                C[i * 4 + 3] = 1.0;
              }
          }

#ifdef HAVEOPENMP
        omp_lock_t lck;
        omp_init_lock(&lck);
#endif
        //we parallelize the gradient calculation by frequency
        //see also the comments for the forward calculation
        //here the explicitly shared variable is Gradient
        //all others are predetermined to be shared
#pragma omp parallel for shared(Gradient) schedule(dynamic,1)
        for (int i = minfreqindex; i < maxfreqindex; ++i)
          {
            try
              {
                ForwardInfo Info;
                Info.Model = Model;
                Info.C = C;
                Info.freqindex = i;
                Info.TempDirName = TempDir.string();
                Info.NameRoot = NameRoot;
                Info.X3DName = X3DName;
                //calculate the gradient for each frequency
                rvec tmp = LQDerivativeFreq(Info, Misfit, RawImpedance);
#ifdef HAVEOPENMP
                omp_set_lock(&lck);
#endif
                //the total gradient is the sum over the gradients for each frequency
                boost::numeric::ublas::subrange(Gradient, 0, nmod) += tmp;
#ifdef HAVEOPENMP
                omp_unset_lock(&lck);
#endif
              } catch (...)
              {
                //we cannot throw exceptions that leave the parallel region
                //so we catch everything and set FatalError to true
                //then outside the parallel region we throw a new exception
                FatalError = true;
              }
            //finished with one frequency
          }
#ifdef HAVEOPENMP
        omp_destroy_lock(&lck);
#endif
        if (WantDistCorr)
          {
            //if we want distortion correct, we calculate the gradient with respect
            //to the distortion parameters and copy the values to the end
            //of the gradient
            jif3D::rvec CGrad = AdaptDist(C, RawImpedance, Misfit);
            std::copy(CGrad.begin(), CGrad.end(), Gradient.begin() + nmod);
          }
        //if we had some exception inside the openmp region, we throw
        // a generic error message
        if (FatalError)
          throw jif3D::FatalException("Problem in MT gradient calculation.");

        return 2.0 * Gradient;
      }
  }
