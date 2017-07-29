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
#include <chrono>
#ifdef HAVEOPENMP
#include <omp.h>
#endif

//#include <boost/multi_array.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/numeric/conversion/cast.hpp>

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
        bool DC, bool Clean) :
        GreenType1(hst), GreenType4(hst), X3DName(x3d), WantDistCorr(DC), CleanFiles(
            Clean)
      {
        //each object gets a unique ID, this way we avoid clashes
        //between the temporary files generated for the calculations with x3d
        NameRoot = ObjectID();
        if (!fs::is_directory(TDir))
          throw FatalException("TDir is not a directory: " + TDir.string(), __FILE__,
          __LINE__);
        TempDir = TDir;
        //we make sure that the .hnk files are there
        //this is a common problem and when we check for use later
        //we are inside an openmp thread and swallow all sensible error messages.
        if (!CheckHNK(fs::path()))
          {
            throw jif3D::FatalException("Cannot find .hnk files in current directory! ",
            __FILE__, __LINE__);
          }
      }

    X3DMTCalculator::~X3DMTCalculator()
      {
        //if we want to clean all temporary files (default)
        if (CleanFiles)
          {
            //remove all the temporary files and directories generated for calculations
            CleanUp();
          }
      }

    rvec X3DMTCalculator::Calculate(X3DModel &Model, size_t minfreqindex,
        size_t maxfreqindex)
      {

        //we define nfreq as int to make the compiler happy in the openmp loop
        assert(minfreqindex <= maxfreqindex);
        maxfreqindex = std::min(maxfreqindex, Model.GetFrequencies().size());
        const size_t nfreq = maxfreqindex - minfreqindex;
        const size_t nmeas = Model.GetMeasPosX().size();

        if (ForwardExecTime.empty())
          {
            for (size_t i = 0; i < nfreq; ++i)
              {
                ForwardExecTime.push_back(std::make_pair(0, minfreqindex + i));
              }
          }
        //if the current model does not contain any ExIndices information
        //generate ExIndices, EyIndices and HIndices as 0:nmeas for each Frequency
        //here we assume that we either have all three indices in the netCDF file or none of them
        std::vector<int> ExIndices(Model.GetExIndices()), EyIndices(Model.GetEyIndices()),
            HIndices(Model.GetHIndices());
        size_t ind_shift = 0;
        if (ExIndices.empty())
          {
            ExIndices.resize(nmeas * nfreq);
            EyIndices.resize(nmeas * nfreq);
            HIndices.resize(nmeas * nfreq);
            for (size_t ifr = 0; ifr < nfreq; ++ifr)
              {
                ind_shift = nmeas * ifr;
                for (size_t i = 0; i < nmeas; ++i)
                  {
                    ExIndices[i + ind_shift] = i;
                  }
              }
            EyIndices = ExIndices;
            HIndices = ExIndices;
          }
        Model.SetFieldIndices(ExIndices, EyIndices, HIndices);

        const size_t nstats = Model.GetExIndices().size() / nfreq;

        std::string ErrorMsg;
        //result will hold the final impedance values with
        //applied distortion correction
        jif3D::rvec result(nstats * nfreq * 8);
        if (!DataTransform)
          {
            DataTransform = boost::make_shared<jif3D::CopyTransform>(result.size());
          }
        bool FatalError = false;
        result.clear();
        //we store the undistorted impedance with the calculator object
        //as we need it later for the gradient calculation
        //and the distortion correction
        RawImpedance.resize(nstats * nfreq * 8);
        RawImpedance.clear();
        //we make a call to the coordinate functions to make sure
        //that we have updated the coordinate information and cached it
        //only then the subsequent calls are thread safe
        Model.GetXCoordinates();
        Model.GetYCoordinates();
        Model.GetZCoordinates();

        //if the current model does not contain any distortion information
        //generate distortion parameters equivalent to an identity matrix
        std::vector<double> C(Model.GetDistortionParameters());
        if (C.size() != nstats * 4)
          {
            C.resize(nstats * 4);
            for (size_t i = 0; i < nstats; ++i)
              {
                C[i * 4] = 1.0;
                C[i * 4 + 1] = 0.0;
                C[i * 4 + 2] = 0.0;
                C[i * 4 + 3] = 1.0;
              }
          }
        //check that the depths to the different background layers match
        //with the depths to grid cell boundaries
        std::vector<double> BGDepths(Model.GetBackgroundThicknesses().size(), 0.0);
        std::partial_sum(Model.GetBackgroundThicknesses().begin(),
            Model.GetBackgroundThicknesses().end(), BGDepths.begin());
        CompareDepths(BGDepths, Model.GetZCoordinates());
        std::vector<std::pair<size_t, size_t>> NewExecTime;
#ifdef HAVEOPENMP
        omp_lock_t lck;
        omp_init_lock(&lck);
        //openmp loop indices have to be int, so we cast out upper limit to int
        //to make the compiler happy
        int maxindex = boost::numeric_cast<int>(maxfreqindex);
        rvec RawImp(RawImpedance.size(),0.0);
        //we parallelize by frequency, this is relatively simple
        //but we can use up to 20 processors for typical MT problems
        // as we do not have the source for x3d, this is our only possibility anyway
        //the const qualified variables above are predetermined to be shared by the openmp standard
#pragma omp parallel for shared(result,RawImp) schedule(dynamic,1)
        for (int i = 0; i < nfreq; ++i)
          {
            //the openmp standard specifies that we cannot leave a parallel construct
            //by throwing an exception, so we catch all exceptions and just
            //generate an error message
            try
              {
                //we want to alternate between items at the beginning of the map and at the end of the map
                const size_t queueindex = (i % 2) == 0 ? i / 2 : nfreq - 1 - i / 2;
                const size_t calcindex = ForwardExecTime.at(queueindex).second;
                ForwardInfo Info(Model,C,calcindex,TempDir.string(),X3DName, NameRoot, GreenType1, GreenType4);
                std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
                ForwardResult freqresult = CalculateFrequency(Info);
                std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
                size_t duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
                NewExecTime.push_back(std::make_pair(duration,calcindex));
                const size_t currindex = calcindex - minfreqindex;
                const size_t startindex = nstats * currindex * 8;

                omp_set_lock(&lck);

                std::copy(freqresult.DistImpedance.begin(), freqresult.DistImpedance.end(),
                    result.begin() + startindex);
                std::copy(freqresult.RawImpedance.begin(),freqresult.RawImpedance.end(),RawImp.begin()+startindex);

                omp_unset_lock(&lck);
              }
            catch (jif3D::FatalException &e)
              {
                ErrorMsg = e.what();
                FatalError = true;
              }
            catch (...)
              {
                FatalError = true;
                std::cerr << "Problem in MT forward calculation.";
              }
            //finished with one frequency
          }
        std::sort(NewExecTime.begin(),NewExecTime.end());
        ForwardExecTime = NewExecTime;
        RawImpedance = RawImp;
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
        for (size_t i = minfreqindex; i < maxfreqindex; ++i)
          {
            ForwardInfo Info(Model,C,i,TempDir.string(),X3DName, NameRoot, GreenType1, GreenType4);
            hpx::naming::id_type const locality_id = localities.at(i % localities.size());
            //std::cout << "Sending frequency: " << i << " to node " << locality_id << std::endl;
            //rvec freqresult = CalculateFrequency(Model, i, TempDir);
            FreqResult.push_back(async(FreqCalc, locality_id, Info));

          }
        wait_all(FreqResult);

        for (size_t i = minfreqindex; i < maxfreqindex; ++i)
          {
            const size_t currindex = i - minfreqindex;
            const size_t startindex = nstats * currindex * 8;
            ForwardResult freqresult = FreqResult[currindex].get();
            //std::cout << " Getting results for frequency " << currindex << std::endl;
            std::copy(freqresult.DistImpedance.begin(), freqresult.DistImpedance.end(),
                result.begin() + startindex);
            std::copy(freqresult.RawImpedance.begin(),freqresult.RawImpedance.end(),RawImpedance.begin()+startindex);
          }

#endif

        //we cannot throw from within the openmp section so if there was an exception
        //inside the parallel region we set FatalErrror to true and throw a new exception here
        if (FatalError)
          throw jif3D::FatalException("Problem in MT forward calculation: " + ErrorMsg,
          __FILE__,
          __LINE__);
        if (DataTransform)
          {
            result = jif3D::ApplyTransform(result, *DataTransform);
          }
        return result;
      }

    rvec X3DMTCalculator::LQDerivative(const X3DModel &Model, const rvec &Misfit,
        size_t minfreqindex, size_t maxfreqindex)
      {
        if (!DataTransform)
          {
            DataTransform = boost::make_shared<jif3D::CopyTransform>(Misfit.size());
          }
        //we define nfreq as int to make the compiler happy in the openmp loop
        maxfreqindex = std::min(maxfreqindex, Model.GetFrequencies().size());
        const int nfreq = maxfreqindex - minfreqindex;
        if (DerivExecTime.empty())
          {
            for (size_t i = 0; i < nfreq; ++i)
              {
                DerivExecTime.push_back(std::make_pair(0, minfreqindex + i));
              }
          }
        std::string ErrorMsg;
        //a few commonly used quantities for shorter notation
        const size_t nmodx = Model.GetConductivities().shape()[0];
        const size_t nmody = Model.GetConductivities().shape()[1];
        const size_t nmodz = Model.GetConductivities().shape()[2];
        //the number of measurement sites
        const size_t nmeas = Model.GetMeasPosX().size();
        const size_t nstats = Model.GetExIndices().size() / nfreq;
        const size_t nmod = nmodx * nmody * nmodz;
        assert(Misfit.size() == nstats * nfreq * 8);
        jif3D::rvec Gradient(nmod);
        if (WantDistCorr)
          {
            //if we want to adapt the distortion parameters, we put
            //the gradient with respect to the distortion parameters at the end
            Gradient.resize(nmod + 4 * nstats);
          }
        bool FatalError = false;
        //we need to initialize all values to zero as we are adding
        //the individual gradients per frequency
        std::fill(Gradient.begin(), Gradient.end(), 0.0);
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
        if (C.size() != nstats * 4)
          {
            C.resize(nstats * 4);
            for (size_t i = 0; i < nstats; ++i)
              {
                C[i * 4] = 1.0;
                C[i * 4 + 1] = 0.0;
                C[i * 4 + 2] = 0.0;
                C[i * 4 + 3] = 1.0;
              }
          }

        jif3D::rvec ProjMisfit(Misfit.size(), 0.0);
        for (size_t i = 0; i < Misfit.size(); i += 2)
          {
            auto GradTrans = DataTransform->Derivative(
                ublas::subrange(RawImpedance, i, i + 2));
            ProjMisfit(i) = GradTrans(0, 0) * Misfit(i) + GradTrans(1, 0) * Misfit(i + 1);
            ProjMisfit(i + 1) = GradTrans(0, 1) * Misfit(i)
                + GradTrans(1, 1) * Misfit(i + 1);
          }
#ifdef HAVEOPENMP
        std::vector<std::pair<size_t, size_t>> NewExecTime;
        omp_lock_t lck;
        omp_init_lock(&lck);
        //openmp loop indices have to be int, so we cast out upper limit to int
        //to make the compiler happy
        int maxindex = boost::numeric_cast<int>(maxfreqindex);
        //we parallelize the gradient calculation by frequency
        //see also the comments for the forward calculation
        //here the explicitly shared variable is Gradient
        //all others are predetermined to be shared
#pragma omp parallel for shared(Gradient) schedule(dynamic,1)
        for (int i = 0; i < nfreq; ++i)
          {
            try
              {
                //we want to alternate between items at the beginning of the map and at the end of the map
                const size_t queueindex = (i % 2) == 0 ? i / 2 : nfreq - 1 - i / 2;
                const size_t calcindex = DerivExecTime.at(queueindex).second;

                ForwardInfo Info(Model,C,calcindex,TempDir.string(),X3DName, NameRoot, GreenType1, GreenType4);
                //calculate the gradient for each frequency
                std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
                GradResult tmp = LQDerivativeFreq(Info, GradInfo(ProjMisfit, RawImpedance));
                std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
                size_t duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
                NewExecTime.push_back(std::make_pair(duration,calcindex));
                omp_set_lock(&lck);
                //the total gradient is the sum over the gradients for each frequency
                std::transform(tmp.Gradient.begin(),tmp.Gradient.end(),
                    Gradient.begin(),Gradient.begin(),std::plus<double>());
                omp_unset_lock(&lck);
              }
            catch (jif3D::FatalException &e)
              {
                ErrorMsg = e.what();
                FatalError = true;
              }
            catch (...)
              {
                //we cannot throw exceptions that leave the parallel region
                //so we catch everything and set FatalError to true
                //then outside the parallel region we throw a new exception
                FatalError = true;
              }
            //finished with one frequency
          }
        std::sort(NewExecTime.begin(),NewExecTime.end());
        DerivExecTime = NewExecTime;
        omp_destroy_lock(&lck);
#endif

#ifdef HAVEHPX
        using hpx::lcos::future;
        using hpx::async;
        using hpx::wait_all;
        std::vector<hpx::naming::id_type> localities =
        hpx::find_all_localities();
        std::vector<hpx::lcos::future<GradResult>> FreqResult;
        FreqResult.reserve(nfreq);
        LQDerivativeFreq_action LQDerivativeFreq;
        for (size_t i = minfreqindex; i < maxfreqindex; ++i)
          {
            ForwardInfo Info(Model,C,i,TempDir.string(),X3DName, NameRoot, GreenType1, GreenType4);

            hpx::naming::id_type const locality_id = localities.at(i % localities.size());
            //rvec freqresult = CalculateFrequency(Model, i, TempDir);
            FreqResult.push_back(async(LQDerivativeFreq, locality_id, Info, GradInfo(ProjMisfit, RawImpedance)));

          }
        wait_all(FreqResult);

        for (int i = 0; i < nfreq; ++i)
          {
            std::vector<double> CurrGrad(FreqResult[i].get().Gradient);
            std::transform(CurrGrad.begin(),CurrGrad.end(),Gradient.begin(),Gradient.begin(),std::plus<double>());
            //boost::numeric::ublas::subrange(Gradient, 0, nmod) += FreqResult[currindex].get().Gradient;

          }

#endif
        if (WantDistCorr)
          {
            //if we want distortion correct, we calculate the gradient with respect
            //to the distortion parameters and copy the values to the end
            //of the gradient
            jif3D::rvec CGrad = AdaptDist(C, RawImpedance, Misfit);
            boost::numeric::ublas::subrange(Gradient, nmod, Gradient.size()) = CGrad;
            //std::copy(CGrad.begin(), CGrad.end(), Gradient.begin() + nmod);
          }
        //if we had some exception inside the openmp region, we throw
        // a generic error message
        if (FatalError)
          throw jif3D::FatalException(
              "Problem in MT gradient calculation. Error message: " + ErrorMsg, __FILE__,
              __LINE__);

        return 2.0 * Gradient;
      }

    rmat X3DMTCalculator::SensitivityMatrix(ModelType &Model, const rvec &Misfit,
        size_t minfreqindex, size_t maxfreqindex)
      {
        if (!DataTransform)
          {
            DataTransform = boost::make_shared<jif3D::CopyTransform>(Misfit.size());
          }
        maxfreqindex = std::min(maxfreqindex, Model.GetFrequencies().size());
        const size_t nfreq = maxfreqindex - minfreqindex;
        const size_t ndata = Misfit.size();
        const size_t nmeas = Model.GetMeasPosX().size();

        const size_t nmodel = Model.GetConductivities().num_elements();
        //const size_t nsites = Model.GetMeasPosX().size();
        const size_t nsites = Model.GetExIndices().size() / nfreq;
        rmat Result;
        if (WantDistCorr)
          {
            //if we want to adapt the distortion parameters, we put
            //the gradient with respect to the distortion parameters at the end
            Result.resize(ndata, nmodel + 4 * nsites);
          }
        else
          {
            Result.resize(ndata, nmodel);
          }
        //set all sensitivities initially to zero
        Result.clear();
        //we read the distortion parameters from the model
        std::vector<double> C(Model.GetDistortionParameters());
        //if they have not been set, we use the identity matrix
        //for each station
        if (C.size() != nsites * 4)
          {
            C.resize(nsites * 4);
            for (size_t i = 0; i < nsites; ++i)
              {
                C[i * 4] = 1.0;
                C[i * 4 + 1] = 0.0;
                C[i * 4 + 2] = 0.0;
                C[i * 4 + 3] = 1.0;
              }
          }

        for (size_t i = 0; i < ndata; ++i)
          {
            size_t freqindex = i / (nsites * 8);
            rvec CurrMisfit(Misfit.size(), 0.0);
            CurrMisfit(i) = 1.0;
            ForwardInfo Info(Model, C, freqindex, TempDir.string(), X3DName, NameRoot,
                GreenType1, GreenType4);
            GradResult CurrGrad = LQDerivativeFreq(Info,
                GradInfo(CurrMisfit, RawImpedance));

            boost::numeric::ublas::matrix_row<rmat> CurrRow(Result, i);

            jif3D::rvec Gradient(CurrGrad.Gradient.size());
            std::copy(CurrGrad.Gradient.begin(), CurrGrad.Gradient.end(),
                Gradient.begin());
            boost::numeric::ublas::subrange(CurrRow, 0, nmodel) = Gradient;
            if (WantDistCorr)
              {
                jif3D::rvec CGrad = AdaptDist(C, RawImpedance, CurrMisfit);
                boost::numeric::ublas::subrange(CurrRow, nmodel, CurrRow.size()) = CGrad;
              }
          }
        return Result;
      }
  }
