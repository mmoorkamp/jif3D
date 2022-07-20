/*
 * X3DTipperCalculator.cpp
 *
 *  Created on: 25 May 2018
 *      Author: mm489
 */
#include "X3DTipperCalculator.h"
#include <complex>
#include <map>
#include <chrono>
#ifdef HAVEOPENMP
#include <omp.h>
#endif
#include <boost/make_shared.hpp>
#include <boost/algorithm/string.hpp>
#include "../Global/ignore.h"
#include "../ModelBase/CellBoundaries.h"
#include "X3DMTCalculator.h"
#include "X3DFreqFunctions.h"
#include "X3DFieldCalculator.h"
#include "ReadWriteX3D.h"
#include "MTEquations.h"
#include "InterpolateField.h"
#include "MTUtils.h"


namespace fs = boost::filesystem;

namespace jif3D
  {

    void X3DTipperCalculator::CleanUp()
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

    X3DTipperCalculator::X3DTipperCalculator(boost::filesystem::path TDir,
        std::string x3d, bool Clean, boost::shared_ptr<jif3D::X3DFieldCalculator> FC) :
        GreenType1(hst), GreenType4(hst), X3DName(x3d), CleanFiles(Clean), FieldCalculator(
            FC)
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
        if (FieldCalculator.get() == nullptr)
          {
            FieldCalculator = boost::make_shared<jif3D::X3DFieldCalculator>(TDir, x3d);
          }
      }

    X3DTipperCalculator::~X3DTipperCalculator()
      {
        CleanUp();
      }

    rvec X3DTipperCalculator::Calculate(const ModelType &Model, const TipperData &Data,
        size_t minfreqindex, size_t maxfreqindex)
      {
        //we define nfreq as int to make the compiler happy in the openmp loop
        assert(minfreqindex <= maxfreqindex);
        maxfreqindex = std::min(maxfreqindex, Data.GetFrequencies().size());

        const size_t nfreq = maxfreqindex - minfreqindex;
        if (ForwardExecTime.empty() || ForwardExecTime.size() != nfreq)
          {
            for (size_t i = 0; i < nfreq; ++i)
              {
                ForwardExecTime.push_back(std::make_pair(0, minfreqindex + i));
              }
          }

        const size_t nstats = Data.GetHxIndices().size() / nfreq;
        const size_t nmodx = Model.GetConductivities().shape()[0];
        const size_t nmody = Model.GetConductivities().shape()[1];


        std::string ErrorMsg;
        //result will hold the final impedance values with
        //applied distortion correction
        jif3D::rvec result(nstats * nfreq * 4);

        bool FatalError = false;
        result.clear();

        //check that the depths to the different background layers match
        //with the depths to grid cell boundaries
        std::vector<double> BGDepths(Model.GetBackgroundThicknesses().size(), 0.0);
        std::partial_sum(Model.GetBackgroundThicknesses().begin(),
            Model.GetBackgroundThicknesses().end(), BGDepths.begin());
        double ZOrigin = Model.GetZCoordinates()[0];
        std::for_each(BGDepths.begin(), BGDepths.end(), [ZOrigin](double& d)
          { d+=ZOrigin;});
        CompareDepths(BGDepths, Model.GetZCoordinates());
        std::vector<std::pair<size_t, size_t>> NewExecTime;
#ifdef HAVEOPENMP
        omp_lock_t lck;
        omp_init_lock(&lck);
#endif
        //we parallelize by frequency, this is relatively simple
        //but we can use up to 20 processors for typical MT problems
        // as we do not have the source for x3d, this is our only possibility anyway
        //the const qualified variables above are predetermined to be shared by the openmp standard
        FieldCalculator->CalculateFields(Model, Data.GetFrequencies(),
            Data.GetMeasPosZ());
#pragma omp parallel for shared(result,NewExecTime) schedule(dynamic,1)
        for (size_t i = 0; i < nfreq; ++i)
          {
            //the openmp standard specifies that we cannot leave a parallel construct
            //by throwing an exception, so we catch all exceptions and just
            //generate an error message
            try
              {
                //we want to alternate between items at the beginning of the map and at the end of the map
                const size_t queueindex = (i % 2) == 0 ? i / 2 : nfreq - 1 - i / 2;
                const size_t calcindex = ForwardExecTime.at(queueindex).second;

                const size_t currindex = calcindex - minfreqindex;
                const size_t startindex = nstats * currindex * 4;
                const size_t ind_shift = nstats * calcindex;

                std::chrono::system_clock::time_point start =
                    std::chrono::system_clock::now();

                std::vector<double> ShiftDepth;
                std::vector<size_t> MeasDepthIndices;
                //construct a vector of indices of unique station depths
                size_t nlevels = ConstructDepthIndices(MeasDepthIndices, ShiftDepth,
                    Model, Data.GetMeasPosZ());
                //we need to call the previous function and take the retun type
                //but we want to suppress compiler warnings about unusued variables
                ignore(nlevels);
                std::complex<double> Tx, Ty;
                double Frequency = Data.GetFrequencies().at(calcindex);
                std::vector<std::complex<double>> Hx1(FieldCalculator->GetHx1(Frequency)),
                    Hx2(FieldCalculator->GetHx2(Frequency));
                std::vector<std::complex<double>> Hy1(FieldCalculator->GetHy1(Frequency)),
                    Hy2(FieldCalculator->GetHy2(Frequency));
                std::vector<std::complex<double>> Hz1(FieldCalculator->GetHz1(Frequency)),
                    Hz2(FieldCalculator->GetHz2(Frequency));
                for (size_t j = 0; j < nstats; ++j)
                  {

                    boost::array<ThreeDModelBase::t3DModelData::index, 3> StationHxIndex =
                        Model.FindAssociatedIndices(
                            Data.GetMeasPosX()[Data.GetHxIndices()[j + ind_shift]],
                            Data.GetMeasPosY()[Data.GetHxIndices()[j + ind_shift]],
                            Data.GetMeasPosZ()[Data.GetHxIndices()[j + ind_shift]]);
                    const size_t offset_Hx = (nmodx * nmody)
                        * MeasDepthIndices[Data.GetHxIndices()[j + ind_shift]]
                        + StationHxIndex[0] * nmody + StationHxIndex[1];

                    boost::array<ThreeDModelBase::t3DModelData::index, 3> StationHyIndex =
                        Model.FindAssociatedIndices(
                            Data.GetMeasPosX()[Data.GetHyIndices()[j + ind_shift]],
                            Data.GetMeasPosY()[Data.GetHyIndices()[j + ind_shift]],
                            Data.GetMeasPosZ()[Data.GetHyIndices()[j + ind_shift]]);
                    const size_t offset_Hy = (nmodx * nmody)
                        * MeasDepthIndices[Data.GetHyIndices()[j + ind_shift]]
                        + StationHyIndex[0] * nmody + StationHyIndex[1];

                    boost::array<ThreeDModelBase::t3DModelData::index, 3> StationHzIndex =
                        Model.FindAssociatedIndices(
                            Data.GetMeasPosX()[Data.GetHzIndices()[j + ind_shift]],
                            Data.GetMeasPosY()[Data.GetHzIndices()[j + ind_shift]],
                            Data.GetMeasPosZ()[Data.GetHzIndices()[j + ind_shift]]);
                    const size_t offset_Hz = (nmodx * nmody)
                        * MeasDepthIndices[Data.GetHzIndices()[j + ind_shift]]
                        + StationHzIndex[0] * nmody + StationHzIndex[1];

                    FieldsToTipper(Hx1[offset_Hx], Hx2[offset_Hx], Hy1[offset_Hy],
                        Hy2[offset_Hy], Hz1[offset_Hz], Hz2[offset_Hz], Tx, Ty);
                    //result is a local array for this frequency
                    //so we can directly use it even in a threaded environment
                    const size_t meas_index = j * 4;

                    result[startindex + meas_index] = Tx.real();
                    result[startindex + meas_index + 1] = Tx.imag();
                    result[startindex + meas_index + 2] = Ty.real();
                    result[startindex + meas_index + 3] = Ty.imag();

                  }
#ifdef HAVEOPENMP
                omp_set_lock(&lck);
#endif
                std::chrono::system_clock::time_point end =
                    std::chrono::system_clock::now();
                size_t duration = std::chrono::duration_cast<std::chrono::seconds>(
                    end - start).count();

                NewExecTime.push_back(std::make_pair(duration, calcindex));
#ifdef HAVEOPENMP
                omp_unset_lock(&lck);
#endif
              } catch (jif3D::FatalException &e)
              {
                ErrorMsg = e.what();
                FatalError = true;
              } catch (...)
              {
                FatalError = true;
                std::cerr << "Problem in MT forward calculation.";
              }

            //finished with one frequency
          }

        std::sort(NewExecTime.begin(), NewExecTime.end());
        ForwardExecTime = NewExecTime;
#ifdef HAVEOPENMP
        omp_destroy_lock(&lck);
#endif
        //we cannot throw from within the openmp section so if there was an exception
        //inside the parallel region we set FatalErrror to true and throw a new exception here
        if (FatalError)
          throw jif3D::FatalException("Problem in MT forward calculation: " + ErrorMsg,
          __FILE__,
          __LINE__);

        return result;

      }

    rvec X3DTipperCalculator::LQDerivative(const ModelType &Model, const TipperData &Data,
        const rvec &Misfit, size_t minfreqindex, size_t maxfreqindex)
      {
        //we define nfreq as int to make the compiler happy in the openmp loop
        maxfreqindex = std::min(maxfreqindex, Data.GetFrequencies().size());
        const int nfreq = maxfreqindex - minfreqindex;
        if (DerivExecTime.empty())
          {
            for (int i = 0; i < nfreq; ++i)
              {
                DerivExecTime.push_back(std::make_pair(0, minfreqindex + i));
              }
          }
        std::string ErrorMsg;
        //a few commonly used quantities for shorter notation

        const size_t nstats = Data.GetHxIndices().size() / nfreq;
        const size_t nmod = Model.GetNModelElements();
        assert(Misfit.size() == nstats * nfreq * 4);

        jif3D::rvec Gradient(nmod);

        bool FatalError = false;
        //we need to initialize all values to zero as we are adding
        //the individual gradients per frequency
        std::fill(Gradient.begin(), Gradient.end(), 0.0);

        std::vector<std::pair<size_t, size_t>> NewExecTime;
#ifdef HAVEOPENMP
        omp_lock_t lck;
        omp_init_lock(&lck);
#endif
        //we parallelize the gradient calculation by frequency
        //see also the comments for the forward calculation
        //here the explicitly shared variable is Gradient
        //all others are predetermined to be shared
#pragma omp parallel for shared(Gradient, NewExecTime) schedule(dynamic,1)
        for (int i = 0; i < nfreq; ++i)
          {
            try
              {
                //we want to alternate between items at the beginning of the map and at the end of the map
                const size_t queueindex = (i % 2) == 0 ? i / 2 : nfreq - 1 - i / 2;
                const size_t calcindex = DerivExecTime.at(queueindex).second;
                std::vector<double> C;
                ForwardInfo Info(Model, C, calcindex, TempDir.string(), X3DName, NameRoot,
                    GreenType1, GreenType4);
                //calculate the gradient for each frequency
                std::chrono::system_clock::time_point start =
                    std::chrono::system_clock::now();
                GradResult tmp = TipperDerivativeFreq(Info, Data, Misfit,
                    *FieldCalculator);
                std::chrono::system_clock::time_point end =
                    std::chrono::system_clock::now();
                size_t duration = std::chrono::duration_cast<std::chrono::seconds>(
                    end - start).count();
#ifdef HAVEOPENMP
                omp_set_lock(&lck);
#endif

                NewExecTime.push_back(std::make_pair(duration, calcindex));
                //the total gradient is the sum over the gradients for each frequency
                std::transform(tmp.Gradient.begin(), tmp.Gradient.end(), Gradient.begin(),
                    Gradient.begin(), std::plus<double>());
#ifdef HAVEOPENMP
                omp_unset_lock(&lck);
#endif
              } catch (jif3D::FatalException &e)
              {
                ErrorMsg = e.what();
                FatalError = true;
              } catch (...)
              {
                //we cannot throw exceptions that leave the parallel region
                //so we catch everything and set FatalError to true
                //then outside the parallel region we throw a new exception
                FatalError = true;
              }
            //finished with one frequency
          }
        std::sort(NewExecTime.begin(), NewExecTime.end());
        DerivExecTime = NewExecTime;
#ifdef HAVEOPENMP
        omp_destroy_lock(&lck);
#endif
        //if we had some exception inside the openmp region, we throw
        // a generic error message
        if (FatalError)
          throw jif3D::FatalException(
              "Problem in Tipper gradient calculation. Error message: " + ErrorMsg,
              __FILE__,
              __LINE__);

        return 2.0 * Gradient;

      }

  } /* namespace jif3D */
