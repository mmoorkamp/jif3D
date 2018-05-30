/*
 * X3DTipperCalculator.cpp
 *
 *  Created on: 25 May 2018
 *      Author: mm489
 */

#include <complex>
#include <map>
#include <chrono>
#ifdef HAVEOPENMP
#include <omp.h>
#endif

#include "X3DTipperCalculator.h"
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

    X3DTipperCalculator::X3DTipperCalculator(boost::filesystem::path TDir,
        std::string x3d, bool Clean , std::vector<boost::shared_ptr<jif3D::X3DFieldCalculator> > FC) :
        GreenType1(hst), GreenType4(hst), X3DName(x3d), CleanFiles(Clean), FieldCalculators(FC)
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

    X3DTipperCalculator::~X3DTipperCalculator()
      {

      }

    rvec X3DTipperCalculator::Calculate(ModelType &Model, size_t minfreqindex,
        size_t maxfreqindex)
      {
        //we define nfreq as int to make the compiler happy in the openmp loop
        assert(minfreqindex <= maxfreqindex);
        maxfreqindex = std::min(maxfreqindex, Model.GetFrequencies().size());

        const int nfreq = maxfreqindex - minfreqindex;
        if (FieldCalculators.empty() || FieldCalculators.size() != nfreq)
          {
            FieldCalculators.resize(nfreq);
            for (auto &fc : FieldCalculators)
              {
                fc = boost::make_shared<jif3D::X3DFieldCalculator>(TempDir, X3DName);
              }
          }
        const size_t nmeas = Model.GetMeasPosX().size();
        if (ForwardExecTime.empty() || ForwardExecTime.size() != nfreq)
          {
            for (int i = 0; i < nfreq; ++i)
              {
                ForwardExecTime.push_back(std::make_pair(0, minfreqindex + i));
              }
          }
        //if the current model does not contain any ExIndices information
        //generate ExIndices, EyIndices and HIndices as 0:nmeas for each Frequency
        //here we assume that we either have all three indices in the netCDF file or none of them
        std::vector<int> ExIndices(Model.GetExIndices()), EyIndices(Model.GetEyIndices()),
            HIndices(Model.GetHIndices());
        size_t ishift = 0;
        if (ExIndices.empty())
          {
            ExIndices.resize(nmeas * nfreq);
            EyIndices.resize(nmeas * nfreq);
            HIndices.resize(nmeas * nfreq);
            for (int ifr = 0; ifr < nfreq; ++ifr)
              {
                ishift = nmeas * ifr;
                for (size_t i = 0; i < nmeas; ++i)
                  {
                    ExIndices[i + ishift] = i;
                  }
              }
            EyIndices = ExIndices;
            HIndices = ExIndices;
          }
        Model.SetFieldIndices(ExIndices, EyIndices, HIndices);

        const size_t nstats = Model.GetExIndices().size() / nfreq;
        const size_t nmodx = Model.GetXCoordinates().size();
        const size_t nmody = Model.GetYCoordinates().size();

        std::string ErrorMsg;
        //result will hold the final impedance values with
        //applied distortion correction
        jif3D::rvec result(nstats * nfreq * 4);

        bool FatalError = false;
        result.clear();

        //we make a call to the coordinate functions to make sure
        //that we have updated the coordinate information and cached it
        //only then the subsequent calls are thread safe
        Model.GetXCoordinates();
        Model.GetYCoordinates();
        Model.GetZCoordinates();

        //check that the depths to the different background layers match
        //with the depths to grid cell boundaries
        std::vector<double> BGDepths(Model.GetBackgroundThicknesses().size(), 0.0);
        std::partial_sum(Model.GetBackgroundThicknesses().begin(),
            Model.GetBackgroundThicknesses().end(), BGDepths.begin());
        CompareDepths(BGDepths, Model.GetZCoordinates());
        std::vector<std::pair<size_t, size_t>> NewExecTime;

        omp_lock_t lck;
        omp_init_lock(&lck);

        //we parallelize by frequency, this is relatively simple
        //but we can use up to 20 processors for typical MT problems
        // as we do not have the source for x3d, this is our only possibility anyway
        //the const qualified variables above are predetermined to be shared by the openmp standard
#pragma omp parallel for shared(result,NewExecTime) schedule(dynamic,1)
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

                const size_t currindex = calcindex - minfreqindex;
                const size_t startindex = nstats * currindex * 4;
                const size_t ind_shift = nstats * calcindex;

                std::chrono::system_clock::time_point start =
                    std::chrono::system_clock::now();

                std::vector<double> ShiftDepth;
                std::vector<size_t> MeasDepthIndices;
                //construct a vector of indices of unique station depths
                size_t nlevels = ConstructDepthIndices(MeasDepthIndices, ShiftDepth,
                    Model);

                FieldCalculators.at(calcindex)->CalculateFields(Model, calcindex);
                std::complex<double> Tx, Ty;

                for (size_t j = 0; j < nstats; ++j)
                  {

                    boost::array<ThreeDModelBase::t3DModelData::index, 3> StationHIndex =
                        Model.FindAssociatedIndices(
                            Model.GetMeasPosX()[Model.GetHIndices()[j + ind_shift]],
                            Model.GetMeasPosY()[Model.GetHIndices()[j + ind_shift]],
                            Model.GetMeasPosZ()[Model.GetHIndices()[j + ind_shift]]);

                    const size_t offset_H = (nmodx * nmody)
                        * MeasDepthIndices[Model.GetHIndices()[j + ind_shift]]
                        + StationHIndex[0] * nmody + StationHIndex[1];

                    FieldsToTipper(FieldCalculators.at(calcindex)->GetHx1()[offset_H],
                        FieldCalculators.at(calcindex)->GetHx2()[offset_H],
                        FieldCalculators.at(calcindex)->GetHy1()[offset_H],
                        FieldCalculators.at(calcindex)->GetHy2()[offset_H],
                        FieldCalculators.at(calcindex)->GetHz1()[offset_H],
                        FieldCalculators.at(calcindex)->GetHz2()[offset_H], Tx, Ty);
                    //result is a local array for this frequency
                    //so we can directly use it even in a threaded environment
                    const size_t meas_index = j * 4;

                    result[startindex + meas_index] = Tx.real();
                    result[startindex + meas_index + 1] = Tx.imag();
                    result[startindex + meas_index + 2] = Ty.real();
                    result[startindex + meas_index + 3] = Ty.imag();

                  }

                omp_set_lock(&lck);
                std::chrono::system_clock::time_point end =
                    std::chrono::system_clock::now();
                size_t duration = std::chrono::duration_cast<std::chrono::seconds>(
                    end - start).count();

                NewExecTime.push_back(std::make_pair(duration, calcindex));
                omp_unset_lock(&lck);
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

        omp_destroy_lock(&lck);

        //we cannot throw from within the openmp section so if there was an exception
        //inside the parallel region we set FatalErrror to true and throw a new exception here
        if (FatalError)
          throw jif3D::FatalException("Problem in MT forward calculation: " + ErrorMsg,
          __FILE__,
          __LINE__);

        return result;

      }

    rvec X3DTipperCalculator::LQDerivative(const ModelType &Model, const rvec &Misfit,
        size_t minfreqindex, size_t maxfreqindex)
      {
        //we define nfreq as int to make the compiler happy in the openmp loop
        maxfreqindex = std::min(maxfreqindex, Model.GetFrequencies().size());
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
        const size_t nmodx = Model.GetConductivities().shape()[0];
        const size_t nmody = Model.GetConductivities().shape()[1];
        const size_t nmodz = Model.GetConductivities().shape()[2];
        const size_t nstats = Model.GetExIndices().size() / nfreq;
        const size_t nmod = nmodx * nmody * nmodz;
        assert(Misfit.size() == nstats * nfreq * 4);

        jif3D::rvec Gradient(nmod);

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

        std::vector<std::pair<size_t, size_t>> NewExecTime;
        omp_lock_t lck;
        omp_init_lock(&lck);

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
                GradResult tmp = TipperDerivativeFreq(Info, Misfit,
                    FieldCalculators.at(calcindex));
                std::chrono::system_clock::time_point end =
                    std::chrono::system_clock::now();
                size_t duration = std::chrono::duration_cast<std::chrono::seconds>(
                    end - start).count();

                omp_set_lock(&lck);
                NewExecTime.push_back(std::make_pair(duration, calcindex));
                //the total gradient is the sum over the gradients for each frequency
                std::transform(tmp.Gradient.begin(), tmp.Gradient.end(), Gradient.begin(),
                    Gradient.begin(), std::plus<double>());
                omp_unset_lock(&lck);
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
        omp_destroy_lock(&lck);

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
