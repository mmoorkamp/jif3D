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
#include <map>
#ifdef HAVEOPENMP
#include <omp.h>
#endif
#include <boost/assign/list_of.hpp> // for 'map_list_of()'
#include <boost/multi_array.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/math/constants/constants.hpp>
#include "../Global/FatalException.h"
#include "../Global/convert.h"
#include "../ModelBase/CellBoundaries.h"
#include "X3DMTCalculator.h"
#include "ReadWriteX3D.h"
#include "MTEquations.h"
#include "InterpolateField.h"

namespace fs = boost::filesystem;

namespace jif3D
  {
    //define some names that are always the same
    //either because x3d uses this convention
    //or because we use them in their own directory
    //and want to keep them simple to make sure x3d can handle them
    const std::string modelfilename("x3d.model");
    const std::string resultfilename("x3d.result");
    const std::string emaname = resultfilename + "0.ema";
    const std::string sourceafilename = modelfilename + "0a.source";
    const std::string sourcebfilename = modelfilename + "0b.source";
    const std::string emoAname = resultfilename + "0a.emo";
    const std::string emoBname = resultfilename + "0b.emo";
    const std::string emaAname = resultfilename + "0a.ema";
    const std::string emaBname = resultfilename + "0b.ema";
    const std::string runext = "_run";
    const std::string dirext = "_dir";
    //associate a type of calculation with a name string
    const std::map<X3DModel::ProblemType, std::string> Extension =
        boost::assign::map_list_of(X3DModel::MT, "MT")(X3DModel::EDIP, "EDIP")(
            X3DModel::MDIP, "MDIP");

    inline void CheckField(const std::vector<std::complex<double> > &Field, size_t nelem)
      {
        if (Field.size() != nelem)
          throw jif3D::FatalException(
              "Number of read in elements in Field: " + jif3D::stringify(Field.size())
                  + " does not match expected: " + jif3D::stringify(nelem));
      }

    //check that the .hnk file for x3d are in a certain directory
    bool CheckHNK(const fs::path &TargetDir)
      {
        return fs::exists(TargetDir / "ndec15.hnk")
            && fs::exists(TargetDir / "ndec20.hnk")
            && fs::exists(TargetDir / "ndec30.hnk")
            && fs::exists(TargetDir / "ndec40.hnk");
      }

    void X3DMTCalculator::CleanUp()
      {
        //under certain conditions we might not be able
        //to delete all files. We don't want the program to stop
        //because of this as we can always delete files afterwards
        //so we pass an error code object to remove_all and ignore the error
        boost::system::error_code ec;
        fs::directory_iterator end_itr; // default construction yields past-the-end
        //go through the directory and delete any file that starts with NameRoot
        for (fs::directory_iterator itr(fs::current_path()); itr != end_itr; ++itr)
          {
            if (boost::algorithm::starts_with(itr->path().filename().string(), NameRoot))
              {
                fs::remove_all(itr->path().filename(), ec);
              }
          }
      }

    X3DMTCalculator::X3DMTCalculator(boost::filesystem::path TDir, bool DC) :
        WantDistCorr(DC)
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

    //copy the .hnk files for x3d from SourceDir to TargetDir
    void CopyHNK(const fs::path &SourceDir, const fs::path &TargetDir)
      {
        //copy file fails with an exception if the target exists
        //so we check before we do the actual copy
        if (!CheckHNK(TargetDir))
          {
            fs::copy_file(SourceDir / "ndec15.hnk", TargetDir / "ndec15.hnk");
            fs::copy_file(SourceDir / "ndec20.hnk", TargetDir / "ndec20.hnk");
            fs::copy_file(SourceDir / "ndec30.hnk", TargetDir / "ndec30.hnk");
            fs::copy_file(SourceDir / "ndec40.hnk", TargetDir / "ndec40.hnk");
          }
      }

    //create a unique ID that we can use to name things and still
    //perform parallel calculations
    std::string X3DMTCalculator::ObjectID()
      {
        //a unique ID created on construction
        boost::uuids::uuid tag = boost::uuids::random_generator()();
        //make a unique filename for the sensitivity file created by this object
        //we use boost uuid to generate a unique identifier tag
        //and translate it to a string to generate the filename
        return "mt" + jif3D::stringify(getpid()) + jif3D::stringify(tag);
      }

    std::string X3DMTCalculator::MakeUniqueName(X3DModel::ProblemType Type,
        const size_t FreqIndex)
      {
        //we assemble the name from the id of the process
        //and the address of the current object
        std::string result(NameRoot);
        //the type of calculation
        result += Extension.find(Type)->second;
        //and the frequency index
        result += jif3D::stringify(FreqIndex);
        return result;
      }
    //create a script that changes to the correct directory
    //and executes x3d in that directory
    void MakeRunFile(const std::string &NameRoot, const std::string DirName)
      {
        std::string RunFileName = NameRoot + runext;
        fs::create_directory(DirName);
        std::ofstream runfile;
        runfile.open(RunFileName.c_str());
        runfile << "#!/bin/bash\n";
        runfile << "cd " << DirName << "\n";
        runfile << "x3d > /dev/null\n";
        runfile << "cd ..\n";
        runfile.close();
        //we also copy the necessary *.hnk files
        //from the current directory to the work directory
        CopyHNK(fs::current_path(), DirName);
      }

    //execute the script that runs x3d
    void RunX3D(const std::string &NameRoot)
      {
        //instead of making the script executable
        //we run a bash with the scriptname as an argument
        //this turns out to be more robust
        const std::string runname = "bash " + NameRoot + runext;
        //it is important to include the std:: namespace specification
        //for the system call, otherwise the GNU compiler picks up
        //a version from the c library that gives trouble in threaded environments
        if (std::system(runname.c_str()))
          throw FatalException("Cannot execute run script: " + runname);
      }

    rvec X3DMTCalculator::CalculateFrequency(const X3DModel &Model,
        const std::vector<double> &C, size_t freqindex)
      {
        const size_t nmeas = Model.GetMeasPosX().size();
        const size_t nmodx = Model.GetXCoordinates().size();
        const size_t nmody = Model.GetYCoordinates().size();

        jif3D::rvec result(nmeas * 8);
        fs::path RootName = TempDir / MakeUniqueName(X3DModel::MT, freqindex);
        fs::path DirName = RootName.string() + dirext;
        std::vector<double> CurrFreq(1, Model.GetFrequencies()[freqindex]);
        std::vector<double> ShiftDepth;
        std::vector<size_t> MeasDepthIndices;
        //construct a vector of indices of unique station depths
        size_t nlevels = ConstructDepthIndices(MeasDepthIndices, ShiftDepth, Model);
        //writing out files causes problems in parallel
        // so we make sure it is done one at a time
#pragma omp critical(forward_write_files)
          {
            MakeRunFile(RootName.string(), DirName.string());
            WriteProjectFile(DirName.string(), CurrFreq, X3DModel::MT, resultfilename,
                modelfilename);
            Write3DModelForX3D((DirName / modelfilename).string(), Model.GetXCellSizes(),
                Model.GetYCellSizes(), Model.GetZCellSizes(), ShiftDepth,
                Model.GetConductivities(), Model.GetBackgroundConductivities(),
                Model.GetBackgroundThicknesses());
          }
        //run x3d in parallel
        RunX3D(RootName.string());
        std::vector<std::complex<double> > Ex1, Ex2, Ey1, Ey2, Hx1, Hx2, Hy1, Hy2;
        std::complex<double> Zxx, Zxy, Zyx, Zyy;
        //read in the electric and magnetic field at the observe sites
#pragma omp critical(forward_read_emo)
          {
            ReadEMO((DirName / emoAname).string(), Ex1, Ey1, Hx1, Hy1);
            ReadEMO((DirName / emoBname).string(), Ex2, Ey2, Hx2, Hy2);
          }
        const size_t nval = (nmodx * nmody * nlevels);
        CheckField(Ex1, nval);
        CheckField(Ex2, nval);
        CheckField(Ey1, nval);
        CheckField(Ey2, nval);
        CheckField(Hx1, nval);
        CheckField(Hx2, nval);
        CheckField(Hy1, nval);
        CheckField(Hy2, nval);
        //calculate impedances from the field spectra for all measurement sites
        for (size_t j = 0; j < nmeas; ++j)
          {
            boost::array<ThreeDModelBase::t3DModelData::index, 3> StationIndex =
                Model.FindAssociatedIndices(Model.GetMeasPosX()[j],
                    Model.GetMeasPosY()[j], Model.GetMeasPosZ()[j]);
// with the current equations we cannot use interpolation as the position
            // of the source in the adjoint calculation is always smeared
            //across the whole cell, this works best for a cell in the centre
            //of the cell and any other position deteriorates convergence
//            std::complex<double> Ex1Inter = InterpolateField(Ex1, Model, j,
//                MeasDepthIndices);
//            std::complex<double> Ex2Inter = InterpolateField(Ex2, Model, j,
//                MeasDepthIndices);
//            std::complex<double> Ey1Inter = InterpolateField(Ey1, Model, j,
//                MeasDepthIndices);
//            std::complex<double> Ey2Inter = InterpolateField(Ey2, Model, j,
//                MeasDepthIndices);
//            std::complex<double> Hx1Inter = InterpolateField(Hx1, Model, j,
//                MeasDepthIndices);
//            std::complex<double> Hx2Inter = InterpolateField(Hx2, Model, j,
//                MeasDepthIndices);
//            std::complex<double> Hy1Inter = InterpolateField(Hy1, Model, j,
//                MeasDepthIndices);
//            std::complex<double> Hy2Inter = InterpolateField(Hy2, Model, j,
//                MeasDepthIndices);
//            FieldsToImpedance(Ex1Inter, Ex2Inter, Ey1Inter, Ey2Inter, Hx1Inter, Hx2Inter,
//                Hy1Inter, Hy2Inter, Zxx, Zxy, Zyx, Zyy);
            const size_t offset = (nmodx * nmody) * MeasDepthIndices[j]
                + StationIndex[0] * nmody + StationIndex[1];
            FieldsToImpedance(Ex1[offset], Ex2[offset], Ey1[offset], Ey2[offset],
                Hx1[offset], Hx2[offset], Hy1[offset], Hy2[offset], Zxx, Zxy, Zyx, Zyy);
            //result is a local array for this frequency
            //so we can directly use it even in a threaded environment
            const size_t meas_index = j * 8;
            const size_t site_index = j * 4;
            const size_t imp_index = freqindex * nmeas * 8;
            result(meas_index) = C[site_index] * Zxx.real()
                + C[site_index + 1] * Zyx.real();
            result(meas_index + 1) = C[site_index] * Zxx.imag()
                + C[site_index + 1] * Zyx.imag();
            result(meas_index + 2) = C[site_index] * Zxy.real()
                + C[site_index + 1] * Zyy.real();
            result(meas_index + 3) = C[site_index] * Zxy.imag()
                + C[site_index + 1] * Zyy.imag();
            result(meas_index + 4) = C[site_index + 3] * Zyx.real()
                + C[site_index + 2] * Zxx.real();
            result(meas_index + 5) = C[site_index + 3] * Zyx.imag()
                + C[site_index + 2] * Zxx.imag();
            result(meas_index + 6) = C[site_index + 3] * Zyy.real()
                + C[site_index + 2] * Zxy.real();
            result(meas_index + 7) = C[site_index + 3] * Zyy.imag()
                + C[site_index + 2] * Zxy.real();
            RawImpedance(imp_index + meas_index) = Zxx.real();
            RawImpedance(imp_index + meas_index + 1) = Zxx.imag();
            RawImpedance(imp_index + meas_index + 2) = Zxy.real();
            RawImpedance(imp_index + meas_index + 3) = Zxy.imag();
            RawImpedance(imp_index + meas_index + 4) = Zyx.real();
            RawImpedance(imp_index + meas_index + 5) = Zyx.imag();
            RawImpedance(imp_index + meas_index + 6) = Zyy.real();
            RawImpedance(imp_index + meas_index + 7) = Zyy.imag();

          }
        return result;
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

    rvec X3DMTCalculator::Calculate(const X3DModel &Model, size_t minfreqindex,
        size_t maxfreqindex)
      {
        //we define nfreq as int to make the compiler happy in the openmp loop
        assert(minfreqindex <= maxfreqindex);
        maxfreqindex = std::min(maxfreqindex, Model.GetFrequencies().size());
        const size_t nfreq = maxfreqindex - minfreqindex;
        const size_t nmeas = Model.GetMeasPosX().size();
        const size_t nmodx = Model.GetXCoordinates().size();
        const size_t nmody = Model.GetYCoordinates().size();

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
#endif
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
                rvec freqresult = CalculateFrequency(Model, C, i);
                size_t startindex = nmeas * i * 8;
#ifdef HAVEOPENMP
                omp_set_lock(&lck);
#endif
                std::copy(freqresult.begin(), freqresult.end(),
                    result.begin() + startindex);
#ifdef HAVEOPENMP
                omp_unset_lock(&lck);
#endif
              } catch (...)
              {
                FatalError = true;
                std::cerr << "Problem in MT forward calculation.";
              }
            //finished with one frequency
          }
#ifdef HAVEOPENMP
        omp_destroy_lock(&lck);
#endif
        //we cannot throw from within the openmp section so if there was an exception
        //inside the parallel region we set FatalErrror to true and throw a new exception here
        if (FatalError)
          throw jif3D::FatalException("Problem in MT forward calculation.");
        return result;

      }

    cmat CalcEExt(const rvec &Misfit, const std::vector<double> &C,
        const size_t startindex, const size_t freq_start_index,
        const std::complex<double> &Hx1, const std::complex<double> &Hx2,
        const std::complex<double> &Hy1, const std::complex<double> &Hy2)
      {
        const size_t siteindex = freq_start_index + startindex * 8;
        cmat AH(2, 2);
        cmat CtAH(2, 2);
        const std::complex<double> magdet = 1. / (Hx1 * Hy2 - Hx2 * Hy1);
        const std::complex<double> A00(Misfit(siteindex), Misfit(siteindex + 1));
        const std::complex<double> A01(Misfit(siteindex + 2), Misfit(siteindex + 3));
        const std::complex<double> A10(Misfit(siteindex + 4), Misfit(siteindex + 5));
        const std::complex<double> A11(Misfit(siteindex + 6), Misfit(siteindex + 7));
        AH(0, 0) = magdet * (conj(A00) * Hy2 - conj(A01) * Hx2);
        AH(0, 1) = magdet * (-conj(A00) * Hy1 + conj(A01) * Hx1);
        AH(1, 0) = magdet * (conj(A10) * Hy2 - conj(A11) * Hx2);
        AH(1, 1) = magdet * (-conj(A10) * Hy1 + conj(A11) * Hx1);
        CtAH(0, 0) = AH(0, 0) * C[startindex * 4] + AH(1, 0) * C[startindex * 4 + 2];
        CtAH(0, 1) = AH(0, 1) * C[startindex * 4] + AH(1, 1) * C[startindex * 4 + 2];
        CtAH(1, 0) = AH(0, 0) * C[startindex * 4 + 1] + AH(1, 0) * C[startindex * 4 + 3];
        CtAH(1, 1) = AH(0, 1) * C[startindex * 4 + 1] + AH(1, 1) * C[startindex * 4 + 3];
        return CtAH;
      }

    void CalcHext(const std::complex<double> &omega_mu, std::complex<double> &Xp1,
        std::complex<double> &Xp2, std::complex<double> &Yp1, std::complex<double> &Yp2,
        const std::complex<double> &Zxx, const std::complex<double> &Zxy,
        const std::complex<double> &Zyx, const std::complex<double> &Zyy)
      {
        //we need temporary variables as we calculate a new  value for Xp1 in the first line
        //but need the old value for Xp1 in the second line
        std::complex<double> temp1 = conj(Zxx) * Xp1 + conj(Zyx) * Yp1;
        std::complex<double> temp2 = conj(Zxy) * Xp1 + conj(Zyy) * Yp1;
        Xp1 = temp1 * omega_mu;
        Yp1 = temp2 * omega_mu;
        //the same remark applies to Xp2
        temp1 = conj(Zxx) * Xp2 + conj(Zyx) * Yp2;
        temp2 = conj(Zxy) * Xp2 + conj(Zyy) * Yp2;
        Xp2 = temp1 * omega_mu;
        Yp2 = temp2 * omega_mu;
      }

    void CalcU(const std::string &RootName,
        const std::vector<std::complex<double> > &XPolMoments,
        const std::vector<std::complex<double> > &YPolMoments,
        std::vector<std::complex<double> > &Ux, std::vector<std::complex<double> > &Uy,
        std::vector<std::complex<double> > &Uz, const std::vector<size_t> &SourceXIndex,
        const std::vector<size_t> &SourceYIndex,
        const std::vector<double> ObservationDepths,
        const jif3D::ThreeDModelBase::t3DModelDim &ZCellBoundaries,
        const jif3D::ThreeDModelBase::t3DModelDim &ZCellSizes, const size_t ncellsx,
        const size_t ncellsy, const size_t ncellsz)
      {
        std::string DirName = RootName + dirext + "/";
#pragma omp critical(calcU_writesource)
          {
            WriteSourceFile(DirName + sourceafilename, SourceXIndex, SourceYIndex,
                ObservationDepths, XPolMoments, YPolMoments, ZCellBoundaries, ZCellSizes,
                ncellsx, ncellsy);
          }
        RunX3D(RootName);
#pragma omp critical(calcU_readema)
          {
            ReadEMA(DirName + emaname, Ux, Uy, Uz, ncellsx, ncellsy, ncellsz);
          }
        //boost::filesystem::remove_all(emaname);
      }

    rvec X3DMTCalculator::LQDerivativeFreq(const X3DModel &Model, const rvec &Misfit,
        const std::vector<double> &C, size_t freqindex)
      {
        //a few commonly used quantities for shorter notation
        const size_t nmodx = Model.GetConductivities().shape()[0];
        const size_t nmody = Model.GetConductivities().shape()[1];
        const size_t nmodz = Model.GetConductivities().shape()[2];
        //the number of observations in the model file, one for each cell in the layer
        const size_t nobs = nmodx * nmody;
        //the number of measurement sites
        const size_t nmeas = Model.GetMeasPosX().size();
        const size_t nmod = nmodx * nmody * nmodz;
        std::vector<double> ShiftDepth;
        //for the controlled source calculations we do not actually
        //need any observe layers as we are only interested in the
        //anomalous fields
        std::vector<double> SourceObserve(1, 0.0);
        std::vector<size_t> MeasDepthIndices;
        size_t nlevels = ConstructDepthIndices(MeasDepthIndices, ShiftDepth, Model);
        jif3D::rvec Gradient(nmod, 0.0);

        fs::path ForwardDirName = TempDir
            / (MakeUniqueName(X3DModel::MT, freqindex) + dirext);
        if (!fs::is_directory(ForwardDirName))
          throw FatalException(
              "In X3D gradient calculation, directory does not exist: "
                  + ForwardDirName.string());

        //read the fields from the forward calculation
        std::vector<std::complex<double> > Ex1_obs, Ex2_obs, Ey1_obs, Ey2_obs, Hx1_obs,
            Hx2_obs, Hy1_obs, Hy2_obs;
#pragma omp critical(gradient_reademo)
          {
            ReadEMO((ForwardDirName / emoAname).string(), Ex1_obs, Ey1_obs, Hx1_obs,
                Hy1_obs);
            ReadEMO((ForwardDirName / emoBname).string(), Ex2_obs, Ey2_obs, Hx2_obs,
                Hy2_obs);
          }
        const size_t nfield = nobs * nlevels;
        CheckField(Ex1_obs, nfield);
        CheckField(Ex2_obs, nfield);
        CheckField(Ey1_obs, nfield);
        CheckField(Ey2_obs, nfield);
        CheckField(Hx1_obs, nfield);
        CheckField(Hx2_obs, nfield);
        CheckField(Hy1_obs, nfield);
        CheckField(Hy2_obs, nfield);

        std::vector<std::complex<double> > Ex1_all, Ex2_all, Ey1_all, Ey2_all, Ez1_all,
            Ez2_all;
        //for the gradient calculation we also need the electric fields
        //at all cells in the model for the two source polarizations of
        //the forward calculations
#pragma omp critical(gradient_readema)
          {
            ReadEMA((ForwardDirName / emaAname).string(), Ex1_all, Ey1_all, Ez1_all,
                nmodx, nmody, nmodz);
            ReadEMA((ForwardDirName / emaBname).string(), Ex2_all, Ey2_all, Ez2_all,
                nmodx, nmody, nmodz);
          }
        //create variables for the adjoint field calculation
        const size_t freq_start_index = nmeas * freqindex * 8;
        std::vector<std::complex<double> > XPolMoments1(nmeas), XPolMoments2(nmeas),
            YPolMoments1(nmeas), YPolMoments2(nmeas), Zeros(nmeas);
        std::vector<size_t> SourceXIndex(nmeas), SourceYIndex(nmeas);
        //we only calculate sources for the observe sites, so
        //we make sure everything else is zero
        std::fill(Zeros.begin(), Zeros.end(), 0.0);

        //make the sources for the electric dipoles
        for (size_t j = 0; j < nmeas; ++j)
          {
            boost::array<ThreeDModelBase::t3DModelData::index, 3> StationIndex =
                Model.FindAssociatedIndices(Model.GetMeasPosX()[j],
                    Model.GetMeasPosY()[j], Model.GetMeasPosZ()[j]);

            const size_t offset = (nmodx * nmody) * MeasDepthIndices[j]
                + StationIndex[0] * nmody + StationIndex[1];
            SourceXIndex.at(j) = StationIndex[0];
            SourceYIndex.at(j) = StationIndex[1];

            //this is an implementation of eq. 12 in Avdeev and Avdeeva
            //we do not have any beta, as this is part of the misfit
//            cmat j_ext = CalcATimesH(MisfitToA(Misfit, siteindex),
//                MakeH(Hx1_obs[offset], Hx2_obs[offset], Hy1_obs[offset], Hy2_obs[offset]));
            cmat j_ext = CalcEExt(Misfit, C, j, freq_start_index, Hx1_obs[offset],
                Hx2_obs[offset], Hy1_obs[offset], Hy2_obs[offset]);
            XPolMoments1.at(j) = conj(j_ext(0, 0));
            YPolMoments1.at(j) = conj(j_ext(1, 0));
            XPolMoments2.at(j) = conj(j_ext(0, 1));
            YPolMoments2.at(j) = conj(j_ext(1, 1));
          }
        //we only want to calculate for one frequency
        //so our vector has just 1 element
        std::vector<double> CurrFreq(1, 0.0);
        CurrFreq[0] = Model.GetFrequencies()[freqindex];
        fs::path EdipName = TempDir / MakeUniqueName(X3DModel::EDIP, freqindex);
        fs::path EdipDirName = EdipName.string() + dirext;
        //again we have to write out some file for the electric
        //dipole calculation with x3d, this shouldn't be done
        //in parallel
#pragma omp critical(gradient_writemodel_edip)
          {
            MakeRunFile(EdipName.string(), EdipDirName.string());
            WriteProjectFile(EdipDirName.string(), CurrFreq, X3DModel::EDIP,
                resultfilename, modelfilename);
            Write3DModelForX3D((EdipDirName / modelfilename).string(),
                Model.GetXCellSizes(), Model.GetYCellSizes(), Model.GetZCellSizes(),
                SourceObserve, Model.GetConductivities(),
                Model.GetBackgroundConductivities(), Model.GetBackgroundThicknesses(),
                true);
            //write an empty source file for the second source polarization
            WriteSourceFile((EdipDirName / sourcebfilename).string(), SourceXIndex,
                SourceYIndex, Model.GetMeasPosZ(), Zeros, Zeros, Model.GetZCoordinates(),
                Model.GetZCellSizes(), nmodx, nmody);
          }
        std::vector<std::complex<double> > Ux1_el, Ux2_el, Uy1_el, Uy2_el, Uz1_el, Uz2_el;
        //calculate the first polarization and read the adjoint fields
        CalcU(EdipName.string(), XPolMoments1, YPolMoments1, Ux1_el, Uy1_el, Uz1_el,
            SourceXIndex, SourceYIndex, Model.GetMeasPosZ(), Model.GetZCoordinates(),
            Model.GetZCellSizes(), nmodx, nmody, nmodz);
        //calculate the second polarization
        CalcU(EdipName.string(), XPolMoments2, YPolMoments2, Ux2_el, Uy2_el, Uz2_el,
            SourceXIndex, SourceYIndex, Model.GetMeasPosZ(), Model.GetZCoordinates(),
            Model.GetZCellSizes(), nmodx, nmody, nmodz);

        //now we calculate the response to magnetic dipole sources
        const std::complex<double> omega_mu =
            -1.0
                / (std::complex<double>(0.0, jif3D::mag_mu) * 2.0
                    * boost::math::constants::pi<double>()
                    * Model.GetFrequencies()[freqindex]);

        fs::path MdipName = TempDir / MakeUniqueName(X3DModel::MDIP, freqindex);
        fs::path MdipDirName = MdipName.string() + dirext;
        //write the files for the magnetic dipole calculation
#pragma omp critical(gradient_writemodel_mdip)
          {
            MakeRunFile(MdipName.string(), MdipDirName.string());
            WriteProjectFile(MdipDirName.string(), CurrFreq, X3DModel::MDIP,
                resultfilename, modelfilename);
            Write3DModelForX3D((MdipDirName / modelfilename).string(),
                Model.GetXCellSizes(), Model.GetYCellSizes(), Model.GetZCellSizes(),
                SourceObserve, Model.GetConductivities(),
                Model.GetBackgroundConductivities(), Model.GetBackgroundThicknesses(),
                true);
            //write an empty source file for the second source polarization
            WriteSourceFile((MdipDirName / sourcebfilename).string(), SourceXIndex,
                SourceYIndex, Model.GetMeasPosZ(), Zeros, Zeros, Model.GetZCoordinates(),
                Model.GetZCellSizes(), nmodx, nmody);
          }
        //make the sources for the magnetic dipoles
        for (size_t j = 0; j < nmeas; ++j)
          {

            size_t offset = freq_start_index + j * 8;
            std::complex<double> Zxx(RawImpedance(offset), RawImpedance(offset + 1)), Zxy(
                RawImpedance(offset + 2), RawImpedance(offset + 3)), Zyx(
                RawImpedance(offset + 4), RawImpedance(offset + 5)), Zyy(
                RawImpedance(offset + 5), RawImpedance(offset + 6));
          }

        std::vector<std::complex<double> > Ux1_mag, Ux2_mag, Uy1_mag, Uy2_mag, Uz1_mag,
            Uz2_mag;
        //calculate the first polarization and read the adjoint fields
        CalcU(MdipName.string(), XPolMoments1, YPolMoments1, Ux1_mag, Uy1_mag, Uz1_mag,
            SourceXIndex, SourceYIndex, Model.GetMeasPosZ(), Model.GetZCoordinates(),
            Model.GetZCellSizes(), nmodx, nmody, nmodz);
        //calculate the second polarization and read the adjoint fields
        CalcU(MdipName.string(), XPolMoments2, YPolMoments2, Ux2_mag, Uy2_mag, Uz2_mag,
            SourceXIndex, SourceYIndex, Model.GetMeasPosZ(), Model.GetZCoordinates(),
            Model.GetZCellSizes(), nmodx, nmody, nmodz);

        const double cell_sizex = Model.GetXCellSizes()[0];
        const double cell_sizey = Model.GetYCellSizes()[0];
        //now we can calculate the gradient for each model cell
        double Volume, gradinc;

        for (size_t j = 0; j < nmod; ++j)
          {
            Volume = cell_sizex * cell_sizey * Model.GetZCellSizes()[j % nmodz];
            //this is an implementation of eq. 14 in Avdeev and Avdeeva
            //we make the update of the gradient atomic, to avoid
            //race conditions
            gradinc = std::real(
                (Ux1_el[j] + Ux1_mag[j]) * Ex1_all[j]
                    + (Uy1_el[j] + Uy1_mag[j]) * Ey1_all[j]
                    + (Uz1_el[j] + Uz1_mag[j]) * Ez1_all[j]
                    + (Ux2_el[j] + Ux2_mag[j]) * Ex2_all[j]
                    + (Uy2_el[j] + Uy2_mag[j]) * Ey2_all[j]
                    + (Uz2_el[j] + Uz2_mag[j]) * Ez2_all[j]) * Volume;

            Gradient(j) += gradinc;
          }
        return Gradient;
      }

    std::vector<double> AdaptDist(const std::vector<double> &C,
        const jif3D::rvec &RawImpedance, const jif3D::rvec &Misfit)
      {
        std::vector<double> result(C.size(), 0.0);
        const size_t nstat = C.size() / 4;
        const size_t nfreq = RawImpedance.size() / (nstat * 8);
        for (size_t i = 0; i < nfreq; ++i)
          {
            for (size_t j = 0; j < nstat; ++j)
              {
                const size_t offset = (i * nstat + j) * 8;
                result[j * 4] += Misfit(offset) * RawImpedance(offset)
                    + Misfit(offset + 1) * RawImpedance(offset + 1)
                    + Misfit(offset + 4) * RawImpedance(offset + 2)
                    + Misfit(offset + 5) * RawImpedance(offset + 3);
                result[j * 4 + 1] += Misfit(offset) * RawImpedance(offset + 4)
                    + Misfit(offset + 1) * RawImpedance(offset + 5)
                    + Misfit(offset + 4) * RawImpedance(offset + 6)
                    + Misfit(offset + 5) * RawImpedance(offset + 7);
                result[j * 4 + 2] += Misfit(offset + 2) * RawImpedance(offset)
                    + Misfit(offset + 3) * RawImpedance(offset + 1)
                    + Misfit(offset + 6) * RawImpedance(offset + 2)
                    + Misfit(offset + 7) * RawImpedance(offset + 3);
                result[j * 4 + 3] += Misfit(offset + 2) * RawImpedance(offset + 4)
                    + Misfit(offset + 3) * RawImpedance(offset + 5)
                    + Misfit(offset + 6) * RawImpedance(offset + 6)
                    + Misfit(offset + 7) * RawImpedance(offset + 7);
              }
          }
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
        //the number of observations in the model file, one for each cell in the layer
        const size_t nobs = nmodx * nmody;
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
                //calculate the gradient for each frequency
                rvec tmp = LQDerivativeFreq(Model, Misfit, C, i);
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
            std::vector<double> CGrad = AdaptDist(C, RawImpedance, Misfit);
            std::copy(CGrad.begin(), CGrad.end(), Gradient.begin() + nmod);
          }
        //if we had some exception inside the openmp region, we throw
        // a generic error message
        if (FatalError)
          throw jif3D::FatalException("Problem in MT gradient calculation.");

        return 2.0 * Gradient;
      }
  }
