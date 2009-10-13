//============================================================================
// Name        : X3DMTCalculator.cpp
// Author      : Jul 7, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================

#include <cassert>
#include <complex>
#include <fstream>
#include <map>
#include <omp.h>
#include <boost/assign/list_of.hpp> // for 'map_list_of()'
#include <boost/multi_array.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/algorithm/string.hpp>
#include "../Global/FatalException.h"
#include "../Global/convert.h"
#include "X3DMTCalculator.h"
#include "ReadWriteX3D.h"
#include "MTEquations.h"

namespace fs = boost::filesystem;

namespace jiba
  {
    //define some names that are alwats the same
    //either because x3d uses this convention
    //or because we use them in their own directory
    //and want to keep them simple to make sure x3d can handle them
    const std::string modelfilename("x3d.model");
    const std::string resultfilename("x3d.result");
    const std::string emaname = resultfilename + "0.ema";
    const std::string sourcefilename = modelfilename + "0a.source";
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

    void X3DMTCalculator::CleanUp()
      {
        std::string NameRoot(ObjectID());
        fs::directory_iterator end_itr; // default construction yields past-the-end
        for (fs::directory_iterator itr(fs::current_path()); itr != end_itr; ++itr)
          {
            if (boost::algorithm::starts_with(itr->leaf(), NameRoot))
              {
                fs::remove_all(itr->leaf());
              }
          }
      }

    X3DMTCalculator::X3DMTCalculator()
      {

      }

    X3DMTCalculator::~X3DMTCalculator()
      {
        CleanUp();
      }

    //check that the .hnk file for x3d are in a certain directory
    bool CheckHNK(const boost::filesystem::path &TargetDir)
      {
        return fs::exists(TargetDir / "ndec15.hnk") && fs::exists(TargetDir
            / "ndec20.hnk") && fs::exists(TargetDir / "ndec30.hnk")
            && fs::exists(TargetDir / "ndec40.hnk");
      }

    //copy the .hnk files for x3d from SourceDir to TargetDir
    void CopyHNK(const boost::filesystem::path &SourceDir,
        const boost::filesystem::path &TargetDir)
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

    std::string X3DMTCalculator::ObjectID()
      {
        return std::string("p" + jiba::stringify(getpid()) + jiba::stringify(
            this));
      }

    std::string X3DMTCalculator::MakeUniqueName(X3DModel::ProblemType Type,
        const size_t FreqIndex)
      {
        //we assemble the name from the id of the process
        //and the address of the current object
        std::string result(ObjectID());
        //the type of calculation
        result += Extension.find(Type)->second;
        //and the frequency index
        result += jiba::stringify(FreqIndex);
        return result;
      }
    //create a script that changes to the correct directory
    //and executes x3d in that directory
    void MakeRunFile(const std::string &NameRoot)
      {
        std::string DirName = NameRoot + dirext;
        std::string RunFileName = NameRoot + runext;
        fs::create_directory(DirName);
        std::ofstream runfile;
        runfile.open(RunFileName.c_str());
        runfile << "#!/bin/bash\n";
        runfile << "cd " << DirName << "\n";
        runfile << "x3d > /dev/null\n";
        runfile << "cd ..\n";
        runfile.close();
        //it is important to include the std:: namespace specification
        //for the system call, otherwise the GNU compiler picks up
        //a version from the c library that gives trouble in threaded environments
        if (std::system((std::string("chmod u+x ./") + RunFileName).c_str()))
          throw FatalException("Cannot make script executable !");
        CopyHNK(fs::current_path(), DirName);
      }

    //execute the script that runs x3d
    void RunX3D(const std::string &NameRoot)
      {
        const std::string runname = "./" + NameRoot + runext;
        if (std::system(runname.c_str()))
          throw FatalException("Cannot execute run script: " + runname);
      }
    //remove the script file and the directory that contains the results
    //produced by x3d
    void CleanFiles(const std::string &NameRoot)
      {
#pragma omp critical
          {
            //fs::remove_all(NameRoot + dirext);
            //fs::remove_all(NameRoot + runext);
          }
      }

    rvec X3DMTCalculator::Calculate(const X3DModel &Model)
      {
        //we define nfreq as int to make the compiler happy in the openmp loop
        const int nfreq = Model.GetFrequencies().size();
        const size_t nmeas = Model.GetMeasPosX().size();
        const size_t nmodx = Model.GetXCoordinates().size();

        jiba::rvec result(nmeas * nfreq * 8);
        result.clear();
        omp_lock_t lck;
        omp_init_lock(&lck);
        //we parallelize by frequency, this is relatively simple
        //but we can use up to 20 processors for typical MT problems
        // as we do not have the source for x3d, this is our only possibility anyway
#pragma omp parallel for default(shared)
        for (int i = 0; i < nfreq; ++i)
          {
            std::string RootName = MakeUniqueName(X3DModel::MT, i);
            std::string DirName = RootName + dirext;
            std::vector<double> CurrFreq(1, Model.GetFrequencies()[i]);
            //writing out files causes problems in parallel
            // so we make sure it is done one at a time
#pragma omp critical
              {
                MakeRunFile(RootName);
                WriteProjectFile(DirName, CurrFreq, X3DModel::MT,
                    resultfilename, modelfilename);
                Write3DModelForX3D(DirName + "/" + modelfilename,
                    Model.GetXCellSizes(), Model.GetYCellSizes(),
                    Model.GetZCellSizes(), Model.GetConductivities(),
                    Model.GetBackgroundConductivities(),
                    Model.GetBackgroundThicknesses());
              }
            //run x3d in parallel
            RunX3D(RootName);
            std::vector<std::complex<double> > Ex1, Ex2, Ey1, Ey2, Hx1, Hx2,
                Hy1, Hy2;
            std::complex<double> Zxx, Zxy, Zyx, Zyy;
            //read in the electric and magnetic field at the observe sites
            ReadEMO(DirName + "/" + emoAname, Ex1, Ey1, Hx1, Hy1);
            ReadEMO(DirName + "/" + emoBname, Ex2, Ey2, Hx2, Hy2);
            const size_t freq_index = nmeas * i * 8;
            //calculate impedances from the field spectra for all measurement sites
            for (size_t j = 0; j < nmeas; ++j)
              {
                //find out where our site is located in the model
                boost::array<ThreeDModelBase::t3DModelData::index, 3>
                    StationIndex = Model.FindAssociatedIndices(
                        Model.GetMeasPosX()[j], Model.GetMeasPosY()[j],
                        Model.GetMeasPosZ()[j]);
                //at the moment we ignore the depth/elevation of the site
                const size_t offset = StationIndex[1] * nmodx + StationIndex[0];
                const size_t meas_index = freq_index + j * 8;
                FieldsToImpedance(Ex1[offset], Ex2[offset], Ey1[offset],
                    Ey2[offset], Hx1[offset], Hx2[offset], Hy1[offset],
                    Hy2[offset], Zxx, Zxy, Zyx, Zyy);
                //in order to guarantee a coherent state of the array
                //we lock the access before writing
                //update the values and then unlock
                omp_set_lock(&lck);
                result(meas_index) = Zxx.real();
                result(meas_index + 1) = Zxx.imag();
                result(meas_index + 2) = Zxy.real();
                result(meas_index + 3) = Zxy.imag();
                result(meas_index + 4) = Zyx.real();
                result(meas_index + 5) = Zyx.imag();
                result(meas_index + 6) = Zyy.real();
                result(meas_index + 7) = Zyy.imag();
                omp_unset_lock(&lck);
              }
            //finished with one frequency
          }
        omp_destroy_lock(&lck);
        return result;
      }

    cmat CalcATimesH(const cmat &A, const cmat &H)
      {
        cmat result(2, 2);
        const std::complex<double> magdet = 1. / (H(0, 0) * H(1, 1) - H(0, 1)
            * H(1, 0));
        result(0, 0) = magdet * (conj(A(0, 0)) * H(1, 1) - conj(A(0, 1)) * H(0,
            1));
        result(0, 1) = magdet * (-conj(A(0, 0)) * H(1, 0) + conj(A(0, 1)) * H(
            0, 0));
        result(1, 0) = magdet * (conj(A(1, 0)) * H(1, 1) - conj(A(1, 1)) * H(0,
            1));
        result(1, 1) = magdet * (-conj(A(1, 0)) * H(1, 0) + conj(A(1, 1)) * H(
            0, 0));
        return result;
      }

    cmat MisfitToA(const rvec &Misfit, const size_t startindex)
      {
        assert(startindex <= Misfit.size()-8);
        cmat result(2, 2);
        result(0, 0) = std::complex<double>(Misfit(startindex), Misfit(
            startindex + 1));
        result(0, 1) = std::complex<double>(Misfit(startindex + 2), Misfit(
            startindex + 3));
        result(1, 0) = std::complex<double>(Misfit(startindex + 4), Misfit(
            startindex + 5));
        result(1, 1) = std::complex<double>(Misfit(startindex + 6), Misfit(
            startindex + 7));
        return result;
      }

    cmat MakeH(const std::complex<double> &Hx1,
        const std::complex<double> &Hx2, const std::complex<double> &Hy1,
        const std::complex<double> &Hy2)
      {
        cmat result(2, 2);
        result(0, 0) = Hx1;
        result(0, 1) = Hx2;
        result(1, 0) = Hy1;
        result(1, 1) = Hy2;
        return result;
      }

    void CalcHext(const std::complex<double> &omega_mu,
        std::complex<double> &Xp1, std::complex<double> &Xp2, std::complex<
            double> &Yp1, std::complex<double> &Yp2,
        const std::complex<double> &Zxx, const std::complex<double> &Zxy,
        const std::complex<double> &Zyx, const std::complex<double> &Zyy)
      {
        std::complex<double> temp1 = conj(Zxx) * Xp1 + conj(Zyx) * Yp1;
        std::complex<double> temp2 = conj(Zxy) * Xp1 + conj(Zyy) * Yp1;
        Xp1 = temp1 * omega_mu;
        Yp1 = temp2 * omega_mu;
        temp1 = conj(Zxx) * Xp2 + conj(Zyx) * Yp2;
        temp2 = conj(Zxy) * Xp2 + conj(Zyy) * Yp2;
        Xp2 = temp1 * omega_mu;
        Yp2 = temp2 * omega_mu;
      }

    void CalcU(const std::string &RootName, const boost::multi_array<
        std::complex<double>, 2> &XPolMoments, const boost::multi_array<
        std::complex<double>, 2> &YPolMoments,
        std::vector<std::complex<double> > &Ux, std::vector<
            std::complex<double> > &Uy, std::vector<std::complex<double> > &Uz,
        const size_t ncellsx, const size_t ncellsy, const size_t ncellsz)
      {
        std::string DirName = RootName + dirext + "/";
#pragma omp critical
          {
            WriteSourceFile(DirName + sourcefilename, 0.0, XPolMoments,
                YPolMoments);
          }
        RunX3D(RootName);
        ReadEMA(DirName + emaname, Ux, Uy, Uz, ncellsx, ncellsy, ncellsz);
        //boost::filesystem::remove_all(emaname);
      }

    rvec X3DMTCalculator::LQDerivative(const X3DModel &Model,
        const rvec &Misfit)
      {
        //we define nfreq as int to make the compiler happy in the openmp loop
        const int nfreq = Model.GetFrequencies().size();
        const size_t ncellsx = Model.GetConductivities().shape()[0];
        const size_t ncellsy = Model.GetConductivities().shape()[1];
        const size_t ncellsz = Model.GetConductivities().shape()[2];
        const size_t nobs = ncellsx * ncellsy;
        const size_t nmeas = Model.GetMeasPosX().size();
        const size_t ndata = nmeas * nfreq * 8;
        const size_t nmod = ncellsx * ncellsy * ncellsz;
        assert(Misfit.size() == ndata);
        jiba::rvec Gradient(nmod);
        Gradient.clear();
        //we parallelize the gradient calculation by frequency
        //see also the comments for the forward calculation
#pragma omp parallel for default(shared)
        for (int i = 0; i < nfreq; ++i)
          {
            std::string ForwardName = MakeUniqueName(X3DModel::MT, i);
            std::string ForwardDirName = ForwardName + dirext + "/";

            //read the fields from the forward calculation
            std::vector<std::complex<double> > Ex1_obs, Ex2_obs, Ey1_obs,
                Ey2_obs, Hx1_obs, Hx2_obs, Hy1_obs, Hy2_obs;
            ReadEMO(ForwardDirName + emoAname, Ex1_obs, Ey1_obs, Hx1_obs,
                Hy1_obs);
            ReadEMO(ForwardDirName + emoBname, Ex2_obs, Ey2_obs, Hx2_obs,
                Hy2_obs);
            assert(Ex1_obs.size()==nobs);
            assert(Ex2_obs.size()==nobs);
            std::vector<std::complex<double> > Ex1_all, Ex2_all, Ey1_all,
                Ey2_all, Ez1_all, Ez2_all;
            //for the gradient calculation we also need the electric fields
            //at all cells in the model for the two source polarizations of
            //the forward calculations
            ReadEMA(ForwardDirName + emaAname, Ex1_all, Ey1_all, Ez1_all,
                ncellsx, ncellsy, ncellsz);
            ReadEMA(ForwardDirName + emaBname, Ex2_all, Ey2_all, Ez2_all,
                ncellsx, ncellsy, ncellsz);

            //create variables for the adjoint field calculation
            const size_t freq_index = nmeas * i * 8;
            boost::multi_array<std::complex<double>, 2> XPolMoments1(
                boost::extents[ncellsx][ncellsy]), YPolMoments1(
                boost::extents[ncellsx][ncellsy]);
            boost::multi_array<std::complex<double>, 2> XPolMoments2(
                boost::extents[ncellsx][ncellsy]), YPolMoments2(
                boost::extents[ncellsx][ncellsy]);
            boost::multi_array<std::complex<double>, 2> Zeros(
                boost::extents[ncellsx][ncellsy]);
            //we only calculate sources for the observe sites, so
            //we make sure everything else is zero
            std::fill(Zeros.origin(), Zeros.origin() + Zeros.num_elements(),
                0.0);
            std::fill(XPolMoments1.origin(), XPolMoments1.origin() + nobs, 0.0);
            std::fill(XPolMoments2.origin(), XPolMoments2.origin() + nobs, 0.0);
            std::fill(YPolMoments1.origin(), YPolMoments1.origin() + nobs, 0.0);
            std::fill(YPolMoments2.origin(), YPolMoments2.origin() + nobs, 0.0);
            //make the sources for the electric dipoles
            for (size_t j = 0; j < nmeas; ++j)
              {
                boost::array<ThreeDModelBase::t3DModelData::index, 3>
                    StationIndex = Model.FindAssociatedIndices(
                        Model.GetMeasPosX()[j], Model.GetMeasPosY()[j],
                        Model.GetMeasPosZ()[j]);
                //at the moment we ignore the depth/elevation of the site
                const size_t offset = StationIndex[1] * ncellsx
                    + StationIndex[0];

                const size_t siteindex = freq_index + j * 8;
                //this is an implementation of eq. 12 in Avdeev and Avdeeva
                //we do not have any beta, as this is part of the misfit
                cmat j_ext = CalcATimesH(MisfitToA(Misfit, siteindex), MakeH(
                    Hx1_obs[offset], Hx2_obs[offset], Hy1_obs[offset],
                    Hy2_obs[offset]));
                XPolMoments1.data()[offset] = conj(j_ext(0, 0));
                YPolMoments1.data()[offset] = conj(j_ext(1, 0));
                XPolMoments2.data()[offset] = conj(j_ext(0, 1));
                YPolMoments2.data()[offset] = conj(j_ext(1, 1));
              }
            //we only want to calculate for one frequency
            //so our vector has just 1 element
            std::vector<double> CurrFreq(1, 0.0);
            CurrFreq[0] = Model.GetFrequencies()[i];
            std::string EdipName = MakeUniqueName(X3DModel::EDIP, i);
            std::string EdipDirName = EdipName + dirext + "/";
            //again we have to write out some file for the electric
            //dipole calculation with x3d, this shouldn't be done
            //in parallel
#pragma omp critical
              {
                MakeRunFile(EdipName);
                WriteProjectFile(EdipDirName, CurrFreq, X3DModel::EDIP,
                    resultfilename, modelfilename);
                Write3DModelForX3D(EdipDirName + modelfilename,
                    Model.GetXCellSizes(), Model.GetYCellSizes(),
                    Model.GetZCellSizes(), Model.GetConductivities(),
                    Model.GetBackgroundConductivities(),
                    Model.GetBackgroundThicknesses());
                //write an empty source file for the second source polarization
                WriteSourceFile(EdipDirName + modelfilename + "0b.source", 0.0,
                    Zeros, Zeros);
              }
            std::vector<std::complex<double> > Ux1_el, Ux2_el, Uy1_el, Uy2_el,
                Uz1_el, Uz2_el;
            //calculate the first polarization and read the adjoint fields
            CalcU(EdipName, XPolMoments1, YPolMoments1, Ux1_el, Uy1_el, Uz1_el,
                ncellsx, ncellsy, ncellsz);
            //calculate the second polarization
            CalcU(EdipName, XPolMoments2, YPolMoments2, Ux2_el, Uy2_el, Uz2_el,
                ncellsx, ncellsy, ncellsz);


            //now we calculate the response to magnetic dipole sources
            const std::complex<double> omega_mu = -1.0 / (std::complex<double>(
                0.0, jiba::mag_mu) * 2.0 * M_PI * Model.GetFrequencies()[i]);

            std::string MdipName = MakeUniqueName(X3DModel::MDIP, i);
            std::string MdipDirName = MdipName + dirext + "/";
            //write the files for the magnetic dipole calculation
#pragma omp critical
              {
                MakeRunFile(MdipName);
                WriteProjectFile(MdipDirName, CurrFreq, X3DModel::MDIP,
                    resultfilename, modelfilename);
                Write3DModelForX3D(MdipDirName + modelfilename,
                    Model.GetXCellSizes(), Model.GetYCellSizes(),
                    Model.GetZCellSizes(), Model.GetConductivities(),
                    Model.GetBackgroundConductivities(),
                    Model.GetBackgroundThicknesses());
              }
            for (size_t j = 0; j < nobs; ++j)
              {
                cmat Z(2, 2);
                FieldsToImpedance(Ex1_obs[j], Ex2_obs[j], Ey1_obs[j],
                    Ey2_obs[j], Hx1_obs[j], Hx2_obs[j], Hy1_obs[j], Hy2_obs[j],
                    Z(0, 0), Z(0, 1), Z(1, 0), Z(1, 1));
                CalcHext(omega_mu, XPolMoments1.data()[j],
                    XPolMoments2.data()[j], YPolMoments1.data()[j],
                    YPolMoments2.data()[j], Z(0, 0), Z(0, 1), Z(1, 0), Z(1, 1));
              }

            std::vector<std::complex<double> > Ux1_mag, Ux2_mag, Uy1_mag,
                Uy2_mag, Uz1_mag, Uz2_mag;
            //calculate the first polarization and read the adjoint fields
            CalcU(MdipName, XPolMoments1, YPolMoments1, Ux1_mag, Uy1_mag,
                Uz1_mag, ncellsx, ncellsy, ncellsz);
            //calculate the second polarization and read the adjoint fields
            CalcU(MdipName, XPolMoments2, YPolMoments2, Ux2_mag, Uy2_mag,
                Uz2_mag, ncellsx, ncellsy, ncellsz);

            const double cell_sizex = Model.GetXCellSizes()[0];
            const double cell_sizey = Model.GetYCellSizes()[0];
            //now we can calculate the gradient for each model cell
            for (size_t j = 0; j < nmod; ++j)
              {
                const double Volume = cell_sizex * cell_sizey
                    * Model.GetZCellSizes()[j % ncellsz];
                //this is an implementation of eq. 14 in Avdeev and Avdeeva
                //we make the update of the gradient atomic, to avoid
                //race conditions
#pragma omp atomic
                Gradient(j) += std::real((Ux1_el[j] + Ux1_mag[j]) * Ex1_all[j]
                    + (Uy1_el[j] + Uy1_mag[j]) * Ey1_all[j] + (Uz1_el[j]
                    + Uz1_mag[j]) * Ez1_all[j] + (Ux2_el[j] + Ux2_mag[j])
                    * Ex2_all[j] + (Uy2_el[j] + Uy2_mag[j]) * Ey2_all[j]
                    + (Uz2_el[j] + Uz2_mag[j]) * Ez2_all[j]) * Volume;
              }
            //clean the files
            CleanFiles(ForwardName);
            CleanFiles(EdipName);
            CleanFiles(MdipName);
            //finished with one frequency
          }
        //plot the gradient for each model cell in a vtk file
        //this is purely for debugging purposes and might be removed in the future
        X3DModel GradMod(Model);
        std::copy(Gradient.begin(), Gradient.end(),
            GradMod.SetConductivities().origin());
        GradMod.WriteVTK("grad.vtk");
        return Gradient;
      }
  }
