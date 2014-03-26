/*
 * X3DFreqFunctions.cpp
 *
 *  Created on: Mar 17, 2014
 *      Author: mmoorkamp
 */
#include "X3DFreqFunctions.h"
#include "../Global/FatalException.h"
#include "../ModelBase/CellBoundaries.h"

#include "ReadWriteX3D.h"
#include "MTUtils.h"
#include "MTEquations.h"
#include <iostream>
namespace fs = boost::filesystem;

using namespace jif3D;
ForwardResult CalculateFrequency(const ForwardInfo &Info)
  {
    const size_t nmeas = Info.Model.GetMeasPosX().size();
    const size_t nmodx = Info.Model.GetXCoordinates().size();
    const size_t nmody = Info.Model.GetYCoordinates().size();
    ForwardResult result;
    result.DistImpedance.resize(nmeas * 8);
    result.RawImpedance.resize(nmeas * 8);

    fs::path TempDir(Info.TempDirName);
    fs::path RootName = TempDir
        / MakeUniqueName(Info.NameRoot, X3DModel::MT, Info.freqindex);
    fs::path DirName = RootName.string() + dirext;
    std::vector<double> CurrFreq(1, Info.Model.GetFrequencies()[Info.freqindex]);
    std::vector<double> ShiftDepth;
    std::vector<size_t> MeasDepthIndices;
    //construct a vector of indices of unique station depths
    size_t nlevels = ConstructDepthIndices(MeasDepthIndices, ShiftDepth, Info.Model);
    //writing out files causes problems in parallel
    // so we make sure it is done one at a time
#pragma omp critical(forward_write_files)
      {
        MakeRunFile(RootName.string(), DirName.string(), Info.X3DName);
        WriteProjectFile(DirName.string(), CurrFreq, X3DModel::MT, resultfilename,
            modelfilename);
        Write3DModelForX3D((DirName / modelfilename).string(), Info.Model.GetXCellSizes(),
            Info.Model.GetYCellSizes(), Info.Model.GetZCellSizes(), ShiftDepth,
            Info.Model.GetConductivities(), Info.Model.GetBackgroundConductivities(),
            Info.Model.GetBackgroundThicknesses());
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
            Info.Model.FindAssociatedIndices(Info.Model.GetMeasPosX()[j],
                Info.Model.GetMeasPosY()[j], Info.Model.GetMeasPosZ()[j]);
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
        FieldsToImpedance(Ex1[offset], Ex2[offset], Ey1[offset], Ey2[offset], Hx1[offset],
            Hx2[offset], Hy1[offset], Hy2[offset], Zxx, Zxy, Zyx, Zyy);
        //result is a local array for this frequency
        //so we can directly use it even in a threaded environment
        const size_t meas_index = j * 8;
        const size_t site_index = j * 4;

        result.DistImpedance(meas_index) = Info.C[site_index] * Zxx.real()
            + Info.C[site_index + 1] * Zyx.real();
        result.DistImpedance(meas_index + 1) = Info.C[site_index] * Zxx.imag()
            + Info.C[site_index + 1] * Zyx.imag();
        result.DistImpedance(meas_index + 2) = Info.C[site_index] * Zxy.real()
            + Info.C[site_index + 1] * Zyy.real();
        result.DistImpedance(meas_index + 3) = Info.C[site_index] * Zxy.imag()
            + Info.C[site_index + 1] * Zyy.imag();
        result.DistImpedance(meas_index + 4) = Info.C[site_index + 3] * Zyx.real()
            + Info.C[site_index + 2] * Zxx.real();
        result.DistImpedance(meas_index + 5) = Info.C[site_index + 3] * Zyx.imag()
            + Info.C[site_index + 2] * Zxx.imag();
        result.DistImpedance(meas_index + 6) = Info.C[site_index + 3] * Zyy.real()
            + Info.C[site_index + 2] * Zxy.real();
        result.DistImpedance(meas_index + 7) = Info.C[site_index + 3] * Zyy.imag()
            + Info.C[site_index + 2] * Zxy.real();
        result.RawImpedance(meas_index) = Zxx.real();
        result.RawImpedance(meas_index + 1) = Zxx.imag();
        result.RawImpedance(meas_index + 2) = Zxy.real();
        result.RawImpedance(meas_index + 3) = Zxy.imag();
        result.RawImpedance(meas_index + 4) = Zyx.real();
        result.RawImpedance(meas_index + 5) = Zyx.imag();
        result.RawImpedance(meas_index + 6) = Zyy.real();
        result.RawImpedance(meas_index + 7) = Zyy.imag();

      }

    return result;
  }

jif3D::rvec LQDerivativeFreq(const ForwardInfo &Info, const jif3D::rvec &Misfit,
    const jif3D::rvec &RawImpedance)
  {
    //a few commonly used quantities for shorter notation
    const size_t nmodx = Info.Model.GetConductivities().shape()[0];
    const size_t nmody = Info.Model.GetConductivities().shape()[1];
    const size_t nmodz = Info.Model.GetConductivities().shape()[2];
    //the number of observations in the model file, one for each cell in the layer
    const size_t nobs = nmodx * nmody;
    //the number of measurement sites
    const size_t nmeas = Info.Model.GetMeasPosX().size();
    const size_t nmod = nmodx * nmody * nmodz;
    std::vector<double> ShiftDepth;
    //for the controlled source calculations we do not actually
    //need any observe layers as we are only interested in the
    //anomalous fields
    std::vector<double> SourceObserve(1, 0.0);
    std::vector<size_t> MeasDepthIndices;
    size_t nlevels = ConstructDepthIndices(MeasDepthIndices, ShiftDepth, Info.Model);
    jif3D::rvec Gradient(nmod, 0.0);

    fs::path TempDir(Info.TempDirName);
    fs::path ForwardDirName = TempDir
        / (MakeUniqueName(Info.NameRoot, X3DModel::MT, Info.freqindex) + dirext);
    if (!fs::is_directory(ForwardDirName))
      throw FatalException(
          "In X3D gradient calculation, directory does not exist: "
              + ForwardDirName.string());

    //read the fields from the forward calculation
    std::vector<std::complex<double> > Ex1_obs, Ex2_obs, Ey1_obs, Ey2_obs, Hx1_obs,
        Hx2_obs, Hy1_obs, Hy2_obs;
#pragma omp critical(gradient_reademo)
      {
        ReadEMO((ForwardDirName / emoAname).string(), Ex1_obs, Ey1_obs, Hx1_obs, Hy1_obs);
        ReadEMO((ForwardDirName / emoBname).string(), Ex2_obs, Ey2_obs, Hx2_obs, Hy2_obs);
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
        ReadEMA((ForwardDirName / emaAname).string(), Ex1_all, Ey1_all, Ez1_all, nmodx,
            nmody, nmodz);
        ReadEMA((ForwardDirName / emaBname).string(), Ex2_all, Ey2_all, Ez2_all, nmodx,
            nmody, nmodz);
      }
    //create variables for the adjoint field calculation
    const size_t freq_start_index = nmeas * Info.freqindex * 8;
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
            Info.Model.FindAssociatedIndices(Info.Model.GetMeasPosX()[j],
                Info.Model.GetMeasPosY()[j], Info.Model.GetMeasPosZ()[j]);

        const size_t offset = (nmodx * nmody) * MeasDepthIndices[j]
            + StationIndex[0] * nmody + StationIndex[1];
        SourceXIndex.at(j) = StationIndex[0];
        SourceYIndex.at(j) = StationIndex[1];

        //this is an implementation of eq. 12 in Avdeev and Avdeeva
        //we do not have any beta, as this is part of the misfit
//            cmat j_ext = CalcATimesH(MisfitToA(Misfit, siteindex),
//                MakeH(Hx1_obs[offset], Hx2_obs[offset], Hy1_obs[offset], Hy2_obs[offset]));
        cmat j_ext = CalcEExt(Misfit, Info.C, j, freq_start_index, Hx1_obs[offset],
            Hx2_obs[offset], Hy1_obs[offset], Hy2_obs[offset]);
        XPolMoments1.at(j) = conj(j_ext(0, 0));
        YPolMoments1.at(j) = conj(j_ext(1, 0));
        XPolMoments2.at(j) = conj(j_ext(0, 1));
        YPolMoments2.at(j) = conj(j_ext(1, 1));
      }
    //we only want to calculate for one frequency
    //so our vector has just 1 element
    std::vector<double> CurrFreq(1, 0.0);
    CurrFreq[0] = Info.Model.GetFrequencies()[Info.freqindex];
    fs::path EdipName = TempDir
        / MakeUniqueName(Info.NameRoot, X3DModel::EDIP, Info.freqindex);
    fs::path EdipDirName = EdipName.string() + dirext;
    //again we have to write out some file for the electric
    //dipole calculation with x3d, this shouldn't be done
    //in parallel
#pragma omp critical(gradient_writemodel_edip)
      {
        MakeRunFile(EdipName.string(), EdipDirName.string(), Info.X3DName);
        WriteProjectFile(EdipDirName.string(), CurrFreq, X3DModel::EDIP, resultfilename,
            modelfilename);
        Write3DModelForX3D((EdipDirName / modelfilename).string(),
            Info.Model.GetXCellSizes(), Info.Model.GetYCellSizes(),
            Info.Model.GetZCellSizes(), SourceObserve, Info.Model.GetConductivities(),
            Info.Model.GetBackgroundConductivities(),
            Info.Model.GetBackgroundThicknesses(), true);
        //write an empty source file for the second source polarization
        WriteSourceFile((EdipDirName / sourcebfilename).string(), SourceXIndex,
            SourceYIndex, Info.Model.GetMeasPosZ(), Zeros, Zeros,
            Info.Model.GetZCoordinates(), Info.Model.GetZCellSizes(), nmodx, nmody);
      }
    std::vector<std::complex<double> > Ux1_el, Ux2_el, Uy1_el, Uy2_el, Uz1_el, Uz2_el;
    //calculate the first polarization and read the adjoint fields
    CalcU(EdipName.string(), XPolMoments1, YPolMoments1, Ux1_el, Uy1_el, Uz1_el,
        SourceXIndex, SourceYIndex, Info.Model.GetMeasPosZ(),
        Info.Model.GetZCoordinates(), Info.Model.GetZCellSizes(), nmodx, nmody, nmodz);
    //calculate the second polarization
    CalcU(EdipName.string(), XPolMoments2, YPolMoments2, Ux2_el, Uy2_el, Uz2_el,
        SourceXIndex, SourceYIndex, Info.Model.GetMeasPosZ(),
        Info.Model.GetZCoordinates(), Info.Model.GetZCellSizes(), nmodx, nmody, nmodz);

    //now we calculate the response to magnetic dipole sources
    fs::path MdipName = TempDir
        / MakeUniqueName(Info.NameRoot, X3DModel::MDIP, Info.freqindex);
    fs::path MdipDirName = MdipName.string() + dirext;
    //write the files for the magnetic dipole calculation
#pragma omp critical(gradient_writemodel_mdip)
      {
        MakeRunFile(MdipName.string(), MdipDirName.string(), Info.X3DName);
        WriteProjectFile(MdipDirName.string(), CurrFreq, X3DModel::MDIP, resultfilename,
            modelfilename);
        Write3DModelForX3D((MdipDirName / modelfilename).string(),
            Info.Model.GetXCellSizes(), Info.Model.GetYCellSizes(),
            Info.Model.GetZCellSizes(), SourceObserve, Info.Model.GetConductivities(),
            Info.Model.GetBackgroundConductivities(),
            Info.Model.GetBackgroundThicknesses(), true);
        //write an empty source file for the second source polarization
        WriteSourceFile((MdipDirName / sourcebfilename).string(), SourceXIndex,
            SourceYIndex, Info.Model.GetMeasPosZ(), Zeros, Zeros,
            Info.Model.GetZCoordinates(), Info.Model.GetZCellSizes(), nmodx, nmody);
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
        SourceXIndex, SourceYIndex, Info.Model.GetMeasPosZ(),
        Info.Model.GetZCoordinates(), Info.Model.GetZCellSizes(), nmodx, nmody, nmodz);
    //calculate the second polarization and read the adjoint fields
    CalcU(MdipName.string(), XPolMoments2, YPolMoments2, Ux2_mag, Uy2_mag, Uz2_mag,
        SourceXIndex, SourceYIndex, Info.Model.GetMeasPosZ(),
        Info.Model.GetZCoordinates(), Info.Model.GetZCellSizes(), nmodx, nmody, nmodz);

    const double cell_sizex = Info.Model.GetXCellSizes()[0];
    const double cell_sizey = Info.Model.GetYCellSizes()[0];
    //now we can calculate the gradient for each model cell
    double Volume, gradinc;

    for (size_t j = 0; j < nmod; ++j)
      {
        Volume = cell_sizex * cell_sizey * Info.Model.GetZCellSizes()[j % nmodz];
        //this is an implementation of eq. 14 in Avdeev and Avdeeva
        //we make the update of the gradient atomic, to avoid
        //race conditions
        gradinc = std::real(
            (Ux1_el[j] + Ux1_mag[j]) * Ex1_all[j] + (Uy1_el[j] + Uy1_mag[j]) * Ey1_all[j]
                + (Uz1_el[j] + Uz1_mag[j]) * Ez1_all[j]
                + (Ux2_el[j] + Ux2_mag[j]) * Ex2_all[j]
                + (Uy2_el[j] + Uy2_mag[j]) * Ey2_all[j]
                + (Uz2_el[j] + Uz2_mag[j]) * Ez2_all[j]) * Volume;

        Gradient(j) += gradinc;
      }
    return Gradient;
  }

#ifdef HAVEHPX

HPX_REGISTER_PLAIN_ACTION(CalculateFrequency_action)
HPX_REGISTER_PLAIN_ACTION( LQDerivativeFreq_action)
#endif
