/*
 * X3DFreqFunctions.cpp
 *
 *  Created on: Mar 17, 2014
 *      Author: mmoorkamp
 */
#include "X3DFreqFunctions.h"
#include "../Global/FatalException.h"
#include "../Global/NumUtil.h"
#include "../ModelBase/CellBoundaries.h"

#include "ReadWriteX3D.h"
#include "MTUtils.h"
#include "MTEquations.h"
#include <iostream>
#include <boost/math/constants/constants.hpp>
namespace fs = boost::filesystem;

using namespace jif3D;
ForwardResult CalculateFrequency(const ForwardInfo &Info, const jif3D::MTData &Data,
    boost::shared_ptr<jif3D::X3DFieldCalculator> Calc)
  {
    //  const size_t nmeas = Info.Model.GetMeasPosX().size();
    const size_t nfreq = Data.GetFrequencies().size();
    const size_t nstats = Data.GetExIndices().size() / nfreq;
    const size_t nmodx = Info.Model.GetData().shape()[0];
    const size_t nmody = Info.Model.GetData().shape()[1];
    const size_t ind_shift = nstats * Info.freqindex;

    ForwardResult result;
    result.DistImpedance.resize(nstats * 8);
    //we store the impedances at the two potential locations for Titan
    //this is needed when inverting Titan data but redundant for MT
    result.RawImpedance.resize(nstats * 8 * 2);

    std::vector<double> ShiftDepth;
    std::vector<size_t> MeasDepthIndices;
    //construct a vector of indices of unique station depths
    size_t nlevels = ConstructDepthIndices(MeasDepthIndices, ShiftDepth, Info.Model,
        Data.GetMeasPosZ());

    Calc->CalculateFields(Info.Model, Data.GetFrequencies(), Data.GetMeasPosZ(),
        Info.freqindex);
    //for Titan24 and other DAS we might have different impedances at the
    //different electric field measurement positions and we need all four of them
    //for MT ZxxEx = ZxxEy and so on, it makes it easier to consider both cases
    //by keeping both impedances and make redundant calculations,
    std::complex<double> ZxxEx, ZxyEx, ZyxEx, ZyyEx;
    std::complex<double> ZxxEy, ZxyEy, ZyxEy, ZyyEy;

    //calculate impedances from the field spectra for all stations
    for (size_t j = 0; j < nstats; ++j)
      {
        boost::array<ThreeDModelBase::t3DModelData::index, 3> StationExIndex =
            Info.Model.FindAssociatedIndices(
                Data.GetMeasPosX()[Data.GetExIndices()[j + ind_shift]],
                Data.GetMeasPosY()[Data.GetExIndices()[j + ind_shift]],
                Data.GetMeasPosZ()[Data.GetExIndices()[j + ind_shift]]);
        boost::array<ThreeDModelBase::t3DModelData::index, 3> StationEyIndex =
            Info.Model.FindAssociatedIndices(
                Data.GetMeasPosX()[Data.GetEyIndices()[j + ind_shift]],
                Data.GetMeasPosY()[Data.GetEyIndices()[j + ind_shift]],
                Data.GetMeasPosZ()[Data.GetEyIndices()[j + ind_shift]]);
        boost::array<ThreeDModelBase::t3DModelData::index, 3> StationHxIndex =
            Info.Model.FindAssociatedIndices(
                Data.GetMeasPosX()[Data.GetHxIndices()[j + ind_shift]],
                Data.GetMeasPosY()[Data.GetHxIndices()[j + ind_shift]],
                Data.GetMeasPosZ()[Data.GetHxIndices()[j + ind_shift]]);
        boost::array<ThreeDModelBase::t3DModelData::index, 3> StationHyIndex =
            Info.Model.FindAssociatedIndices(
                Data.GetMeasPosX()[Data.GetHyIndices()[j + ind_shift]],
                Data.GetMeasPosY()[Data.GetHyIndices()[j + ind_shift]],
                Data.GetMeasPosZ()[Data.GetHyIndices()[j + ind_shift]]);
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
        const size_t offset_Ex = (nmodx * nmody)
            * MeasDepthIndices[Data.GetExIndices()[j + ind_shift]]
            + StationExIndex[0] * nmody + StationExIndex[1];
        const size_t offset_Ey = (nmodx * nmody)
            * MeasDepthIndices[Data.GetEyIndices()[j + ind_shift]]
            + StationEyIndex[0] * nmody + StationEyIndex[1];
        const size_t offset_Hx = (nmodx * nmody)
            * MeasDepthIndices[Data.GetHxIndices()[j + ind_shift]]
            + StationHxIndex[0] * nmody + StationHxIndex[1];
        const size_t offset_Hy = (nmodx * nmody)
            * MeasDepthIndices[Data.GetHyIndices()[j + ind_shift]]
            + StationHyIndex[0] * nmody + StationHyIndex[1];
        //calculate the full impedance at the location of the Ex field measurements
        FieldsToImpedance(Calc->GetEx1()[offset_Ex], Calc->GetEx2()[offset_Ex],
            Calc->GetEy1()[offset_Ex], Calc->GetEy2()[offset_Ex],
            Calc->GetHx1()[offset_Hx], Calc->GetHx2()[offset_Hx],
            Calc->GetHy1()[offset_Hy], Calc->GetHy2()[offset_Hy], ZxxEx, ZxyEx, ZyxEx,
            ZyyEx);
        //calculate the full impedance at the location of the Ey field measurements
        FieldsToImpedance(Calc->GetEx1()[offset_Ey], Calc->GetEx2()[offset_Ey],
            Calc->GetEy1()[offset_Ey], Calc->GetEy2()[offset_Ey],
            Calc->GetHx1()[offset_Hx], Calc->GetHx2()[offset_Hx],
            Calc->GetHy1()[offset_Hy], Calc->GetHy2()[offset_Hy], ZxxEy, ZxyEy, ZyxEy,
            ZyyEy);
        //result is a local array for this frequency
        //so we can directly use it even in a threaded environment
        const size_t meas_index = j * 8;
        const size_t site_index = j * 4;
        const double angle = Data.GetRotAngles().at(j) / 180.0
            * boost::math::constants::pi<double>();
        const double c2 = jif3D::pow2(std::cos(angle));
        const double s2 = jif3D::pow2(std::sin(angle));
        const double sc = std::sin(angle) * cos(angle);
        const double Cxx = Info.C[site_index];
        const double Cxy = Info.C[site_index + 1];
        const double Cyx = Info.C[site_index + 2];
        const double Cyy = Info.C[site_index + 3];
        std::complex<double> Zxx = (c2 * Cxx - sc * Cxy) * ZxxEx
            + (sc * Cxx + c2 * Cxy) * ZyxEx + (s2 * Cyy - sc * Cyx) * ZxxEy
            - (s2 * Cyx + sc * Cyy) * ZyxEy;
        std::complex<double> Zxy = (c2 * Cxx - sc * Cxy) * ZxyEx
            + (sc * Cxx + c2 * Cxy) * ZyyEx + (s2 * Cyy - sc * Cyx) * ZxyEy
            - (s2 * Cyx + sc * Cyy) * ZyyEy;
        std::complex<double> Zyx = (sc * Cxx - s2 * Cxy) * ZxxEx
            + (s2 * Cxx + sc * Cxy) * ZyxEx + (c2 * Cyx - sc * Cyy) * ZxxEy
            + (sc * Cyx + c2 * Cyy) * ZyxEy;
        std::complex<double> Zyy = (sc * Cxx - s2 * Cxy) * ZxyEx
            + (s2 * Cxx + sc * Cxy) * ZyyEx + (c2 * Cyx - sc * Cyy) * ZxyEy
            + (sc * Cyx + c2 * Cyy) * ZyyEy;
        result.DistImpedance[meas_index] = Zxx.real();
        result.DistImpedance[meas_index + 1] = Zxx.imag();
        result.DistImpedance[meas_index + 2] = Zxy.real();
        result.DistImpedance[meas_index + 3] = Zxy.imag();
        result.DistImpedance[meas_index + 4] = Zyx.real();
        result.DistImpedance[meas_index + 5] = Zyx.imag();
        result.DistImpedance[meas_index + 6] = Zyy.real();
        result.DistImpedance[meas_index + 7] = Zyy.imag();

        const size_t rindex = j * 16;

        result.RawImpedance[rindex] = ZxxEx.real();
        result.RawImpedance[rindex + 1] = ZxxEx.imag();
        result.RawImpedance[rindex + 2] = ZxyEx.real();
        result.RawImpedance[rindex + 3] = ZxyEx.imag();
        result.RawImpedance[rindex + 4] = ZyxEx.real();
        result.RawImpedance[rindex + 5] = ZyxEx.imag();
        result.RawImpedance[rindex + 6] = ZyyEx.real();
        result.RawImpedance[rindex + 7] = ZyyEx.imag();
        result.RawImpedance[rindex + 8] = ZxxEy.real();
        result.RawImpedance[rindex + 9] = ZxxEy.imag();
        result.RawImpedance[rindex + 10] = ZxyEy.real();
        result.RawImpedance[rindex + 11] = ZxyEy.imag();
        result.RawImpedance[rindex + 12] = ZyxEy.real();
        result.RawImpedance[rindex + 13] = ZyxEy.imag();
        result.RawImpedance[rindex + 14] = ZyyEy.real();
        result.RawImpedance[rindex + 15] = ZyyEy.imag();
      }

    return result;
  }

GradResult LQDerivativeFreq(const ForwardInfo &Info, const jif3D::MTData &Data,
    const GradInfo &GI, boost::shared_ptr<jif3D::X3DFieldCalculator> Calc)
  {
    //a few commonly used quantities for shorter notation
    const size_t nmodx = Info.Model.GetConductivities().shape()[0];
    const size_t nmody = Info.Model.GetConductivities().shape()[1];
    const size_t nmodz = Info.Model.GetConductivities().shape()[2];
    //the number of observations in the fields files, one for each cell in the layer
    const size_t nobs = nmodx * nmody;
    //the number of measurement sites
    const size_t nmeas = Data.GetMeasPosX().size();
    //the number of transfer functions
    const size_t nstats = Data.GetExIndices().size() / Data.GetFrequencies().size();

    const size_t nmod = nmodx * nmody * nmodz;
    std::vector<double> ShiftDepth;
    //for the controlled source calculations we do not actually
    //need any observe layers as we are only interested in the
    //anomalous fields
    std::vector<double> SourceObserve(1, 0.0);
    std::vector<size_t> MeasDepthIndices;
    size_t nlevels = ConstructDepthIndices(MeasDepthIndices, ShiftDepth, Info.Model,
        Data.GetMeasPosZ());
    jif3D::rvec Gradient(nmod, 0.0);
    jif3D::rvec Misfit(GI.Misfit.size());
    std::copy(GI.Misfit.begin(), GI.Misfit.end(), Misfit.begin());
    fs::path TempDir(Info.TempDirName);
    fs::path ForwardDirName = Calc->GetForwardDirName();

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
    const size_t freq_start_index = nstats * Info.freqindex * 8;
    const size_t ind_shift = nstats * Info.freqindex;

    std::vector<std::complex<double> > EXPolMoments1(2 * nstats, 0.0), EXPolMoments2(
        2 * nstats, 0.0), EYPolMoments1(2 * nstats, 0.0), EYPolMoments2(2 * nstats, 0.0),
        Zeros(2 * nstats, 0.0), HZeros(nstats, 0.0);
    std::vector<std::complex<double> > HXPolMoments1(nstats, 0.0), HXPolMoments2(nstats,
        0.0), HYPolMoments1(nstats, 0.0), HYPolMoments2(nstats, 0.0);
    std::vector<double> XSourceDepth(2 * nstats, 0.0), YSourceDepth(2 * nstats, 0.0),
        HxSourceDepth(nstats), HySourceDepth(nstats, 0.0);
    std::vector<size_t> XSourceXIndex(2 * nstats, 0), XSourceYIndex(2 * nstats, 0),
        YSourceXIndex(2 * nstats, 0), YSourceYIndex(2 * nstats, 0), HxSourceXIndex(nstats,
            0), HxSourceYIndex(nstats, 0), HySourceXIndex(nstats, 0), HySourceYIndex(
            nstats, 0), ZeroIndex(2 * nstats, 0);

    //make the sources for the electric dipoles
    for (size_t j = 0; j < nstats; ++j)
      {
        const int ExInd = Data.GetExIndices()[j + ind_shift];
        const int EyInd = Data.GetEyIndices()[j + ind_shift];
        const int HxInd = Data.GetHxIndices()[j + ind_shift];
        const int HyInd = Data.GetHyIndices()[j + ind_shift];

        boost::array<ThreeDModelBase::t3DModelData::index, 3> StationExIndex =
            Info.Model.FindAssociatedIndices(Data.GetMeasPosX()[ExInd],
                Data.GetMeasPosY()[ExInd], Data.GetMeasPosZ()[ExInd]);

        boost::array<ThreeDModelBase::t3DModelData::index, 3> StationEyIndex =
            Info.Model.FindAssociatedIndices(Data.GetMeasPosX()[EyInd],
                Data.GetMeasPosY()[EyInd], Data.GetMeasPosZ()[EyInd]);

        boost::array<ThreeDModelBase::t3DModelData::index, 3> StationHxIndex =
            Info.Model.FindAssociatedIndices(Data.GetMeasPosX()[HxInd],
                Data.GetMeasPosY()[HxInd], Data.GetMeasPosZ()[HxInd]);

        boost::array<ThreeDModelBase::t3DModelData::index, 3> StationHyIndex =
            Info.Model.FindAssociatedIndices(Data.GetMeasPosX()[HyInd],
                Data.GetMeasPosY()[HyInd], Data.GetMeasPosZ()[HyInd]);
        /*        const size_t offset_Ex = (nmodx * nmody) * MeasDepthIndices[Info.Model.GetExIndices()[j + ind_shift]]
         + StationExIndex[0] * nmody + StationExIndex[1];
         const size_t offset_Ey = (nmodx * nmody) * MeasDepthIndices[Info.Model.GetEyIndices()[j + ind_shift]]
         + StationEyIndex[0] * nmody + StationEyIndex[1]; */
        const size_t offset_Hx = (nmodx * nmody) * MeasDepthIndices[HxInd]
            + StationHxIndex[0] * nmody + StationHxIndex[1];
        const size_t offset_Hy = (nmodx * nmody) * MeasDepthIndices[HyInd]
            + StationHyIndex[0] * nmody + StationHyIndex[1];

        XSourceXIndex.at(2 * j) = StationExIndex[0];
        XSourceYIndex.at(2 * j) = StationExIndex[1];
        XSourceDepth.at(2 * j) = Data.GetMeasPosZ()[ExInd];
        YSourceXIndex.at(2 * j) = XSourceXIndex.at(2 * j);
        YSourceYIndex.at(2 * j) = XSourceYIndex.at(2 * j);
        YSourceDepth.at(2 * j) = XSourceDepth.at(2 * j);

        XSourceXIndex.at(2 * j + 1) = StationEyIndex[0];
        XSourceYIndex.at(2 * j + 1) = StationEyIndex[1];
        XSourceDepth.at(2 * j + 1) = Data.GetMeasPosZ()[EyInd];
        YSourceXIndex.at(2 * j + 1) = XSourceXIndex.at(2 * j + 1);
        YSourceYIndex.at(2 * j + 1) = XSourceYIndex.at(2 * j + 1);
        YSourceDepth.at(2 * j + 1) = XSourceDepth.at(2 * j + 1);

        HxSourceXIndex.at(j) = StationHxIndex[0];
        HxSourceYIndex.at(j) = StationHxIndex[1];

        HySourceXIndex.at(j) = StationHyIndex[0];
        HySourceYIndex.at(j) = StationHyIndex[1];

        HxSourceDepth.at(j) = Data.GetMeasPosZ()[HxInd];
        HySourceDepth.at(j) = Data.GetMeasPosZ()[HyInd];

        //this is an implementation of eq. 12 in Avdeev and Avdeeva
        //we do not have any beta, as this is part of the misfit
        cmat AH = CalcEExt(Misfit, Info.C, j, freq_start_index, Calc->GetHx1()[offset_Hx],
            Calc->GetHx2()[offset_Hx], Calc->GetHy1()[offset_Hy],
            Calc->GetHy2()[offset_Hy]);
        rmat R(2, 2), C1(2, 2), C2(2, 2);
        const double angle = Data.GetRotAngles().at(j) / 180.0 * M_PI;
        R(0, 0) = std::cos(angle);
        R(0, 1) = -1.0 * std::sin(angle);
        R(1, 0) = std::sin(angle);
        R(1, 1) = std::cos(angle);
        C1(0, 0) = Info.C[j * 4];
        C1(0, 1) = 0.0;
        C1(1, 0) = Info.C[j * 4 + 1];
        C1(1, 1) = 0.0;
        C2(0, 1) = Info.C[j * 4 + 2];
        C2(0, 0) = 0.0;
        C2(1, 1) = Info.C[j * 4 + 3];
        C2(1, 0) = 0.0;
        cmat RAH = ublas::prod(ublas::trans(R), AH);
        rmat RC1 = ublas::prod(R, C1);
        rmat RC2 = ublas::prod(R, C2);
        //std::cout << "Matrices: " << std::endl;
        //std::cout << R << " " << C1 << " " << C2 << " " << AH << " " << RAH << " " << RC1 << " " << RC2 << std::endl;
        cmat E1eff = ublas::prod(RC1, RAH); //ublas::prod(ublas::prod(ublas::prod(R , C1), ublas::trans(R)),AH);
        cmat E2eff = ublas::prod(RC2, RAH); //ublas::prod(ublas::prod(ublas::prod(R , C2), ublas::trans(R)),AH)
        EXPolMoments1.at(2 * j) = E1eff(0, 0); //AH(0, 0) * Info.C[j * 4];
        EXPolMoments2.at(2 * j) = E1eff(0, 1); //AH(0, 1) * Info.C[j * 4];
        EYPolMoments1.at(2 * j) = E1eff(1, 0); //AH(0, 0) * Info.C[j * 4 + 1];
        EYPolMoments2.at(2 * j) = E1eff(1, 1); //AH(0, 1) * Info.C[j * 4 + 1];
        EXPolMoments1.at(2 * j + 1) = E2eff(0, 0); //AH(1, 0) * Info.C[j * 4 + 2];
        EXPolMoments2.at(2 * j + 1) = E2eff(0, 1); //AH(1, 1) * Info.C[j * 4 + 2];
        EYPolMoments1.at(2 * j + 1) = E2eff(1, 0); //AH(1, 0) * Info.C[j * 4 + 3];
        EYPolMoments2.at(2 * j + 1) = E2eff(1, 1); //AH(1, 1) * Info.C[j * 4 + 3];
      }
    //we only want to calculate for one frequency
    //so our vector has just 1 element
    std::vector<double> CurrFreq =
      { Data.GetFrequencies()[Info.freqindex] };
    fs::path EdipName = TempDir
        / MakeUniqueName(Info.NameRoot, X3DModel::EDIP, Info.freqindex);
    std::vector<std::string> PolExt =
      { "a", "b" };
    //again we have to write out some file for the electric
    //dipole calculation with x3d, this shouldn't be done
    //in parallel
#pragma omp critical(gradient_writemodel_edip)
      {
        for (auto Ext : PolExt)
          {
            std::string CurrName = EdipName.string() + Ext;
            fs::path CurrDir = EdipName.string() + Ext + dirext;
            MakeRunFile(CurrName, CurrDir.string(), Info.X3DName);
            WriteProjectFile(CurrDir.string(), CurrFreq, X3DModel::EDIP, resultfilename,
                modelfilename, Info.GreenStage1, Info.GreenStage4);
            Write3DModelForX3D((CurrDir / modelfilename).string(),
                Info.Model.GetXCellSizes(), Info.Model.GetYCellSizes(),
                Info.Model.GetZCellSizes(), SourceObserve, Info.Model.GetConductivities(),
                Info.Model.GetBackgroundConductivities(),
                Info.Model.GetBackgroundThicknesses(), true);
            //write an empty source file for the second source polarization
            WriteSourceFile((CurrDir / sourcebfilename).string(), XSourceXIndex,
                XSourceYIndex, XSourceDepth, YSourceXIndex, YSourceYIndex, YSourceDepth,
                YSourceXIndex, YSourceYIndex, YSourceDepth, Zeros, Zeros, Zeros,
                Info.Model.GetZCoordinates(), Info.Model.GetZCellSizes(), nmodx, nmody);
          }
      }
    std::vector<std::complex<double> > Ux1_el, Ux2_el, Uy1_el, Uy2_el, Uz1_el, Uz2_el;
    //calculate the first polarization and read the adjoint fields
#pragma omp task default(shared)
    CalcU(EdipName.string() + PolExt[0], EXPolMoments1, EYPolMoments1, Zeros, Ux1_el,
        Uy1_el, Uz1_el, XSourceXIndex, XSourceYIndex, XSourceDepth, YSourceXIndex,
        YSourceYIndex, YSourceDepth, ZeroIndex, ZeroIndex, YSourceDepth,
        Info.Model.GetZCoordinates(), Info.Model.GetZCellSizes(), nmodx, nmody, nmodz);

    //calculate the second polarization
#pragma omp task default(shared)
    CalcU(EdipName.string() + PolExt[1], EXPolMoments2, EYPolMoments2, Zeros, Ux2_el,
        Uy2_el, Uz2_el, XSourceXIndex, XSourceYIndex, XSourceDepth, YSourceXIndex,
        YSourceYIndex, YSourceDepth, ZeroIndex, ZeroIndex, YSourceDepth,
        Info.Model.GetZCoordinates(), Info.Model.GetZCellSizes(), nmodx, nmody, nmodz);

    //now we calculate the response to magnetic dipole sources
    fs::path MdipName = TempDir
        / MakeUniqueName(Info.NameRoot, X3DModel::MDIP, Info.freqindex);
    //write the files for the magnetic dipole calculation
#pragma omp critical(gradient_writemodel_mdip)
      {
        for (auto Ext : PolExt)
          {
            std::string CurrName = MdipName.string() + Ext;
            fs::path CurrDir = MdipName.string() + Ext + dirext;
            MakeRunFile(CurrName, CurrDir.string(), Info.X3DName);
            WriteProjectFile(CurrDir.string(), CurrFreq, X3DModel::MDIP, resultfilename,
                modelfilename, Info.GreenStage1, Info.GreenStage4);
            Write3DModelForX3D((CurrDir / modelfilename).string(),
                Info.Model.GetXCellSizes(), Info.Model.GetYCellSizes(),
                Info.Model.GetZCellSizes(), SourceObserve, Info.Model.GetConductivities(),
                Info.Model.GetBackgroundConductivities(),
                Info.Model.GetBackgroundThicknesses(), true);
            //write an empty source file for the second source polarization
            WriteSourceFile((CurrDir / sourcebfilename).string(), HxSourceXIndex,
                HxSourceYIndex, HxSourceDepth, HySourceXIndex, HySourceYIndex,
                HySourceDepth, HxSourceXIndex, HxSourceYIndex, HxSourceDepth, HZeros,
                HZeros, HZeros, Info.Model.GetZCoordinates(), Info.Model.GetZCellSizes(),
                nmodx, nmody);
          }
      }
    //make the sources for the magnetic dipoles
    //now we calculate the response to magnetic dipole sources
    const std::complex<double> omega_mu =
        1.0
            / (std::complex<double>(0.0, jif3D::mag_mu) * 2.0
                * boost::math::constants::pi<double>()
                * Data.GetFrequencies()[Info.freqindex]);

    std::complex<double> ZxxEx, ZxyEx, ZyxEx, ZyyEx;
    std::complex<double> ZxxEy, ZxyEy, ZyxEy, ZyyEy;
    //std::cout << " Freq: " << CurrFreq.at(0) << std::endl;
    for (size_t j = 0; j < nstats; ++j)
      {

        const size_t offset_Ex = (nmodx * nmody)
            * MeasDepthIndices[Data.GetExIndices()[j + ind_shift]]
            + XSourceXIndex.at(2 * j) * nmody + XSourceYIndex.at(2 * j);
        const size_t offset_Ey = (nmodx * nmody)
            * MeasDepthIndices[Data.GetEyIndices()[j + ind_shift]]
            + XSourceXIndex.at(2 * j + 1) * nmody + XSourceYIndex.at(2 * j + 1);
        const size_t offset_Hx = (nmodx * nmody)
            * MeasDepthIndices[Data.GetHxIndices()[j + ind_shift]]
            + HxSourceXIndex.at(j) * nmody + HxSourceYIndex.at(j);
        const size_t offset_Hy = (nmodx * nmody)
            * MeasDepthIndices[Data.GetHyIndices()[j + ind_shift]]
            + HySourceXIndex.at(j) * nmody + HySourceYIndex.at(j);
        //calculate the full impedance at the location of the Ex field measurements
        FieldsToImpedance(Calc->GetEx1()[offset_Ex], Calc->GetEx2()[offset_Ex],
            Calc->GetEy1()[offset_Ex], Calc->GetEy2()[offset_Ex],
            Calc->GetHx1()[offset_Hx], Calc->GetHx2()[offset_Hx],
            Calc->GetHy1()[offset_Hy], Calc->GetHy2()[offset_Hy], ZxxEx, ZxyEx, ZyxEx,
            ZyyEx);
        //calculate the full impedance at the location of the Ey field measurements
        FieldsToImpedance(Calc->GetEx1()[offset_Ey], Calc->GetEx2()[offset_Ey],
            Calc->GetEy1()[offset_Ey], Calc->GetEy2()[offset_Ey],
            Calc->GetHx1()[offset_Hx], Calc->GetHx2()[offset_Hx],
            Calc->GetHy1()[offset_Hy], Calc->GetHy2()[offset_Hy], ZxxEy, ZxyEy, ZyxEy,
            ZyyEy);

        size_t offset = freq_start_index + j * 8;

        //project the electric dipole to magnetic dipole moments using the undistorted impedance
        HXPolMoments1[j] = omega_mu
            * (ZxxEx * EXPolMoments1[2 * j] + ZyxEx * EYPolMoments1[2 * j]
                + ZxxEy * EXPolMoments1[2 * j + 1] + ZyxEy * EYPolMoments1[2 * j + 1]);
        HYPolMoments1[j] = omega_mu
            * (ZxyEx * EXPolMoments1[2 * j] + ZyyEx * EYPolMoments1[2 * j]
                + ZxyEy * EXPolMoments1[2 * j + 1] + ZyyEy * EYPolMoments1[2 * j + 1]);
        HXPolMoments2[j] = omega_mu
            * (ZxxEx * EXPolMoments2[2 * j] + ZyxEx * EYPolMoments2[2 * j]
                + ZxxEy * EXPolMoments2[2 * j + 1] + ZyxEy * EYPolMoments2[2 * j + 1]);
        HYPolMoments2[j] = omega_mu
            * (ZxyEx * EXPolMoments2[2 * j] + ZyyEx * EYPolMoments2[2 * j]
                + ZxyEy * EXPolMoments2[2 * j + 1] + ZyyEy * EYPolMoments2[2 * j + 1]);
        /*std::cout << "ZEx:" << ZxxEx << " " << ZxyEx << " " << ZyxEx << " " << ZyyEx
         << std::endl;
         std::cout << "ZEy:" << ZxxEy << " " << ZxyEy << " " << ZyxEy << " " << ZyyEy
         << std::endl;
         std::cout << "EPol1:" << EXPolMoments1[2 * j] << " " << EXPolMoments1[2 * j + 1]
         << " " << EYPolMoments1[2 * j] << EYPolMoments1[2 * j + 1] << std::endl;
         std::cout << "EPol2:" << EXPolMoments2[2 * j] << " " << EXPolMoments2[2 * j + 1]
         << " " << EYPolMoments2[2 * j] << EYPolMoments2[2 * j + 1] << std::endl;
         std::cout << j << " " << HXPolMoments1[j] << " " << HYPolMoments1[j] << " "
         << HXPolMoments2[j] << HYPolMoments2[j] << std::endl << std::endl;*/
      }

    std::vector<std::complex<double> > Ux1_mag, Ux2_mag, Uy1_mag, Uy2_mag, Uz1_mag,
        Uz2_mag;
    //calculate the first polarization and read the adjoint fields
#pragma omp task default(shared)
    CalcU(MdipName.string() + PolExt[0], HXPolMoments1, HYPolMoments1, Zeros, Ux1_mag,
        Uy1_mag, Uz1_mag, HxSourceXIndex, HxSourceYIndex, HxSourceDepth, HySourceXIndex,
        HySourceYIndex, HySourceDepth, ZeroIndex, ZeroIndex, HxSourceDepth,
        Info.Model.GetZCoordinates(), Info.Model.GetZCellSizes(), nmodx, nmody, nmodz);
    //calculate the second polarization and read the adjoint fields
#pragma omp task default(shared)
    CalcU(MdipName.string() + PolExt[1], HXPolMoments2, HYPolMoments2, Zeros, Ux2_mag,
        Uy2_mag, Uz2_mag, HxSourceXIndex, HxSourceYIndex, HxSourceDepth, HySourceXIndex,
        HySourceYIndex, HySourceDepth, ZeroIndex, ZeroIndex, HySourceDepth,
        Info.Model.GetZCoordinates(), Info.Model.GetZCellSizes(), nmodx, nmody, nmodz);

    const double cell_sizex = Info.Model.GetXCellSizes()[0];
    const double cell_sizey = Info.Model.GetYCellSizes()[0];
    //now we can calculate the gradient for each model cell
    double Volume, gradinc;
#pragma omp taskwait
    for (size_t j = 0; j < nmod; ++j)
      {
        Volume = cell_sizex * cell_sizey * Info.Model.GetZCellSizes()[j % nmodz];
        //this is an implementation of eq. 14 in Avdeev and Avdeeva
        gradinc = std::real(
            (Ux1_el[j] + Ux1_mag[j]) * Ex1_all[j] + (Uy1_el[j] + Uy1_mag[j]) * Ey1_all[j]
                + (Uz1_el[j] + Uz1_mag[j]) * Ez1_all[j]
                + (Ux2_el[j] + Ux2_mag[j]) * Ex2_all[j]
                + (Uy2_el[j] + Uy2_mag[j]) * Ey2_all[j]
                + (Uz2_el[j] + Uz2_mag[j]) * Ez2_all[j]) * Volume;
        Gradient(j) += gradinc;
      }
    return GradResult(Gradient);
  }

GradResult TipperDerivativeFreq(const ForwardInfo &Info, const jif3D::TipperData &Data,
    const jif3D::rvec &Misfit, boost::shared_ptr<jif3D::X3DFieldCalculator> Calc)
  {
    //a few commonly used quantities for shorter notation
    const size_t nmodx = Info.Model.GetConductivities().shape()[0];
    const size_t nmody = Info.Model.GetConductivities().shape()[1];
    const size_t nmodz = Info.Model.GetConductivities().shape()[2];
    //the number of observations in the fields files, one for each cell in the layer
    const size_t nobs = nmodx * nmody;
    //the number of measurement sites
    const size_t nmeas = Data.GetMeasPosX().size();
    //the number of transfer functions
    const size_t nstats = Data.GetHxIndices().size() / Data.GetFrequencies().size();

    const size_t nmod = nmodx * nmody * nmodz;
    std::vector<double> ShiftDepth;
    //for the controlled source calculations we do not actually
    //need any observe layers as we are only interested in the
    //anomalous fields
    std::vector<double> SourceObserve(1, 0.0);
    std::vector<size_t> MeasDepthIndices;
    size_t nlevels = ConstructDepthIndices(MeasDepthIndices, ShiftDepth, Info.Model,
        Data.GetMeasPosZ());
    jif3D::rvec Gradient(nmod, 0.0);

    fs::path TempDir(Info.TempDirName);
    fs::path ForwardDirName = Calc->GetForwardDirName();

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
    const size_t freq_start_index = nstats * Info.freqindex * 4;
    const size_t ind_shift = nstats * Info.freqindex;

    std::vector<std::complex<double> > HXPolMoments1(nstats, 0.0), HXPolMoments2(nstats,
        0.0), HYPolMoments1(nstats, 0.0), HYPolMoments2(nstats, 0.0), HZPolMoments1(
        nstats, 0.0), HZPolMoments2(nstats, 0.0), Zeros(nstats, 0.0);
    std::vector<double> HxSourceDepth(nstats, 0.0), HySourceDepth(nstats, 0.0),
        HzSourceDepth(nstats, 0.0);
    std::vector<size_t> HxSourceXIndex(nstats, 0), HxSourceYIndex(nstats, 0),
        HySourceXIndex(nstats, 0), HySourceYIndex(nstats, 0), HzSourceXIndex(nstats, 0),
        HzSourceYIndex(nstats, 0);

    //we only want to calculate for one frequency
    //so our vector has just 1 element
    std::vector<double> CurrFreq =
      { Data.GetFrequencies()[Info.freqindex] };
    std::vector<std::string> PolExt =
      { "a", "b" };

    //now we calculate the response to magnetic dipole sources
    fs::path MdipPath = TempDir
        / MakeUniqueName(Info.NameRoot, X3DModel::MDIP, Info.freqindex);
    std::string MdipName = MdipPath.string() + "tip";
    //write the files for the magnetic dipole calculation
#pragma omp critical(gradient_writemodel_mdip)
      {
        for (auto Ext : PolExt)
          {
            std::string CurrName = MdipName + Ext;
            fs::path CurrDir = MdipName + Ext + dirext;
            MakeRunFile(CurrName, CurrDir.string(), Info.X3DName);
            WriteProjectFile(CurrDir.string(), CurrFreq, X3DModel::MDIP, resultfilename,
                modelfilename, Info.GreenStage1, Info.GreenStage4);
            Write3DModelForX3D((CurrDir / modelfilename).string(),
                Info.Model.GetXCellSizes(), Info.Model.GetYCellSizes(),
                Info.Model.GetZCellSizes(), SourceObserve, Info.Model.GetConductivities(),
                Info.Model.GetBackgroundConductivities(),
                Info.Model.GetBackgroundThicknesses(), true);
            //write an empty source file for the second source polarization
            WriteSourceFile((CurrDir / sourcebfilename).string(), HxSourceXIndex,
                HxSourceYIndex, HxSourceDepth, HySourceXIndex, HySourceYIndex,
                HySourceDepth, HzSourceXIndex, HzSourceYIndex, HzSourceDepth, Zeros,
                Zeros, Zeros, Info.Model.GetZCoordinates(), Info.Model.GetZCellSizes(),
                nmodx, nmody);
          }
      }
    //make the sources for the magnetic dipoles
    //now we calculate the response to magnetic dipole sources
    const std::complex<double> omega_mu = 1.0
        / (std::complex<double>(0.0, jif3D::mag_mu) * 2.0
            * boost::math::constants::pi<double>()
            * Data.GetFrequencies()[Info.freqindex]);
    for (size_t j = 0; j < nstats; ++j)
      {
        const size_t hxind = Data.GetHxIndices()[j + ind_shift];
        const size_t hyind = Data.GetHyIndices()[j + ind_shift];
        const size_t hzind = Data.GetHzIndices()[j + ind_shift];
        boost::array<ThreeDModelBase::t3DModelData::index, 3> StationHxIndex =
            Info.Model.FindAssociatedIndices(Data.GetMeasPosX()[hxind],
                Data.GetMeasPosY()[hxind], Data.GetMeasPosZ()[hxind]);
        const size_t offset_Hx = (nmodx * nmody) * MeasDepthIndices[hxind]
            + StationHxIndex[0] * nmody + StationHxIndex[1];

        boost::array<ThreeDModelBase::t3DModelData::index, 3> StationHyIndex =
            Info.Model.FindAssociatedIndices(Data.GetMeasPosX()[hyind],
                Data.GetMeasPosY()[hyind], Data.GetMeasPosZ()[hyind]);
        const size_t offset_Hy = (nmodx * nmody) * MeasDepthIndices[hyind]
            + StationHyIndex[0] * nmody + StationHyIndex[1];

        boost::array<ThreeDModelBase::t3DModelData::index, 3> StationHzIndex =
            Info.Model.FindAssociatedIndices(Data.GetMeasPosX()[hzind],
                Data.GetMeasPosY()[hzind], Data.GetMeasPosZ()[hzind]);
        const size_t offset_Hz = (nmodx * nmody) * MeasDepthIndices[hzind]
            + StationHzIndex[0] * nmody + StationHzIndex[1];

        HxSourceXIndex.at(j) = StationHxIndex[0];
        HxSourceYIndex.at(j) = StationHxIndex[1];
        HxSourceDepth.at(j) = Data.GetMeasPosZ()[hxind];

        HySourceXIndex.at(j) = StationHyIndex[0];
        HySourceYIndex.at(j) = StationHyIndex[1];
        HySourceDepth.at(j) = Data.GetMeasPosZ()[hyind];

        HzSourceXIndex.at(j) = StationHzIndex[0];
        HzSourceYIndex.at(j) = StationHzIndex[1];
        HzSourceDepth.at(j) = Data.GetMeasPosZ()[hzind];

        std::complex<double> Hx1 = Calc->GetHx1()[offset_Hx];
        std::complex<double> Hx2 = Calc->GetHx2()[offset_Hx];
        std::complex<double> Hy1 = Calc->GetHy1()[offset_Hy];
        std::complex<double> Hy2 = Calc->GetHy2()[offset_Hy];
        std::complex<double> Hz1 = Calc->GetHz1()[offset_Hz];
        std::complex<double> Hz2 = Calc->GetHz2()[offset_Hz];

        std::complex<double> Tx, Ty;
        FieldsToTipper(Hx1, Hx2, Hy1, Hy2, Hz1, Hz2, Tx, Ty);

        const size_t siteindex = freq_start_index + j * 4;
        const std::complex<double> magdet = 1. / (Hx1 * Hy2 - Hx2 * Hy1);
        const std::complex<double> A00(Misfit(siteindex), -Misfit(siteindex + 1));
        const std::complex<double> A01(Misfit(siteindex + 2), -Misfit(siteindex + 3));

        HZPolMoments1.at(j) = -omega_mu * magdet * (A00 * Hy2 - A01 * Hx2);
        HZPolMoments2.at(j) = -omega_mu * magdet * (-A00 * Hy1 + A01 * Hx1);

        HXPolMoments1.at(j) = -Tx * HZPolMoments1[j];
        HYPolMoments1.at(j) = -Ty * HZPolMoments1[j];
        HXPolMoments2.at(j) = -Tx * HZPolMoments2[j];
        HYPolMoments2.at(j) = -Ty * HZPolMoments2[j];

      }

    std::vector<std::complex<double> > Ux1_mag, Ux2_mag, Uy1_mag, Uy2_mag, Uz1_mag,
        Uz2_mag;
    //calculate the first polarization and read the adjoint fields
#pragma omp task default(shared)
    CalcU(MdipName + PolExt[0], HXPolMoments1, HYPolMoments1, HZPolMoments1, Ux1_mag,
        Uy1_mag, Uz1_mag, HxSourceXIndex, HxSourceYIndex, HxSourceDepth, HySourceXIndex,
        HySourceYIndex, HySourceDepth, HzSourceXIndex, HzSourceYIndex, HzSourceDepth,
        Info.Model.GetZCoordinates(), Info.Model.GetZCellSizes(), nmodx, nmody, nmodz);
    //calculate the second polarization and read the adjoint fields
#pragma omp task default(shared)
    CalcU(MdipName + PolExt[1], HXPolMoments2, HYPolMoments2, HZPolMoments2, Ux2_mag,
        Uy2_mag, Uz2_mag, HxSourceXIndex, HxSourceYIndex, HxSourceDepth, HySourceXIndex,
        HySourceYIndex, HySourceDepth, HzSourceXIndex, HzSourceYIndex, HzSourceDepth,
        Info.Model.GetZCoordinates(), Info.Model.GetZCellSizes(), nmodx, nmody, nmodz);

    const double cell_sizex = Info.Model.GetXCellSizes()[0];
    const double cell_sizey = Info.Model.GetYCellSizes()[0];
    //now we can calculate the gradient for each model cell
    double Volume, gradinc;
#pragma omp taskwait
    for (size_t j = 0; j < nmod; ++j)
      {
        Volume = cell_sizex * cell_sizey * Info.Model.GetZCellSizes()[j % nmodz];
        //this is an implementation of eq. 14 in Avdeev and Avdeeva
        Gradient(j) = std::real(
            Ux1_mag[j] * Ex1_all[j] + Uy1_mag[j] * Ey1_all[j] + Uz1_mag[j] * Ez1_all[j]
                + Ux2_mag[j] * Ex2_all[j] + Uy2_mag[j] * Ey2_all[j]
                + Uz2_mag[j] * Ez2_all[j]) * Volume;
      }
    return GradResult(Gradient);
  }

#ifdef HAVEHPX
HPX_REGISTER_ACTION(CalculateFrequency_action)
HPX_REGISTER_ACTION(LQDerivativeFreq_action)
#endif
