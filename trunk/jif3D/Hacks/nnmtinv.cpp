//============================================================================
// Name        : onedmtinv.cpp
// Author      : 13 Jun 2012
// Version     : 
// Copyright   : 2012, mm489
//============================================================================

#include <iostream>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/program_options.hpp>
#include "../Global/FileUtil.h"
#include "../Global/convert.h"
#include "../Global/Noise.h"
#include "../Global/VecMat.h"
#include "../MT/X3DModel.h"
#include "../MT/ReadWriteImpedances.h"
#include "../MT/OneDMTObjective.h"
#include "../MT/MTEquations.h"
#include "../Inversion/LimitedMemoryQuasiNewton.h"
#include "../Inversion/JointObjective.h"
#include "../Inversion/ModelTransforms.h"
#include "../Regularization/OneDRegularization.h"
#include "../ModelBase/VTKTools.h"

namespace po = boost::program_options;

int FindNearest(double modx, double mody, const std::vector<double> &measx,
    const std::vector<double> &measy)
  {
    const size_t nstat = measx.size();
    double mindist = 1e55;
    int index = 0;
    for (size_t i = 0; i < nstat; ++i)
      {
        double dist = std::pow(modx - measx[i], 2) + std::pow(mody - measy[i], 2);
        if (dist < mindist)
          {
            mindist = dist;
            index = i;
          }
      }
    return index;
  }

int main(int argc, char *argv[])
  {
    double mincond = 1e-6;
    double maxcond = 10;
    po::options_description desc("General options");
    desc.add_options()("help", "produce help message")("threads", po::value<int>(),
        "The number of openmp threads")("mincond",
        po::value(&mincond)->default_value(1e-6),
        "The minimum value for conductivity in S/m")("maxcond",
        po::value(&maxcond)->default_value(10),
        "The maximum value for conductivity in S/m");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    std::string modelfilename = jif3D::AskFilename("Model file: ");
    std::string datafilename = jif3D::AskFilename("Data file: ");
    std::string dataextension = jif3D::GetFileExtension(modelfilename);
    jif3D::X3DModel Model;
    Model.ReadNetCDF(modelfilename);
    jif3D::rvec Impedances, ImpError;
    std::vector<double> Frequencies, StatXCoord, StatYCoord, StatZCoord, GoodXCoord,
        GoodYCoord;

    std::vector<double> C;
    jif3D::ReadImpedancesFromNetCDF(datafilename, Frequencies, StatXCoord, StatYCoord,
        StatZCoord, Impedances, ImpError, C);

    double reglambda = 1.0;
    std::cout << "Lambda: ";
    std::cin >> reglambda;

    const size_t nfreq = Frequencies.size();
    const size_t nstat = StatXCoord.size();
    const size_t nlayers = Model.GetZCellSizes().size();

    std::ofstream misfitfile("misfit.out");
    std::ofstream modelfile("model.out");
    std::ofstream impfile("imp.out");

    std::vector<std::vector<double> > InvResult;

    for (size_t i = 0; i < nstat; ++i)
      {
        boost::shared_ptr<jif3D::OneDMTObjective> MTObjective(new jif3D::OneDMTObjective);
        jif3D::rvec OneDImp(nfreq * 2, 0.0), OneDErr(nfreq * 2, 0.0);

        for (size_t j = 0; j < nfreq; ++j)
          {

            size_t siteindex = (j * nstat + i) * 8;
            OneDImp(j * 2) = (Impedances(siteindex + 2) - Impedances(siteindex + 4)) / 2.0;
            OneDImp(j * 2 + 1) = (Impedances(siteindex + 3) - Impedances(siteindex + 5))/2.0;

            OneDErr(j * 2) = std::max(ImpError(siteindex + 2), ImpError(siteindex + 4));
            OneDErr(j * 2 + 1) = OneDErr(j * 2);
          }

        Model.SetFrequencies() = Frequencies;
        MTObjective->SetModelGeometry(Model);
        MTObjective->SetObservedData(OneDImp);
        jif3D::rvec AvgImp(OneDImp.size());

        for (size_t i = 0; i < nfreq; ++i)
          {
            AvgImp(2 * i) = (OneDImp(2 * i) + OneDImp(2 * i + 1)) / 2.0;
            AvgImp(2 * i + 1) = AvgImp(2 * i);
          }
        impfile << i << " ";
        std::copy(AvgImp.begin(), AvgImp.end(),
            std::ostream_iterator<double>(impfile, " "));
        impfile << std::endl;
        jif3D::rvec DataError(jif3D::ConstructError(AvgImp, OneDErr, 0.05, 0.0));
        MTObjective->SetDataError(DataError);

        jif3D::rvec InvModel(nlayers);
        std::copy(Model.GetBackgroundConductivities().begin(),
            Model.GetBackgroundConductivities().end(), InvModel.begin());
        boost::shared_ptr<jif3D::OneDRegularization> Regularization(
            new jif3D::OneDRegularization(nlayers));
        //Regularization->SetReferenceModel(InvModel);

        boost::shared_ptr<jif3D::ChainedTransform> ConductivityTransform(
            new jif3D::ChainedTransform);

        jif3D::rvec RefModel(InvModel);
        std::fill(RefModel.begin(), RefModel.end(), 1.0);
        //because the tanh transform is used inside a logarithmic transform
        //we need to take the natural logarithm of the actual minimum and maximum
        ConductivityTransform->AppendTransform(
            boost::shared_ptr<jif3D::GeneralModelTransform>(
                new jif3D::TanhTransform(std::log(mincond), std::log(maxcond))));
        ConductivityTransform->AppendTransform(
            boost::shared_ptr<jif3D::GeneralModelTransform>(
                new jif3D::LogTransform(RefModel)));
        // ConductivityTransform->AddTransform(
        //    boost::shared_ptr<jif3D::GeneralModelTransform>(new jif3D::ModelCopyTransform));
        InvModel = ConductivityTransform->PhysicalToGeneralized(InvModel);

        boost::shared_ptr<jif3D::ModelCopyTransform> RegTrans(
            new jif3D::ModelCopyTransform);
        boost::shared_ptr<jif3D::JointObjective> Objective(
            new jif3D::JointObjective(false));
        Objective->AddObjective(MTObjective, ConductivityTransform, 1.0, "MT",jif3D::JointObjective::datafit);
        Objective->AddObjective(Regularization, RegTrans, reglambda, "Reg");

        jif3D::LimitedMemoryQuasiNewton Optimizer(Objective);
        //Optimizer.SetModelCovDiag(InvModel);
        Objective->CalcMisfit(InvModel);
        size_t iteration = 0;
        double maxiter = 100.0;
        const size_t ndata = OneDImp.size();
        bool terminate = false;
        while (iteration < maxiter && !terminate
            && Objective->GetIndividualFits()[0] > ndata)
          {
            try
              {
                jif3D::rvec CondInvModel = ConductivityTransform->GeneralizedToPhysical(
                    InvModel);
                std::vector<double> Tmp(CondInvModel.begin(), CondInvModel.end());
                //Model.SetBackgroundConductivities(Tmp);
                Optimizer.MakeStep(InvModel);
                /*Model.WriteNetCDF(
                 modelfilename + jif3D::stringify(iteration) + ".mt.inv.nc");
                 misfitfile << iteration + 1 << " " << Optimizer.GetMisfit() << " ";
                 std::copy(Objective->GetIndividualFits().begin(),
                 Objective->GetIndividualFits().end(),
                 std::ostream_iterator<double>(misfitfile, " "));

                 misfitfile << std::endl;*/
                iteration++;
                terminate = jif3D::WantAbort();
              } catch (jif3D::FatalException &e)
              {
                std::cerr << e.what() << std::endl;
                terminate = true;
              }
          }

        double misfit = Objective->CalcMisfit(InvModel);
        misfitfile << i << " " << misfit << " ";

        std::copy(Objective->GetIndividualFits().begin(),
            Objective->GetIndividualFits().end(),
            std::ostream_iterator<double>(misfitfile, " "));
        misfitfile << " " << Objective->GetNEval() << " " << Objective->GetRMS()[0];
        misfitfile << std::endl;

        if (sqrt(Objective->GetRMS()[0]) < 3.0)
          {
            InvModel = ConductivityTransform->GeneralizedToPhysical(InvModel);
            GoodXCoord.push_back(StatXCoord[i]);
            GoodYCoord.push_back(StatYCoord[i]);
            std::vector<double> Result(InvModel.begin(),InvModel.end());
            InvResult.push_back(Result);

            modelfile << i << " ";
            std::copy(InvModel.begin(), InvModel.end(),
                std::ostream_iterator<double>(modelfile, " "));
            modelfile << std::endl;

          }

      }

    const double xsize2 = Model.GetXCellSizes()[0] / 2.0;
    const double ysize2 = Model.GetYCellSizes()[0] / 2.0;

    for (size_t i = 0; i < Model.GetXCellSizes().size(); ++i)
      for (size_t j = 0; j < Model.GetYCellSizes().size(); ++j)
        {
          int index = FindNearest(Model.GetXCoordinates()[i] + xsize2,
              Model.GetYCoordinates()[j] + ysize2, GoodXCoord, GoodYCoord);
          for (size_t k = 0; k < Model.GetZCellSizes().size(); ++k)
            {
              Model.SetData()[i][j][k] = InvResult.at(index).at(k);
            }
        }

    //calculate the predicted data
    /*    std::cout << "Calculating response of inversion model." << std::endl;
     jif3D::rvec InvData(jif3D::OneDMTCalculator().Calculate(Model));
     jif3D::rvec FullImp(nfreq * 8, 0.0), FullErr(nfreq * 8, 0.0);
     for (size_t i = 0; i < nfreq; ++i)
     {
     FullImp(i * 8 + 2) = InvData(i * 2);
     FullImp(i * 8 + 3) = InvData(i * 2 + 1);
     FullImp(i * 8 + 4) = -InvData(i * 2);
     FullImp(i * 8 + 5) = -InvData(i * 2 + 1);
     FullErr(i * 8 + 2) = DataError(i * 2);
     FullErr(i * 8 + 3) = DataError(i * 2);
     FullErr(i * 8 + 4) = DataError(i * 2);
     FullErr(i * 8 + 5) = DataError(i * 2);
     }
     std::vector<double> Coord =
     { 0.0 };
     jif3D::WriteImpedancesToNetCDF(modelfilename + ".inv_imp.nc", Frequencies, Coord,
     Coord, Coord, FullImp, FullErr);*/

    //and write out the data and model
    //here we have to distinguish again between scalar and ftg data
    std::cout << "Writing out inversion results." << std::endl;

    jif3D::rvec si(StatXCoord.size());
    std::iota(si.begin(), si.end(), 1);
    jif3D::Write3DDataToVTK(datafilename + ".vtk", "MTSites", si, StatXCoord, StatYCoord,
        StatZCoord);

    Model.WriteVTK(modelfilename + ".inv.vtk");
    Model.WriteNetCDF(modelfilename + ".inv.nc");
    std::cout << std::endl;
  }
