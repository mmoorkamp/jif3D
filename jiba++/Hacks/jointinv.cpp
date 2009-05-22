//============================================================================
// Name        : jointinv.cpp
// Author      : May 12, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


//============================================================================
// Name        : gravinv.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================


#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <boost/bind.hpp>
#include "../Global/convert.h"
#include "../Global/FatalException.h"
#include "../Global/NumUtil.h"
#include "../Global/VectorTransform.h"
#include "../Global/FileUtil.h"
#include "../ModelBase/VTKTools.h"
#include "../ModelBase/NetCDFTools.h"
#include "../Inversion/LimitedMemoryQuasiNewton.h"
#include "../Inversion/NonLinearConjugateGradient.h"
#include "../Inversion/JointObjective.h"
#include "../Inversion/MinDiffRegularization.h"
#include "../Inversion/ModelTransforms.h"
#include "../Inversion/ConstructError.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include "../Tomo/ReadWriteTomographyData.h"
#include "../Tomo/TomographyObjective.h"
#include "../Tomo/TomographyCalculator.h"
#include "../Gravity/GravityObjective.h"
#include "../Gravity/ReadWriteGravityData.h"
#include "../Gravity/ThreeDGravityCalculator.h"
#include "../Gravity/MinMemGravityCalculator.h"
#include "../Gravity/DepthWeighting.h"

namespace ublas = boost::numeric::ublas;

int main(int argc, char *argv[])
{
	//these objects hold information about the measurements and their geometry
	jiba::rvec TomoData, ScalGravData, FTGData;

	//first we read in the starting model and the measured data
	std::string modelfilename = jiba::AskFilename("Starting model Filename: ");
	//we read in the starting modelfile
	jiba::ThreeDSeismicModel TomoModel;
	TomoModel.ReadNetCDF(modelfilename);
	TomoModel.WriteVTK(modelfilename + ".vtk");
	//get the name of the file containing the data and read it in
	std::string tomodatafilename = jiba::AskFilename(
			"Tomography Data Filename: ");

	//read in data
	jiba::ReadTraveltimes(tomodatafilename, TomoData, TomoModel);

	std::string scalgravdatafilename = jiba::AskFilename(
			"Scalar Gravity Data Filename: ");
	std::string ftgdatafilename = jiba::AskFilename("FTG Data Filename: ");
	std::string gravmodelfilename = jiba::AskFilename(
			"Gravity Model Filename: ");
	jiba::ThreeDGravityModel GravModel;
	GravModel.ReadNetCDF(gravmodelfilename);
	GravModel = TomoModel;

	jiba::ThreeDGravityModel::tMeasPosVec PosX, PosY, PosZ;
	jiba::ReadScalarGravityMeasurements(scalgravdatafilename, ScalGravData,
			PosX, PosY, PosZ);
	jiba::ReadTensorGravityMeasurements(ftgdatafilename, FTGData, PosX, PosY,
			PosZ);
	GravModel.ClearMeasurementPoints();
	for (size_t i = 0; i < PosX.size(); ++i)
	{
		GravModel.AddMeasurementPoint(PosX.at(i), PosY.at(i), PosZ.at(i));
	}

	//if we don't have data inversion doesn't make sense;
	if (TomoData.empty() || ScalGravData.empty())
	{
		std::cerr << "No measurements defined" << std::endl;
		exit(100);
	}

	jiba::rvec InvModel(TomoModel.GetSlownesses().num_elements());
	std::copy(TomoModel.GetSlownesses().origin(),
			TomoModel.GetSlownesses().origin()
					+ TomoModel.GetSlownesses().num_elements(),
			InvModel.begin());

	jiba::rvec RefModel(InvModel);

	boost::shared_ptr<jiba::GeneralModelTransform> DensityTransform(
			new jiba::LogDensityTransform(RefModel));
	boost::shared_ptr<jiba::GeneralModelTransform> SlownessTransform(
			new jiba::LogTransform(RefModel));
	//double average = std::accumulate(InvModel.begin(),InvModel.end(),0.0)/InvModel.size();
	//std::fill(RefModel.begin(),RefModel.end(),average);
	InvModel = SlownessTransform->PhysicalToGeneralized(InvModel);
	jiba::rvec
			DensStartModel(DensityTransform->GeneralizedToPhysical(InvModel));
	std::cout << "Background layers: "
			<< GravModel.GetBackgroundDensities().size() << std::endl;
	std::copy(DensStartModel.begin(), DensStartModel.end(),
			GravModel.SetDensities().origin());
	GravModel.WriteNetCDF("out_dens.nc");

	boost::shared_ptr<jiba::TomographyObjective> TomoObjective(
			new jiba::TomographyObjective());
	TomoObjective->SetObservedData(TomoData);
	TomoObjective->SetModelGeometry(TomoModel);
	TomoObjective->SetDataCovar(jiba::ConstructError(TomoData, sqrt(0.02)));

	boost::shared_ptr<jiba::GravityObjective> ScalGravObjective(
			new jiba::GravityObjective());
	ScalGravObjective->SetObservedData(ScalGravData);
	ScalGravObjective->SetModelGeometry(GravModel);
	ScalGravObjective->SetDataCovar(jiba::ConstructError(ScalGravData, sqrt(
			0.02)));

	boost::shared_ptr<jiba::GravityObjective> FTGObjective(
			new jiba::GravityObjective(true));
	FTGObjective->SetObservedData(FTGData);
	FTGObjective->SetModelGeometry(GravModel);
	FTGObjective->SetDataCovar(jiba::ConstructError(FTGData, sqrt(0.02)));

	const double z0 = 5.0;
	const double DepthExponent = -2.0;
	jiba::rvec WeightVector, ModelWeight(InvModel.size());
	//calculate the depth scaling
	jiba::ConstructDepthWeighting(GravModel.GetZCellSizes(), z0, WeightVector,
			jiba::WeightingTerm(DepthExponent));
	for (size_t i = 0; i < ModelWeight.size(); ++i)
	{
		ModelWeight(i) = WeightVector(i % GravModel.GetZCellSizes().size());
	}

	boost::shared_ptr<jiba::JointObjective> Objective(
			new jiba::JointObjective());
	boost::shared_ptr<jiba::MinDiffRegularization> Regularization(
			new jiba::MinDiffRegularization());

	Regularization->SetReferenceModel(RefModel);
	//Regularization->SetDataCovar(RefModel);
	double tomolambda = 1.0;
	double scalgravlambda = 1.0;
	double ftglambda = 1.0;
	double reglambda = 1.0;
	std::cout << "Tomography Lambda: ";
	std::cin >> tomolambda;
	std::cout << "Scalar Gravimetry Lambda: ";
	std::cin >> scalgravlambda;
	std::cout << "FTG Lambda: ";
	std::cin >> ftglambda;
	std::cout << "Regularization Lambda: ";
	std::cin >> reglambda;
	if (tomolambda > 0.0)
	{
		Objective->AddObjective(TomoObjective, SlownessTransform, tomolambda);
	}
	if (scalgravlambda > 0.0)
	{
		Objective->AddObjective(ScalGravObjective, DensityTransform,
				scalgravlambda);
	}
	if (ftglambda > 0.0)
	{
		Objective->AddObjective(FTGObjective, DensityTransform, ftglambda);
	}
	Objective->AddObjective(Regularization, SlownessTransform, reglambda);

	std::cout << "Performing inversion." << std::endl;

	jiba::LimitedMemoryQuasiNewton LBFGS(Objective, 5);
	//LBFGS.SetModelCovDiag(ModelWeight);

	const size_t ndata = TomoData.size() + ScalGravData.size() + FTGData.size();
	size_t iteration = 0;
	size_t maxiter = 30;
	jiba::rvec TomoInvModel(SlownessTransform->GeneralizedToPhysical(InvModel));
	do
	{
		std::cout << "Iteration" << iteration << std::endl;
		LBFGS.MakeStep(InvModel);

		++iteration;
		TomoInvModel = SlownessTransform->GeneralizedToPhysical(InvModel);
		std::copy(TomoInvModel.begin(), TomoInvModel.begin()
				+ TomoModel.GetSlownesses().num_elements(),
				TomoModel.SetSlownesses().origin());
		std::cout << "Gradient Norm: " << LBFGS.GetGradNorm() << std::endl;
		TomoModel.WriteVTK(modelfilename + jiba::stringify(iteration)
				+ ".tomo.inv.vtk");
		std::cout << "Currrent Misfit: " << LBFGS.GetMisfit() << std::endl;
		std::cout << "Currrent Gradient: " << LBFGS.GetGradNorm() << std::endl;
	} while (iteration < maxiter && LBFGS.GetMisfit() > ndata
			&& LBFGS.GetGradNorm() > 1e-6);

	jiba::rvec DensInvModel(DensityTransform->GeneralizedToPhysical(InvModel));

	std::copy(TomoInvModel.begin(), TomoInvModel.begin()
			+ TomoModel.GetSlownesses().num_elements(),
			TomoModel.SetSlownesses().origin());
	std::copy(DensInvModel.begin(), DensInvModel.begin()
			+ GravModel.SetDensities().num_elements(),
			GravModel.SetDensities().origin());

	//calculate the predicted refraction data
	std::cout << "Calculating response of inversion model." << std::endl;
	jiba::rvec TomoInvData(jiba::TomographyCalculator().Calculate(TomoModel));
	jiba::SaveTraveltimes(modelfilename + ".inv_tt.nc", TomoInvData, TomoModel);

	boost::shared_ptr<jiba::MinMemGravityCalculator>
			GravityCalculator =
					boost::shared_ptr<jiba::MinMemGravityCalculator>(
							jiba::CreateGravityCalculator<
									jiba::MinMemGravityCalculator>::MakeScalar());
	jiba::rvec GravInvData(GravityCalculator->Calculate(GravModel));
	jiba::SaveScalarGravityMeasurements(modelfilename + ".inv_sgd.nc",
			GravInvData, GravModel.GetMeasPosX(), GravModel.GetMeasPosY(),
			GravModel.GetMeasPosZ());
	//and write out the data and model
	//here we have to distinguish again between scalar and ftg data
	std::cout << "Writing out inversion results." << std::endl;

	TomoModel.WriteVTK(modelfilename + ".tomo.inv.vtk");
	GravModel.WriteVTK(modelfilename + ".grav.inv.vtk");
	GravModel.WriteNetCDF(modelfilename + ".grav.inv.nc");
	std::cout << std::endl;
}
