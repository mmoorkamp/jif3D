//============================================================================
// Name        : GravityObjective.cpp
// Author      : Apr 16, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#include "GravityObjective.h"
#include "MinMemGravityCalculator.h"
#include <cassert>
namespace jiba
  {
	void VectorToModel(const jiba::rvec &InVector, jiba::ThreeDGravityModel &Model)
	{
		const size_t ngrid = Model.GetDensities().num_elements();
	assert(InVector.size() ==  ngrid ||
		   InVector.size() == ngrid+ Model.GetBackgroundDensities().size());
	   std::copy(InVector.begin(), InVector.begin()+ngrid,
		            Model.SetDensities().origin());
	   if (InVector.size() == ngrid+ Model.GetBackgroundDensities().size())
	   {
	     std::vector<double> Background(InVector.begin()+ngrid,InVector.end());
         Model.SetBackgroundDensities(Background);
	   }
	}

    GravityObjective::GravityObjective(bool ftg, bool cuda)
      {
        if (ftg)
          {
            Calculator
                = boost::shared_ptr<jiba::ThreeDGravityCalculator>(
                    jiba::CreateGravityCalculator<jiba::MinMemGravityCalculator>::MakeTensor(
                        cuda));
          }
        else
          {
            Calculator
                = boost::shared_ptr<jiba::ThreeDGravityCalculator>(
                    jiba::CreateGravityCalculator<jiba::MinMemGravityCalculator>::MakeScalar(
                        cuda));
          }
      }

    GravityObjective::~GravityObjective()
      {

      }

    void GravityObjective::ImplDataDifference(const jiba::rvec &Model,
        jiba::rvec &Diff)
      {
        assert(DensityModel.GetMeasPosX().size() * Calculator->GetDataPerMeasurement() == ObservedData.size() );
        VectorToModel(Model,DensityModel);
        jiba::rvec SynthData(Calculator->Calculate(DensityModel));
        Diff.resize(ObservedData.size());
        std::transform(SynthData.begin(), SynthData.end(),
            ObservedData.begin(), Diff.begin(), std::minus<double>());
      }

    jiba::rvec GravityObjective::ImplGradient(const jiba::rvec &Model,
        const jiba::rvec &Diff)
      {
        assert(DensityModel.GetMeasPosX().size()* Calculator->GetDataPerMeasurement() == Diff.size() );
        VectorToModel(Model,DensityModel);
        return Calculator->LQDerivative(DensityModel, Diff);
      }
  }
