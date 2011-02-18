//============================================================================
// Name        : calccross.cpp
// Author      : Sep 29, 2010
// Version     :
// Copyright   : 2010, mmoorkamp
//============================================================================

#include "../Global/FileUtil.h"
#include "../Regularization/CrossGradient.h"
#include "../Tomo/ThreeDSeismicModel.h"

/*! \file calccross.cpp
 * Calculate the cross-gradient constraint for each model cell between two models.
 * For simplicity the two input models have to be seismic models and we do not
 * consider any model covariances.
 */

int main()
  {
    std::string ModelName1 = jiba::AskFilename("Model file 1: ");
    std::string ModelName2 = jiba::AskFilename("Model file 2: ");

    jiba::ThreeDSeismicModel Model1, Model2;

    Model1.ReadNetCDF(ModelName1);
    Model2.ReadNetCDF(ModelName2);

    const size_t nparam = Model1.GetSlownesses().num_elements();

    if (nparam != Model2.GetSlownesses().num_elements())
      {
        std::cerr << "Model sizes do not match !";
        return 100;
      }
    jiba::rvec ModelVec(2 * nparam);
    std::copy(Model1.GetSlownesses().origin(), Model1.GetSlownesses().origin()
        + nparam, ModelVec.begin());
    std::copy(Model2.GetSlownesses().origin(), Model2.GetSlownesses().origin()
        + nparam, ModelVec.begin() + nparam);

    jiba::CrossGradient CGObjective(Model1);

    CGObjective.CalcMisfit(ModelVec);

    jiba::rvec DiffVec(nparam);
    for (size_t i = 0; i < nparam; ++i)
      {
        DiffVec(i) = std::pow(CGObjective.GetDataDifference()(i), 2)
            + std::pow(CGObjective.GetDataDifference()(i + nparam), 2)
            + std::pow(CGObjective.GetDataDifference()(i + 2 * nparam), 2);
      }

    std::copy(DiffVec.begin(), DiffVec.end(), Model1.SetSlownesses().origin());

    std::string outfilename = jiba::AskFilename("Outfile name: ", false);
    Model1.WriteNetCDF(outfilename);

  }

