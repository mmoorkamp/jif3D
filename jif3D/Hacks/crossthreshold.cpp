//============================================================================
// Name        : calccross.cpp
// Author      : Sep 29, 2010
// Version     :
// Copyright   : 2010, mmoorkamp
//============================================================================

#include "../Global/FileUtil.h"
#include "../Global/NumUtil.h"
#include "../ModelBase/VTKTools.h"
#include "../ModelBase/ThreeDModelBase.h"
#include "../Gravity/ThreeDGravityModel.h"

/*! \file calccross.cpp
 * Calculate the cross-gradient constraint for each model cell between two models.
 * For simplicity the two input models have to be seismic models and we do not
 * consider any model covariances.
 */

int main()
  {
    std::string ModelName1 = jif3D::AskFilename("Model file 1: ");
    jif3D::ThreeDGravityModel Model;

    std::vector<double> XCoord, YCoord, ZCoord;
    jif3D::ThreeDModelBase::t3DModelData CompX, CompY, CompZ;
    jif3D::Read3DVectorModelFromVTK(ModelName1, XCoord, YCoord, ZCoord, CompX, CompY,
        CompZ);
    jif3D::ThreeDModelBase::t3DModelData Abs(CompX), Tear(CompX);
    Model.SetMeshSize(XCoord.size() - 1, YCoord.size() - 1, ZCoord.size() - 1);
    Model.SetXCoordinates(XCoord);
    Model.SetYCoordinates(YCoord);
    Model.SetZCoordinates(ZCoord);

    double threshold = 0.0;
    std::cout << "Threshold: ";
    std::cin >> threshold;
    double sum = 0;
    for (size_t i = 0; i < Model.GetModelShape()[0]; ++i)
      {
        for (size_t j = 0; j < Model.GetModelShape()[1]; ++j)
          {
            for (size_t k = 0; k < Model.GetModelShape()[2]; ++k)
              {
                Abs[i][j][k] = jif3D::pow2(CompX[i][j][k]) + jif3D::pow2(CompY[i][j][k])
                    + jif3D::pow2(CompZ[i][j][k]);
                sum += Abs[i][j][k];
                if (Abs[i][j][k] > threshold)
                  {
                    Tear[i][j][k] = 0.0;
                  }
                else
                  {
                    Tear[i][j][k] = 1.0;
                  }
              }
          }
      }
    std::cout << CompX[0][0][0] << " " << CompY[0][0][0] << " " << CompZ[0][0][0] << std::endl;
    std::cout << "CG: " << sum << std::endl;
    Model.SetData() = Abs;
    Model.WriteNetCDF(ModelName1 + ".cross.nc");
    Model.WriteVTK(ModelName1 + ".cross.vtk");
    Model.SetData() = Tear;

    Model.WriteNetCDF(ModelName1 + ".ctear.nc");
    Model.WriteVTK(ModelName1 + ".ctear.vtk");
  }

