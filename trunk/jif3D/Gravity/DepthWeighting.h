//============================================================================
// Name        : DepthWeighting.h
// Author      : Sep 18, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================

#ifndef DEPTHWEIGHTING_H_
#define DEPTHWEIGHTING_H_

#include "ThreeDGravityModel.h"

namespace jif3D
  {
    /** \addtogroup gravity Gravity forward modeling, display and inversion */
    /* @{ */
    //! This class implements the depth weighting function necessary to counter-act the decay of the sensitivities for gravity data
    /*! Both scalar and tensor gravity data have sensitivity kernels that strongly decay with depth. The result
     * is that in an inversion all structure is concentrated at the top. This class implements functions of the form
     * \f$ (z-z_0)^n \f$ and its derivatives. For scalar data we choose \f$n=-2\f$ and fit \f$ z_0\f$ so we match
     * the decay of the kernel with depth.
     */
    class WeightingTerm
      {
    private:
      //! The exponent in the depth weighting, should match the decay of the sensitivities
      double n;
    public:
      //! Given the exponent specified in the constructor calculate the function value for given z and z0
      double operator()(const double z, const double z0) const;
      //! Calculate the derivative of the function for a given z and z0
      double deriv(const double z, const double z0) const;
      //! Calculate the average function value between z1 and z2 for a given z0
      double average(const double z1, const double z2, const double z0) const;
      //! The constructor takes the exponent n
      WeightingTerm(const double exponent) :
          n(exponent)
        {
        }
      };

    //! Given the values of the sensitivity kernel with depth, find z0 that matches the decay
    double FitZ0(const jif3D::rvec &SensProfile,
        const ThreeDModelBase::t3DModelDim &ZSizes,
        const jif3D::WeightingTerm &WeightFunction);

    //! Given a z0 and the model geometry, construct a vector of weights
    void ConstructDepthWeighting(const ThreeDModelBase::t3DModelDim &ZSizes,
        const double z0, rvec &WeightVector, const jif3D::WeightingTerm &WeightFunction);

    //! Extract sensitivities for a site that is closest to the middle of the modeling domain
    /*! For depth weighting we need the sensitivities below a site to match with our weighting function.
     * We somewhat arbitrarily chose the site closest to the middle, as here the effects of the finite modeling domain
     * should be smallest. The function extracts the row from the matrix that corresponds to this location
     * @param Model The model object, needed for the grid information
     * @param Sensitivities The sensitivity matrix that we want to extract the profile from
     * @param MeasPerPos How many measurements per site,e.g. 1 for only scalar and 9 for only FTG
     * @param SensProfile The depth profile of the sensitivity below the site
     */
    template<class ModelType>
    void ExtractMiddleSens(const ModelType &Model, const jif3D::rmat &Sensitivities,
        const size_t MeasPerPos, jif3D::rvec &SensProfile)
      {
        const size_t nmeas = Model.GetMeasPosX().size();
        const double midx = Model.GetXCoordinates()[Model.GetXCoordinates().size() - 1]
            / 2.0;
        const double midy = Model.GetYCoordinates()[Model.GetYCoordinates().size() - 1]
            / 2.0;

        jif3D::rvec distances(nmeas);
        for (size_t i = 0; i < nmeas; ++i)
          {
            distances(i) = sqrt(
                pow(Model.GetMeasPosX()[i] - midx, 2)
                    + pow(Model.GetMeasPosY()[i] - midy, 2));
          }
        const size_t midindex = distance(distances.begin(),
            std::min_element(distances.begin(), distances.end()));
        boost::array<jif3D::ThreeDModelBase::t3DModelData::index, 3> modelindex(
            Model.FindAssociatedIndices(Model.GetMeasPosX()[midindex],
                Model.GetMeasPosY()[midindex], 0.0));
        //we store the sensitivities for the background at the end of the matrix
        //so we can ignore it here
        boost::numeric::ublas::matrix_row<const jif3D::rmat> MiddleSens(Sensitivities,
            midindex * MeasPerPos);

        const size_t ysize = Model.GetModelShape()[1];
        const size_t zsize = Model.GetModelShape()[2];

        SensProfile.resize(zsize);
        //the same here, if we operate on the first ngrid elements
        //the background does not matter
        const size_t startindex = (zsize * ysize) * modelindex[0] + zsize * modelindex[1];
        std::copy(MiddleSens.begin() + startindex,
            MiddleSens.begin() + startindex + zsize, SensProfile.begin());

      }

    template<class ModelType, class CalculatorType>
    void CalculateMiddleSens(const ModelType &Model, CalculatorType Calculator,
        jif3D::rvec &SensProfile)
      {
        ModelType LocalModel(Model);
        const size_t nx = LocalModel.GetXCellSizes().size();
        const size_t ny = LocalModel.GetYCellSizes().size();
        const double sizex = LocalModel.GetXCoordinates()[nx - 1]
            + LocalModel.GetXCellSizes()[nx - 1];
        const double sizey = LocalModel.GetYCoordinates()[ny - 1]
            + LocalModel.GetYCellSizes()[ny - 1];
        double zpos = -0.01;
        if (LocalModel.GetMeasPosZ().size() > 0)
          {
            zpos = LocalModel.GetMeasPosZ()[0];
          }
        LocalModel.ClearMeasurementPoints();
        LocalModel.AddMeasurementPoint(sizex / 2.0, sizey / 2.0, zpos);

        Calculator.Calculate(LocalModel);
        rvec Misfit(Calculator.GetDataPerMeasurement(), 1.0);
        rvec Deriv = Calculator.LQDerivative(LocalModel, Misfit);
        rmat Sens(1.0, Deriv.size());
        std::copy(Deriv.begin(),Deriv.end(),Sens.data().begin());
        ExtractMiddleSens(LocalModel,Sens,1,SensProfile);

        std::copy(Sens.data().begin(), Sens.data().end(), LocalModel.SetData().origin());
        LocalModel.WriteVTK("sens.vtk");

      }
  /* @} */
  }
#endif /* DEPTHWEIGHTING_H_ */
