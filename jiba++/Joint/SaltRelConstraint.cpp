//============================================================================
// Name        : SaltRelConstraint.cpp
// Author      : May 3, 2011
// Version     : 
// Copyright   : 2011, mmoorkamp
//============================================================================

#include "SaltRelConstraint.h"
#include <fstream>
#include <boost/math/special_functions/pow.hpp>

using boost::math::pow;

namespace jiba
  {

    const double refslow = 1.0 / 3000.0;
    const double refdens = 2.0;
    const double refcond = 0.1;
    SaltRelConstraint::SaltRelConstraint(
        boost::shared_ptr<jiba::GeneralModelTransform> DensTrans,
        boost::shared_ptr<jiba::GeneralModelTransform> CondTrans) :
        DensityTransform(DensTrans), ConductivityTransform(CondTrans)
      {

      }

    SaltRelConstraint::~SaltRelConstraint()
      {
      }

    double SaltTerm(const jiba::rvec &Model, const size_t index)
      {
        const size_t ncells = Model.size() / 3;
        return pow<2>((Model(index) - 1.0 / saltvel) / refslow)
            + pow<2>((Model(index + ncells) - saltdens) / refdens)
            + pow<2>((Model(index + 2 * ncells) - 1.0 / saltres) / refcond);
      }

    void SaltRelConstraint::ImplDataDifference(const jiba::rvec &Model, jiba::rvec &Diff)
      {
        const size_t ncells = Model.size() / 3;
        Diff.resize(ncells);
        //predict densities from slownesses, the first ncells vector entries contain the slownesses
        jiba::rvec PredDens(
            DensityTransform->GeneralizedToPhysical(ublas::subrange(Model, 0, ncells)));
        jiba::rvec PredCond(
            ConductivityTransform->GeneralizedToPhysical(
                ublas::subrange(Model, 0, ncells)));
        std::ofstream condrel("condrel.out");
        std::ofstream gravrel("gravrel.out");
        for (size_t i = 0; i < ncells; ++i)
          {
            int xi, yi, zi;
            Geometry.OffsetToIndex(i, xi, yi, zi);
            if (Geometry.GetSlownesses()[xi][yi][zi] > 0.0)
              {
                //we assume the order slowness, density, conductivity in the model vector
                Diff(i) = pow<2>((PredDens(i) - Model(i + ncells)) / refdens)
                    + pow<2>((PredCond(i) - Model(i + 2 * ncells)) / refcond);
                condrel << 1.0 / Model(i) << " " << 1.0 / Model(i + 2 * ncells) << " "
                    << 1.0 / PredCond(i) << " " << Diff(i) << "\n";
                gravrel << 1.0 / Model(i) << " " << Model(i + ncells) << " "
                    << PredDens(i) << "\n";
                // std::cout << Model(i) << " " << Model(i + ncells) << " "
                //    << Model(i + 2 * ncells) << " " << PredDens(i) << " "
                //    << PredCond(i);
                Diff(i) *= SaltTerm(Model, i);
                //std::cout << " " << SaltTerm(Model, i) << std::endl;
                //we take the square root, as the base class squares each entry of the difference
                //vector before summing up
                Diff(i) = std::sqrt(Diff(i));
              }
            else
              {
                Diff(i) = 0.0;
              }
          }
      }

    jiba::rvec SaltRelConstraint::ImplGradient(const jiba::rvec &Model,
        const jiba::rvec &Diff)
      {
        const size_t nparam = Model.size();
        const size_t ncells = Model.size() / 3;
        jiba::rvec Gradient(nparam, 0.0);
        jiba::rvec DensDiff(ncells, 0.0), CondDiff(ncells, 0.0), SaltDiff(ncells, 0.0);
        //predict densities from slownesses, the first ncells vector entries contain the slownesses
        jiba::rvec PredDens(
            DensityTransform->GeneralizedToPhysical(ublas::subrange(Model, 0, ncells)));
        jiba::rvec PredCond(
            ConductivityTransform->GeneralizedToPhysical(
                ublas::subrange(Model, 0, ncells)));
        for (size_t i = 0; i < ncells; ++i)
          {
            int xi, yi, zi;
            Geometry.OffsetToIndex(i, xi, yi, zi);
            if (Geometry.GetSlownesses()[xi][yi][zi] > 0.0)
              {
                DensDiff(i) = (PredDens(i) - Model(i + ncells)) / pow<2>(refdens);
                CondDiff(i) = (PredCond(i) - Model(i + 2 * ncells)) / pow<2>(refcond);
                SaltDiff(i) = SaltTerm(Model, i);
              }
            else
              {
                DensDiff(i) = 0.0;
                CondDiff(i) = 0.0;
                SaltDiff(i) = 0.0;
              }
          }

        //the gradient for the slownesses
        ublas::subrange(Gradient, 0, ncells) = 2.0
            * (DensityTransform->Derivative(ublas::subrange(Model, 0, ncells), DensDiff)
                + ConductivityTransform->Derivative(ublas::subrange(Model, 0, ncells),
                    CondDiff));
        ublas::subrange(Gradient, 0, ncells) = ublas::element_prod(
            ublas::subrange(Gradient, 0, ncells), SaltDiff);

        //gradient for density
        ublas::subrange(Gradient, ncells, 2 * ncells) = -2.0
            * ublas::element_prod(DensDiff, SaltDiff);
        //gradient for conductivity
        ublas::subrange(Gradient, 2 * ncells, 3 * ncells) = -2.0
            * ublas::element_prod(CondDiff, SaltDiff);

        for (size_t i = 0; i < ncells; ++i)
          {
            const double misterm = 2.0
                * (pow<2>(DensDiff(i) * refdens) + pow<2>(CondDiff(i) * refcond));
            Gradient(i) += misterm * (Model(i) - 1.0 / saltvel) / pow<2>(refslow);
            Gradient(i + ncells) += misterm * (Model(i + ncells) - saltdens)
                / pow<2>(refdens);
            Gradient(i + 2 * ncells) += misterm * (Model(i + 2 * ncells) - 1.0 / saltres)
                / pow<2>(refcond);
          }

        return Gradient;
      }
  }
