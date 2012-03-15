//============================================================================
// Name        : plotsaltrel.cpp
// Author      : 5 Mar 2012
// Version     : 
// Copyright   : 2012, mm489
//============================================================================

#include "../Inversion/ModelTransforms.h"
#include "../Joint/SaltRelConstraint.h"
#include <fstream>

class DensTransform: public jiba::GeneralModelTransform
  {
public:
  DensTransform()
    {
    }
  virtual ~DensTransform()
    {
    }
  virtual jiba::rvec GeneralizedToPhysical(const jiba::rvec &FullModel) const
    {
      jiba::rvec result(1);
      result(0) = 2.2;
      return result;
    }
  virtual jiba::rvec PhysicalToGeneralized(const jiba::rvec &FullModel) const
    {
      jiba::rvec result(1);
      result(0) = 2.2;
      return result;
    }
  virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
      const jiba::rvec &Derivative) const
    {
      jiba::rvec result(1);
      result(0) = 0.0;
      return result;
    }
  };

int main()
  {

    boost::shared_ptr<jiba::GeneralModelTransform> DensTrans(new DensTransform());
    boost::shared_ptr<jiba::GeneralModelTransform> CondTrans(
        new jiba::ConductivityTransform(
            boost::shared_ptr<jiba::GeneralModelTransform>(
                new jiba::ModelCopyTransform())));
    jiba::SaltRelConstraint SaltRel(DensTrans, CondTrans);

    double minvel = 2000;
    double maxvel = 6000;
    double velstep = 100;
    double logmincond = -4;
    double logmaxcond = 0.5;
    double logcondstep = 0.05;
    std::ofstream outfile("saltrel.dat");
    for (double vel = minvel; vel < maxvel; vel += velstep)
      {
        for (double logcond = logmincond; logcond < logmaxcond; logcond += logcondstep)
          {
            jiba::rvec Model(3);
            Model(0) = 1.0 / vel;
            Model(1) = 2.2;
            Model(2) = pow(10, logcond);
            double misfit = SaltRel.CalcMisfit(Model);
            outfile << vel << " " << logcond << " " << misfit << std::endl;
          }
      }

  }
