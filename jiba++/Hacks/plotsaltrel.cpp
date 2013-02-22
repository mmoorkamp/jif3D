//============================================================================
// Name        : plotsaltrel.cpp
// Author      : 5 Mar 2012
// Version     : 
// Copyright   : 2012, mm489
//============================================================================

#include "../Inversion/ModelTransforms.h"
#include "../Joint/SaltRelConstraint.h"
#include "../Joint/SaltRelTrans.h"
#include "../ModelBase/VTKTools.h"
#include <fstream>

//class DensTransform: public jiba::GeneralModelTransform
//  {
//public:
//  DensTransform()
//    {
//    }
//  virtual ~DensTransform()
//    {
//    }
//  virtual jiba::rvec GeneralizedToPhysical(const jiba::rvec &FullModel) const
//    {
//      jiba::rvec result(1);
//      result(0) = 2.2;
//      return result;
//    }
//  virtual jiba::rvec PhysicalToGeneralized(const jiba::rvec &FullModel) const
//    {
//      jiba::rvec result(1);
//      result(0) = 2.2;
//      return result;
//    }
//  virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
//      const jiba::rvec &Derivative) const
//    {
//      jiba::rvec result(1);
//      result(0) = 0.0;
//      return result;
//    }
//  };

int main()
  {

    boost::shared_ptr<jiba::GeneralModelTransform> DensTrans(
        new jiba::DensityTransform(
            boost::shared_ptr<jiba::GeneralModelTransform>(
                new jiba::ModelCopyTransform())));

    boost::shared_ptr<jiba::GeneralModelTransform> CondTrans(
        new jiba::SlowCondTrans(
            boost::shared_ptr<jiba::GeneralModelTransform>(
                new jiba::ModelCopyTransform())));
    jiba::SaltRelConstraint SaltRel(DensTrans, CondTrans);

    double minvel = 2000;
    double maxvel = 6000;
    double velstep = 100;
    double logmincond = -4;
    double logmaxcond = 0.5;
    double logcondstep = 0.05;
    double mindens = 1.5;
    double maxdens = 4.0;
    double densstep = 0.05;

    const size_t nvel = (maxvel - minvel) / velstep + 1;
    const size_t nlogcond = (logmaxcond - logmincond) / logcondstep + 1;
    const size_t ndens = (maxdens - mindens) / densstep + 1;

    boost::multi_array<double, 1> VelAxis(boost::extents[nvel]), LogCondAxis(
        boost::extents[nlogcond]), DensAxis(boost::extents[ndens]);
    boost::multi_array<double, 3> Misfit(boost::extents[nvel][nlogcond][ndens]);

    std::fill_n(VelAxis.origin(), nvel, velstep);
    std::fill_n(LogCondAxis.origin(), nlogcond, logcondstep);
    std::fill_n(DensAxis.origin(), ndens, densstep);
    size_t currvelindex = 0, currcondindex = 0, currdensindex = 0;
    for (double vel = minvel; vel < maxvel; vel += velstep)
      {
        currcondindex = 0;

        for (double logcond = logmincond; logcond < logmaxcond; logcond += logcondstep)
          {

            currdensindex = 0;

            for (double dens = mindens; dens < maxdens; dens += densstep)
              {

                jiba::rvec Model(3);
                Model(0) = 1.0 / vel;
                Model(1) = dens;
                Model(2) = pow(10, logcond);
                Misfit[currvelindex][currcondindex][currdensindex] = SaltRel.CalcMisfit(
                    Model);
                ++currdensindex;
              }
            ++currcondindex;
          }
        ++currvelindex;
      }
    jiba::Write3DModelToVTK("saltrel.vtk", "SaltReal", VelAxis, LogCondAxis, DensAxis,
        Misfit);

    jiba::SlowCondTrans SlowCond(
        boost::shared_ptr<jiba::GeneralModelTransform>(new jiba::ModelCopyTransform()));
    jiba::rvec Model(nvel);
    size_t i = 0;
    for (double vel = minvel; vel < maxvel; vel += velstep)
      {
        Model(i) = 1.0 / vel;
        ++i;
      }
    jiba::rvec Cond = SlowCond.GeneralizedToPhysical(Model);
    std::ofstream outfile("velres.dat");
    for (i = 0; i < nvel; ++i)
      {
        outfile << 1.0 / Model(i) << " " << 1.0 / Cond(i) << "\n";
      }
  }