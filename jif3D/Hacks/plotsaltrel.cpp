//============================================================================
// Name        : plotsaltrel.cpp
// Author      : 5 Mar 2012
// Version     : 
// Copyright   : 2012, mm489
//============================================================================

#include "../Inversion/ModelTransforms.h"
#include "../Gravity/ThreeDGravityModel.h"
#include "../Joint/SaltRelConstraint.h"
#include "../ModelBase/VTKTools.h"
#include <fstream>

//class DensTransform: public jif3D::GeneralModelTransform
//  {
//public:
//  DensTransform()
//    {
//    }
//  virtual ~DensTransform()
//    {
//    }
//  virtual jif3D::rvec GeneralizedToPhysical(const jif3D::rvec &FullModel) const
//    {
//      jif3D::rvec result(1);
//      result(0) = 2.2;
//      return result;
//    }
//  virtual jif3D::rvec PhysicalToGeneralized(const jif3D::rvec &FullModel) const
//    {
//      jif3D::rvec result(1);
//      result(0) = 2.2;
//      return result;
//    }
//  virtual jif3D::rvec Derivative(const jif3D::rvec &FullModel,
//      const jif3D::rvec &Derivative) const
//    {
//      jif3D::rvec result(1);
//      result(0) = 0.0;
//      return result;
//    }
//  };

int main()
  {
    jif3D::ThreeDGravityModel RModel;
    RModel.SetDensities().resize(boost::extents[1][1][1]);
    boost::shared_ptr<jif3D::GeneralModelTransform> DensTrans(
        new jif3D::DensityTransform(
            boost::shared_ptr<jif3D::GeneralModelTransform>(
                new jif3D::ModelCopyTransform()),RModel));

    boost::shared_ptr<jif3D::GeneralModelTransform> CondTrans(
        new jif3D::ConductivityTransform(
            boost::shared_ptr<jif3D::GeneralModelTransform>(
                new jif3D::ModelCopyTransform()),RModel));
    jif3D::SaltRelConstraint SaltRel(DensTrans, CondTrans);

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

                jif3D::rvec Model(3);
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
    jif3D::Write3DModelToVTK("saltrel.vtk", "SaltReal", VelAxis, LogCondAxis, DensAxis,
        Misfit);

    jif3D::ConductivityTransform SlowCond(
        boost::shared_ptr<jif3D::GeneralModelTransform>(new jif3D::ModelCopyTransform()),RModel);
    jif3D::rvec Model(nvel);
    size_t i = 0;
    for (double vel = minvel; vel < maxvel; vel += velstep)
      {
        Model(i) = 1.0 / vel;
        ++i;
      }
    jif3D::rvec Cond = SlowCond.GeneralizedToPhysical(Model);
    std::ofstream outfile("velres.dat");
    for (i = 0; i < nvel; ++i)
      {
        outfile << 1.0 / Model(i) << " " << 1.0 / Cond(i) << "\n";
      }
  }
