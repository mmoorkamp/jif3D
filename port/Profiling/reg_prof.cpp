//============================================================================
// Name        : reg_prof.cpp
// Author      : 16 Aug 2012
// Version     : 
// Copyright   : 2012, mm489
//============================================================================

#include "../Gravity/test_common.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include "../Regularization/GradientRegularization.h"
#include "../Global/VecMat.h"
#include <boost/date_time/posix_time/posix_time.hpp>
#include <cstdlib>

int main()
  {
    for (size_t nelem = 10; nelem < 100; nelem += 1)
      {

        jif3D::ThreeDGravityModel GravModel;
        const size_t nx = nelem;
        const size_t ny = nelem;
        const size_t nz = nelem;
        GravModel.SetMeshSize(nx, ny, nz);

        const size_t msize = GravModel.GetDensities().num_elements();
        jif3D::rvec StartModel(msize), PertModel(msize);
        jif3D::rvec ConstMod(msize);
        std::fill(ConstMod.begin(), ConstMod.end(), 1.0);
        std::generate(StartModel.begin(), StartModel.end(), rand);
        std::generate(PertModel.begin(), PertModel.end(), rand);

        boost::posix_time::ptime firststarttime =
            boost::posix_time::microsec_clock::local_time();
        jif3D::GradientRegularization Regularization(GravModel, 0.0);
        boost::posix_time::ptime firstendtime =
            boost::posix_time::microsec_clock::local_time();
        Regularization.SetReferenceModel(StartModel);
        Regularization.SetDataError(StartModel);
        Regularization.SetXWeight(5.0);
        Regularization.SetYWeight(4.0);
        Regularization.SetZWeight(3.0);
        boost::posix_time::ptime secondstarttime =
            boost::posix_time::microsec_clock::local_time();
        double zero = Regularization.CalcMisfit(StartModel + ConstMod);
        boost::posix_time::ptime secondendtime =
            boost::posix_time::microsec_clock::local_time();
        std::cout << nelem << " " << (firstendtime - firststarttime).total_microseconds()
            << " " << (secondendtime - secondstarttime).total_microseconds() << std::endl;
      }
  }
