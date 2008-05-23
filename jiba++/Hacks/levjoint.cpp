#include <iostream>
#include "Iso1DMTObjective.h"
#include "MTStation.h"
#include "gentypes.h"
#include "lm.h"
#include "MTFitSetup.h"
#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/assign/std/vector.hpp> 
#include <numeric>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

using namespace std;
using namespace boost::assign;
C1DMTSynthData MTInputData;
MTStation Best;
Iso1DMTObjective MTObjective(MTInputData);
double inputgrav = 0.0;

static const double Grav_const = 6.67428e-8; // in units cm^3/g s

struct CalcInfSheet : public std::binary_function<double, double, double>
  {
  double operator()(const double thick, const double density) const
    {
      return 2.0 * M_PI * Grav_const * thick * density * 1000.0;
    }
  };

struct ResDensRel : public std::unary_function<double, double>
  {
  double operator()(const double res) const
    {
      return 2.26 + 0.147 * log10(res) - 0.0156 * pow(log10(res), 2); // some funny relationship
    }
  };

template <class VectorType> double CalcGrav(const VectorType &density,
    const VectorType &thickness)
  {
    return inner_product(thickness.begin(), thickness.end(), density.begin(),
        0.0, plus<double>(), CalcInfSheet());

  }

void AddNoise(std::complex<double> &impelement, double &noiseest,
    const double noiselevel, boost::lagged_fibonacci607 &generator)
  {
    complex<double> oldimp(impelement);
    impelement = boost::variate_generator<boost::lagged_fibonacci607&,
    boost::normal_distribution<> >(generator,
        boost::normal_distribution<>(impelement.real(),
            fabs(impelement.real()* noiselevel)))();
    impelement += complex<double>(0.0,1.0)*boost::variate_generator<boost::lagged_fibonacci607&,
    boost::normal_distribution<> >(generator,
        boost::normal_distribution<>(oldimp.imag(),
            fabs(oldimp.imag()* noiselevel)))();
    noiseest = noiselevel * abs(oldimp);
  }

void misfit(double *p, double *x, int m, int n, void *data)
  {
    const unsigned int nlayers = m/2;
    const unsigned int nparam = m;
    ttranscribed mtmember(nparam), densities(nlayers), thicknesses(nlayers);
    for (unsigned int i = 0; i < nparam; ++i)
      {
        mtmember(i) = p[i];
      }
    copy(p+nlayers, p+m, thicknesses.begin());
    transform(p, p+nlayers, densities.begin(), ResDensRel());
    MTObjective.CalcPerformance(mtmember);
    for (size_t i = 0; i < MTObjective.GetMisfit().size(); ++i)
      {
        x[i] = MTObjective.GetMisfit()(i);
      }
    double currgrav = CalcGrav(densities, thicknesses);
    x[MTObjective.GetMisfit().size() +1] = pow(inputgrav - currgrav, 2);
    //x[MTObjective.GetMisfit().size() +1] = 0;
  }

int main()
  {
    //create input model
    trealdata resistivities, thicknesses, densities, frequencies;
    resistivities += 1.0, 20.0, 2.0, 100.0;
    thicknesses += 1.2, 2.2, 2.5, 2.0;
    frequencies += 1.8600e+000, 8.7100e-001, 4.0700e-001, 1.9100e-001,
        8.9100e-002, 4.1700e-002, 1.9500e-002, 9.1200e-003, 4.2700e-003,
        2.0000e-003, 9.3300e-004, 4.3700e-004, 2.0400e-004, 9.5500e-005,
        4.4700e-005, 2.0900e-005, 9.7700e-006, 4.5700e-006, 2.1400e-006,
        1.0000e-006;
    transform(resistivities.begin(), resistivities.end(),
        back_inserter(densities), ResDensRel());

    //calculate MT input data
    MTInputData.SetResistivities(resistivities);
    MTInputData.SetThicknesses(thicknesses);
    MTInputData.SetFrequencies(frequencies);
    MTInputData.GetData();
    
    //add noise to the data
    const double noiselevel = 0.05;
    boost::lagged_fibonacci607
        generator(static_cast<unsigned int>(std::time(0)));
    double zxxerr, zxyerr, zyxerr, zyyerr;
    for (unsigned int i = 0; i < MTInputData.GetMTData().size(); ++i)
      {
        AddNoise(MTInputData.SetMTData().at(i).SetZxx(), zxxerr, noiselevel, generator);
        AddNoise(MTInputData.SetMTData().at(i).SetZxy(), zxyerr, noiselevel, generator);
        AddNoise(MTInputData.SetMTData().at(i).SetZyx(), zyxerr, noiselevel, generator);
        AddNoise(MTInputData.SetMTData().at(i).SetZyy(), zyyerr, noiselevel, generator);
        MTInputData.SetMTData().at(i).SetErrors(zxxerr, zxyerr, zyxerr, zyyerr);
      }
    
    // create objective function
    MTObjective = Iso1DMTObjective(MTInputData);
    //calculate Gravity input data
    cout << "Input Model densities: ";
    copy(densities.begin(), densities.end(),
        ostream_iterator<double>(cout, " "));
    cout << endl;
    inputgrav = CalcGrav(thicknesses, densities);
    cout << "True Gravity: " << inputgrav << endl;
    //steup fit for MT
    MTObjective.AppendFitParameters(&MTTensor::GetRhoxy, &MTTensor::GetdRhoxy,
        noiselevel);
    MTObjective.AppendFitParameters(&MTTensor::GetPhixy, &MTTensor::GetdPhixy,
        noiselevel);

    const int nlayers = 4;
    const int nparams = nlayers * 2;
    ttranscribed mtmember(nparams);

    const double startres = 2.0;
    const double minres = -1;
    const double maxres = 5;
    const double startthick = 1.0;
    const double minthick = 0.1;
    const double maxthick = 100.0;
    fill_n(mtmember.begin(), nlayers, startres);
    fill_n(mtmember.begin()+nlayers, nlayers, startthick);

    MTObjective.CalcPerformance(mtmember);

    //double p[nparams],lb[nparams],ub[nparams]; 
    double *p = new double[nparams];
    double *lb = new double[nparams];
    double *ub = new double[nparams];
    double *x;

    int n = MTObjective.GetSynthData().size()+1;
    //we have n data for the MT plus one for gravity
    x = new double[n];
    double opts[LM_OPTS_SZ], info[LM_INFO_SZ];

    opts[0]=LM_INIT_MU;
    opts[1]=1E-15;
    opts[2]=1E-15;
    opts[3]=1E-20;
    opts[4]=LM_DIFF_DELTA; // relevant only if the finite difference jacobian version is used 

    for (int i = 0; i < nlayers; ++i)
      {
        p[i] = startres;
        lb[i] = minres;
        ub[i] = maxres;
      }
    for (int i = 0; i < nlayers; ++i)
      {
        p[i+nlayers] = startthick;
        lb[i+nlayers] = minthick;
        ub[i+nlayers] = maxthick;
      }

    for (int i=0; i<n; i++)
      x[i]=0.0;

    double ret=dlevmar_bc_dif(misfit, p, x, nparams, n, lb, ub, 200, opts,
        info, NULL, NULL, NULL); // no jacobian
    cout << "Levenberg-Marquardt returned " << ret << " in " << info[5]
        << "iter, reason " << info[6] << endl;

    cout << endl << " Minimization info:"<< endl;
    for (int i=0; i<LM_INFO_SZ; ++i)
      cout << info[i] << " ";
    cout << endl;

    ttranscribed mtbest(nlayers * 2);
    for (int i = 0; i < nlayers * 2; ++i)
      {
        mtbest(i) = p[i];
      }
    cout << endl;
    double mtdiff = MTObjective.CalcPerformance(mtbest);

    cout << "MT-Fit: " << mtdiff << endl;
    ostringstream filename;
    filename << "best_lev";

    cout << " Saved as : " << filename.str() << endl;
    MTObjective.WriteData(filename.str());
    MTObjective.WriteModel(filename.str()+"_mt.mod");
    MTObjective.WritePlot(filename.str()+"_mt.plot");

    delete []x;
    delete []p;
    delete []lb;
    delete []ub;
    return 0;
  }
