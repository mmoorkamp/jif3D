#include "../Global/FileUtil.h"
#include <iostream>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/sum.hpp>
#include <boost/accumulators/statistics/weighted_sum.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/p_square_quantile.hpp>
#include "../MT/ReadWriteImpedances.h"

using namespace boost::accumulators;

int main()
  {
    std::string filename = jif3D::AskFilename("MT data:");
    std::vector<double> Freq, StatX, StatY, StatZ;
    jif3D::rvec Imp, Err;
    jif3D::ReadImpedancesFromNetCDF(filename, Freq, StatX, StatY, StatZ, Imp, Err);
    const size_t nstat = StatX.size();
    const size_t nfreq = Freq.size();
    typedef accumulator_set<double, stats<tag::p_square_quantile> > accumulator_t;
    accumulator_set<double, stats<tag::median(with_p_square_quantile)> > acc;
    for (size_t i = 0; i < nfreq; ++i)
      {
        accumulator_t acc0(quantile_probability = 0.001);
        accumulator_t acc1(quantile_probability = 0.999);
        accumulator_t acc2(quantile_probability = 0.5);
        for (size_t k = 2; k < 6; ++k)
          {
            for (size_t j = i * nstat * 8; j < (i + 1) * nstat * 8; j += 8)
              {
                acc0(Imp(j + k));
                acc1(Imp(j + k));
                acc2(Imp(j + k));
              }

            double minthresh = p_square_quantile(acc0);
            double maxthresh = p_square_quantile(acc1);
            double median = p_square_quantile(acc2);
            std::cout << "Frequency: " << Freq.at(i) << " " << minthresh << " " << median
                << " " << maxthresh << std::endl;
          }
        std::cout << std::endl;
      }
  }
