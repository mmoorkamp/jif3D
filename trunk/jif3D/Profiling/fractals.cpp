////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2012 Andrew Kemp
//  with additions from Max Moorkamp
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <hpx/hpx_init.hpp>
#include <hpx/include/actions.hpp>
#include <hpx/include/util.hpp>
#include <hpx/include/lcos.hpp>
#include <hpx/hpx_fwd.hpp>
#include <hpx/parallel/numeric.hpp>
#include <hpx/parallel/algorithms/transform.hpp>
#include <hpx/parallel/execution_policy.hpp>
#include <vector>
#include <algorithm>
#include <parallel/algorithm>

#include <boost/cstdint.hpp>
#include <boost/format.hpp>
#include <boost/serialization/serialization.hpp>

int const sizeY = 512;
int const sizeX = sizeY;

struct my_complex
  {
  double real;
  double imag;
  my_complex(double r, double i) :
      real(r), imag(i)
    {
    }
  my_complex() :
      real(), imag()
    {
    }
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
    {
      ar & real;
      ar & imag;
    }
  };

int fractals(double x0, double y0, int max_iteration);
int fractals_complex(my_complex P, int max_iteration);
std::vector<int> fractals_chunk(int i, int max_iteration);

// This is to generate the required boilerplate we need for the remote
// invocation to work.
HPX_PLAIN_ACTION(fractals_chunk, fractals_action)

///////////////////////////////////////////////////////////////////////////////
int fractals(double x0, double y0, int max_iteration)
  {
    double x = 0, y = 0;
    int iteration = 0;
    while (x * x + y * y < 2 * 2 && iteration < max_iteration)
      {
        double xtemp = x * x - y * y + x0;
        y = 2 * x * y + y0;
        x = xtemp;

        ++iteration;
      }
    return iteration;
  }

int fractals_complex(my_complex P, int max_iteration)
  {
    return fractals(P.real, P.imag, max_iteration);
  }

std::vector<int> fractals_chunk(int i, int max_iteration)
  {
    std::vector<int> Result;
    Result.reserve(sizeY);
    for (int j = 0; j < sizeY; j++)
      {
        double x0 = (double) i * 3.5f / (double) sizeX - 2.5f;
        double y0 = (double) j * 2.0f / (double) sizeY - 1.0f;
        my_complex P(x0, y0);
        Result.push_back(fractals_complex(P, max_iteration));
      }
    return Result;
  }

///////////////////////////////////////////////////////////////////////////////
int hpx_main()
  {

    hpx::util::high_resolution_timer t;

      {
        using namespace std;

        using hpx::future;
        using hpx::async;
        using hpx::wait_all;

        int const max_iteration = 255;

        vector<future<vector<int> > > iteration;
        iteration.reserve(sizeX);

        hpx::id_type const here = hpx::find_here();
        fractals_action fractal_line;
        std::vector<hpx::naming::id_type> localities = hpx::find_all_localities();

        std::cout << "Running on " << localities.size() << " locations \n";
        std::cout << "Using " << hpx::get_num_worker_threads() << " threads \n";

        std::cout << "Initial setup completed in " << t.elapsed()
            << "s. Initializing and running futures...\n";
        t.restart();

        for (int i = 0; i < sizeX; i++)
          {
            hpx::naming::id_type const locality_id = localities.at(i % localities.size());
            iteration.push_back(async(fractal_line, locality_id, i, max_iteration));
          }
        wait_all(iteration);

        double hpxcalctime = t.elapsed();
        std::cout << sizeX * sizeY << " calculations run in " << hpxcalctime
            << "s.\n Transferring from futures to general memory...\n";
        t.restart();

        vector<int> HPXActionResult;
        for (int i = 0; i < sizeX; ++i)
          {
            vector<int> it = iteration[i].get();
            std::copy(it.begin(), it.end(), back_inserter(HPXActionResult));
          }
        double hpxtransfertime = t.elapsed();
        std::cout << "Transfer process completed in " << hpxtransfertime << "s. \n";

        std::cout << "Setting up new input \n" << std::endl;
        t.restart();
        std::vector<my_complex> InValues;
        std::vector<int> OutValues(sizeX * sizeY), OpenMPValues(sizeX * sizeY);
        InValues.reserve(sizeX * sizeY);
        for (int i = 0; i < sizeX; i++)
          {
            for (int j = 0; j < sizeY; j++)
              {
                double x0 = (double) i * 3.5f / (double) sizeX - 2.5f;
                double y0 = (double) j * 2.0f / (double) sizeY - 1.0f;
                InValues.push_back(my_complex(x0, y0));
              }
          }
        std::cout << InValues.size() << " Values initialized in " << t.elapsed()
            << "s. \n";

        t.restart();
        hpx::parallel::transform(hpx::parallel::parallel_execution_policy(),
            InValues.begin(), InValues.end(), OutValues.begin(),
            [max_iteration](my_complex P) -> int
              { return fractals_complex(P,max_iteration);});
        double hpxalgotime = t.elapsed();
        std::cout << "HPX parallel calculation finished in " << hpxalgotime << "s. \n";
        bool IsEqualHPX = std::equal(HPXActionResult.begin(), HPXActionResult.end(),
            OutValues.begin());

        t.restart();
        __gnu_parallel::transform(InValues.begin(), InValues.end(), OutValues.begin(),
            [max_iteration](my_complex P) -> int
              { return fractals_complex(P,max_iteration);});
        double gnualgotime = t.elapsed();
        std::cout << "OpenMP parallel calculation finished in " << gnualgotime << "s. \n";
        bool IsEqualOpenMP = std::equal(HPXActionResult.begin(), HPXActionResult.end(),
            OutValues.begin());

        t.restart();
        std::transform(InValues.begin(), InValues.end(), OutValues.begin(),
            [max_iteration](my_complex P) -> int
              { return fractals_complex(P,max_iteration);});
        double serialalgotime = t.elapsed();
        bool IsEqualSerial = std::equal(HPXActionResult.begin(), HPXActionResult.end(),
            OutValues.begin());
        std::cout << "Serial calculation finished in " << serialalgotime << "s. \n";

        std::cout << "HPX action and hpx::parallel results are equal " << IsEqualHPX
            << " \n";
        std::cout << "HPX action and gcc::parallel results are equal " << IsEqualOpenMP
            << " \n";
        std::cout << "HPX action and serial results are equal " << IsEqualSerial << " \n";

        std::cout << "\n \n Times: \n";
        std::cout << "Hpx actions: " << hpxcalctime + hpxtransfertime << "\n";
        std::cout << "Hpx algorithm: " << hpxalgotime << "\n";
        std::cout << "GNP algorithm: " << gnualgotime << "\n";
        std::cout << "Serial algorithm: " << serialalgotime << "\n";
      }

    return hpx::finalize(); // Handles HPX shutdown
  }

int main(int argc, char* argv[])
  {
    return hpx::init(argc, argv);
  }

