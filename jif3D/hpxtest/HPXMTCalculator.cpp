//============================================================================
// Name        : X3DMTCalculator.cpp
// Author      : Jul 7, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================

#include <hpx/config.hpp>
#include <hpx/include/lcos.hpp>
#include <cassert>
#include <complex>
#include <fstream>
#include <iomanip>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/multi_array.hpp>

#include "../Global/FatalException.h"
#include "../Global/convert.h"
#include "../MT/MTUtils.h"
#include "HPXMTCalculator.h"
#include "CalcFreq.h"



namespace fs = boost::filesystem;



namespace jif3D
  {






    void HPXMTCalculator::CleanUp()
      {
        //std::string NameRoot(ObjectID());
        /*		fs::directory_iterator end_itr; // default construction yields past-the-end
         //go through the directory and delete any file that starts with NameRoot
         for (fs::directory_iterator itr(fs::current_path()); itr != end_itr;
         ++itr)
         {
         if (boost::algorithm::starts_with(itr->path().filename().string(),
         NameRoot))
         {
         fs::remove_all(itr->path().filename());
         }
         }*/
      }

    HPXMTCalculator::HPXMTCalculator(boost::filesystem::path TDir, std::string x3d):
        X3DName(x3d)
      {
        if (!fs::is_directory(TDir))
          throw FatalException("TDir is not a directory: " + TDir.string());
        TempDirName = TDir.string();
        //we make sure that the .hnk files are there
        //this is a common problem and when we check for use later
        //we are inside an openmp thread and swallow all sensible error messages.
        if (!CheckHNK(fs::path()))
          {
            throw jif3D::FatalException("Cannot find .hnk files in current directory! ");
          }
      }

    HPXMTCalculator::~HPXMTCalculator()
      {
        CleanUp();
      }




    rvec HPXMTCalculator::Calculate(const X3DModel &Model, size_t minfreqindex,
        size_t maxfreqindex)
      {
        using hpx::lcos::future;
        using hpx::async;
        using hpx::wait_all;

        assert(minfreqindex <= maxfreqindex);
        const size_t maxindex = std::min(maxfreqindex, Model.GetFrequencies().size());
        const size_t nfreq = maxindex - minfreqindex;
        const size_t nmeas = Model.GetMeasPosX().size();
        const size_t nmodx = Model.GetXCoordinates().size();
        const size_t nmody = Model.GetYCoordinates().size();

        jif3D::rvec result(nmeas * nfreq * 8);
        result.clear();
        //we make a call to the coordinate functions to make sure
        //that we have updated the coordinate information and cached it
        //only then the subsequent calls are thread safe
        Model.GetXCoordinates();
        Model.GetYCoordinates();
        Model.GetZCoordinates();
        std::vector<double> BGDepths(Model.GetBackgroundThicknesses().size(), 0.0);
        std::partial_sum(Model.GetBackgroundThicknesses().begin(),
            Model.GetBackgroundThicknesses().end(), BGDepths.begin());
        CompareDepths(BGDepths, Model.GetZCoordinates());

        std::vector<hpx::naming::id_type> localities =
                hpx::find_all_localities();
        std::cout << "Number of localities: " << localities.size() << std::endl;



        // Get the number of worker OS-threads in use by this locality.
        std::size_t const os_threads = hpx::get_os_thread_count();

        // Populate a set with the OS-thread numbers of all OS-threads on this
        // locality. When the hello world message has been printed on a particular
        // OS-thread, we will remove it from the set.
        std::set<std::size_t> attendance;
        for (std::size_t os_thread = 0; os_thread < os_threads; ++os_thread)
            attendance.insert(os_thread);

        std::cout << "Number of threads: " << attendance.size() << std::endl;

        std::vector<future<jif3D::rvec>> FreqResult;
        FreqResult.reserve(nfreq);
        CalculateFrequency_action FreqCalc;
    	boost::posix_time::ptime calcstarttime =
    			boost::posix_time::microsec_clock::local_time();
        for (int i = minfreqindex; i < maxindex; ++i)
          {
        	hpx::naming::id_type const locality_id = localities.at(i % localities.size());
            //rvec freqresult = CalculateFrequency(Model, i, TempDir);
            FreqResult.push_back(async(FreqCalc, locality_id, Model, i, TempDirName, X3DName));

          }
        wait_all(FreqResult);
    	boost::posix_time::ptime calcendtime =
    			boost::posix_time::microsec_clock::local_time();
        std::cout << "Nfreq: " << FreqResult.size() << std::endl;
        for (int i = minfreqindex; i < maxindex; ++i)
          {
            size_t currindex = i - minfreqindex;
            size_t startindex = nmeas * currindex * 8;
            jif3D::rvec imp = FreqResult[currindex].get();
            std::copy(imp.begin(), imp.end(), result.begin() + startindex);
          }
    	boost::posix_time::ptime gatherendtime =
    			boost::posix_time::microsec_clock::local_time();

    	std::cout << " Times: " << std::endl;
    	std::cout << "Calc: " << (calcendtime - calcstarttime).total_microseconds()
    			<< std::endl;
    	std::cout << "gather: " << (gatherendtime - calcendtime).total_microseconds()
    			<< std::endl << std::endl;
        return result;

      }





  }

