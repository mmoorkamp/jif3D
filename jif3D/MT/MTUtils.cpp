//============================================================================
// Name        : MTUtils.cpp
// Author      : 10 Mar 2014
// Version     : 
// Copyright   : 2014, mm489
//============================================================================

#include "MTUtils.h"
#include "../Global/NumUtil.h"
#include "../Global/convert.h"
#include "../Global/FatalException.h"
#include "ReadWriteX3D.h"

#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/math/constants/constants.hpp>

namespace jif3D
  {
    const std::string runext = "_run";

//associate a type of calculation with a name string
    const std::map<X3DModel::ProblemType, std::string> Extension =
      {
        { X3DModel::MT, "MT" },
        { X3DModel::EDIP, "EDIP" },
        { X3DModel::MDIP, "MDIP" } };

    namespace fs = boost::filesystem;

    std::string MakeUniqueName(const std::string &NameRoot, X3DModel::ProblemType Type,
        const size_t FreqIndex)
      {
        //we assemble the name from the id of the process
        //and the address of the current object
        std::string result(NameRoot);
        //the type of calculation
        result += Extension.find(Type)->second;
        //and the frequency index
        result += jif3D::stringify(FreqIndex);
        return result;
      }
//create a script that changes to the correct directory
//and executes x3d in that directory
    void MakeRunFile(const std::string &NameRoot, const std::string &DirName,
        const std::string &X3DName)
      {
        std::string RunFileName = NameRoot + runext;
        fs::create_directory(DirName);
        //boost::iostreams::stream<boost::iostreams::file_descriptor_sink> runfile(
        //    RunFileName.c_str());
        std::ofstream runfile(RunFileName.c_str());
        runfile << "#!/bin/bash\n";
        runfile << " echo " "Executing: $BASH_SOURCE" "\n ";
        runfile << "cd " << DirName << "\n";
        runfile << X3DName << " \n";
        runfile << "cd ..\n";
        /*        if (runfile.rdbuf())
         {
         runfile.rdbuf()->pubsync();
         }
         ::fdatasync(runfile->handle());*/
        //we also copy the necessary *.hnk files
        //from the current directory to the work directory
        CopyHNK(fs::current_path(), DirName);
      }

//copy the .hnk files for x3d from SourceDir to TargetDir
    void CopyHNK(const fs::path &SourceDir, const fs::path &TargetDir)
      {
        //copy file fails with an exception if the target exists
        //so we check before we do the actual copy
        if (!CheckHNK(TargetDir))
          {
            fs::copy_file(SourceDir / "ndec15.hnk", TargetDir / "ndec15.hnk");
            fs::copy_file(SourceDir / "ndec20.hnk", TargetDir / "ndec20.hnk");
            fs::copy_file(SourceDir / "ndec30.hnk", TargetDir / "ndec30.hnk");
            fs::copy_file(SourceDir / "ndec40.hnk", TargetDir / "ndec40.hnk");
          }
      }

    void CheckField(const std::vector<std::complex<double> > &Field, size_t nelem)
      {
        if (Field.size() != nelem)
          throw jif3D::FatalException(
              "Number of read in elements in Field: " + jif3D::stringify(Field.size())
                  + " does not match expected: " + jif3D::stringify(nelem), __FILE__,
              __LINE__);
      }

    void CompareDepths(const std::vector<double> &BGDepths,
        const jif3D::ThreeDModelBase::t3DModelDim &ModelDepths)
      {
        size_t mindex = 0;
        for (size_t i = 0; i < BGDepths.size(); ++i)
          {
            double rd = round(BGDepths[i]);
            while (mindex < ModelDepths.size() && round(ModelDepths[mindex]) < rd)
              {
                ++mindex;
              }
            if (mindex < ModelDepths.size() && !jif3D::roughlyEqual<float,float>()(round(ModelDepths[mindex]),rd))
              {
                throw jif3D::FatalException(
                    "Depth to background layer: " + jif3D::stringify(BGDepths[i])
                        + " does not match grid cell depth: "
                        + jif3D::stringify(ModelDepths[mindex]), __FILE__, __LINE__);
              }
          }
      }

//check that the .hnk file for x3d are in a certain directory
    bool CheckHNK(const fs::path &TargetDir)
      {
        return fs::exists(TargetDir / "ndec15.hnk")
            && fs::exists(TargetDir / "ndec20.hnk")
            && fs::exists(TargetDir / "ndec30.hnk")
            && fs::exists(TargetDir / "ndec40.hnk");
      }

//execute the script that runs x3d
    void RunX3D(const std::string &NameRoot)
      {
        const std::string runname = NameRoot + runext;
        int tries = 0;
        while (!fs::exists(runname) && tries < 10)
          {
            ++tries;
          }
        if (tries >= 10)
          {
            throw FatalException("Cannot find run script: " + runname, __FILE__,
            __LINE__);
          }

        //instead of making the script executable
        //we run a bash with the scriptname as an argument
        //this turns out to be more robust
        const std::string execname = "bash " + runname + " > " + runname + "_log";
        /*        namespace bp = ::boost::process;
         const std::string execname = bp::find_executable_in_path("bash");
         std::vector<std::string> args;
         args.push_back("bash");
         args.push_back(runname);

         bp::context ctx;
         ctx.stdout_behavior = bp::capture_stream();
         ctx.environment = bp::self::get_environment();

         bp::child c = bp::launch(execname, args, ctx);
         bp::pistream &is = c.get_stdout();
         std::string line;
         while (std::getline(is, line))
         BOOST_LOG_TRIVIAL(debug)<< line << std::endl;

         bp::posix_status s = c.wait();
         if (s.exited())
         BOOST_LOG_TRIVIAL(debug)<< "Program returned exit code " << s.exit_status() << std::endl;
         else if (s.signaled())
         {
         BOOST_LOG_TRIVIAL(debug) << "Program received signal " << s.term_signal() << std::endl;
         if (s.dumped_core())
         BOOST_LOG_TRIVIAL(debug) << "Program also dumped core" << std::endl;
         }
         else if (s.stopped())
         BOOST_LOG_TRIVIAL(debug) << "Program stopped by signal " << s.stop_signal() << std::endl;
         else
         BOOST_LOG_TRIVIAL(debug) << "Unknown termination reason" << std::endl;*/

        //it is important to include the std:: namespace specification
        //for the system call, otherwise the GNU compiler picks up
        //a version from the c library that gives trouble in threaded environments
        int result = std::system(execname.c_str());
        if (result)
          throw FatalException(
              "Cannot execute run script: " + runname + " Error code: "
                  + jif3D::stringify(result), __FILE__, __LINE__);
      }

    jif3D::rvec AdaptDist(const std::vector<double> &C, const jif3D::rvec &RawImpedance,
        const jif3D::rvec &Misfit, const std::vector<double> &RotAngles)
      {
        jif3D::rvec result(C.size(), 0.0);
        const size_t nstat = C.size() / 4;
        const size_t nfreq = RawImpedance.size() / (nstat * 16);
        //This is an implementation of equation 10 in Avdeeva et al. 2015
        //we calculate the partial derivative of the objective function
        //with respect to the distortion parameters

        for (size_t i = 0; i < nfreq; ++i)
          {
            for (size_t j = 0; j < nstat; ++j)
              {
                const double angle = RotAngles.at(j) / 180.0
                    * boost::math::constants::pi<double>();
                const double c2 = jif3D::pow2(std::cos(angle));
                const double s2 = jif3D::pow2(std::sin(angle));
                const double sc = std::sin(angle) * cos(angle);

//                result(j * 4) += Misfit(offset) * RawImpedance(offset)
//                    + Misfit(offset + 1) * RawImpedance(offset + 1)
//                    + Misfit(offset + 2) * RawImpedance(offset + 2)
//                    + Misfit(offset + 3) * RawImpedance(offset + 3);
//                result(j * 4 + 1) += Misfit(offset) * RawImpedance(offset + 4)
//                    + Misfit(offset + 1) * RawImpedance(offset + 5)
//                    + Misfit(offset + 2) * RawImpedance(offset + 6)
//                    + Misfit(offset + 3) * RawImpedance(offset + 7);
//                result(j * 4 + 2) += Misfit(offset + 4) * RawImpedance(offset)
//                    + Misfit(offset + 5) * RawImpedance(offset + 1)
//                    + Misfit(offset + 6) * RawImpedance(offset + 2)
//                    + Misfit(offset + 7) * RawImpedance(offset + 3);
//                result(j * 4 + 3) += Misfit(offset + 4) * RawImpedance(offset + 4)
//                    + Misfit(offset + 5) * RawImpedance(offset + 5)
//                    + Misfit(offset + 6) * RawImpedance(offset + 6)
//                    + Misfit(offset + 7) * RawImpedance(offset + 7);

                const size_t moffset = (i * nstat + j) * 8;
                const size_t roffset = moffset * 2;
                result(j * 4) +=
                    Misfit(moffset)
                        * (c2 * RawImpedance(roffset) + sc * RawImpedance(roffset + 4))
                        + Misfit(moffset + 1)
                            * (c2 * RawImpedance(roffset + 1)
                                + sc * RawImpedance(roffset + 5))
                        + Misfit(moffset + 2)
                            * (c2 * RawImpedance(roffset + 2)
                                + sc * RawImpedance(roffset + 6))
                        + Misfit(moffset + 3)
                            * (c2 * RawImpedance(roffset + 3)
                                + sc * RawImpedance(roffset + 7))
                        + Misfit(moffset + 4)
                            * (sc * RawImpedance(roffset) + s2 * RawImpedance(roffset + 4))
                        + Misfit(moffset + 5)
                            * (sc * RawImpedance(roffset + 1)
                                + s2 * RawImpedance(roffset + 5))
                        + Misfit(moffset + 6)
                            * (sc * RawImpedance(roffset + 2)
                                + s2 * RawImpedance(roffset + 6))
                        + Misfit(moffset + 7)
                            * (sc * RawImpedance(roffset + 3)
                                + s2 * RawImpedance(roffset + 7));

                result(j * 4 + 1) += Misfit(moffset)
                    * (-sc * RawImpedance(roffset) + c2 * RawImpedance(roffset + 4))
                    + Misfit(moffset + 1)
                        * (-sc * RawImpedance(roffset + 1)
                            + c2 * RawImpedance(roffset + 5))
                    + Misfit(moffset + 2)
                        * (-sc * RawImpedance(roffset + 2)
                            + c2 * RawImpedance(roffset + 6))
                    + Misfit(moffset + 3)
                        * (-sc * RawImpedance(roffset + 3)
                            + c2 * RawImpedance(roffset + 7))
                    + Misfit(moffset + 4)
                        * (-s2 * RawImpedance(roffset) + sc * RawImpedance(roffset + 4))
                    + Misfit(moffset + 5)
                        * (-s2 * RawImpedance(roffset + 1)
                            + sc * RawImpedance(roffset + 5))
                    + Misfit(moffset + 6)
                        * (-s2 * RawImpedance(roffset + 2)
                            + sc * RawImpedance(roffset + 6))
                    + Misfit(moffset + 7)
                        * (-s2 * RawImpedance(roffset + 3)
                            + sc * RawImpedance(roffset + 7));

                result(j * 4 + 2) += Misfit(moffset)
                    * (-sc * RawImpedance(roffset + 8) - s2 * RawImpedance(roffset + 12))
                    + Misfit(moffset + 1)
                        * (-sc * RawImpedance(roffset + 9)
                            - s2 * RawImpedance(roffset + 13))
                    + Misfit(moffset + 2)
                        * (-sc * RawImpedance(roffset + 10)
                            - s2 * RawImpedance(roffset + 14))
                    + Misfit(moffset + 3)
                        * (-sc * RawImpedance(roffset + 11)
                            - s2 * RawImpedance(roffset + 15))
                    + Misfit(moffset + 4)
                        * (c2 * RawImpedance(roffset + 8)
                            + sc * RawImpedance(roffset + 12))
                    + Misfit(moffset + 5)
                        * (c2 * RawImpedance(roffset + 9)
                            + sc * RawImpedance(roffset + 13))
                    + Misfit(moffset + 6)
                        * (c2 * RawImpedance(roffset + 10)
                            + sc * RawImpedance(roffset + 14))
                    + Misfit(moffset + 7)
                        * (c2 * RawImpedance(roffset + 11)
                            + sc * RawImpedance(roffset + 15));
                result(j * 4 + 3) += Misfit(moffset)
                    * (s2 * RawImpedance(roffset + 8) - sc * RawImpedance(roffset + 12))
                    + Misfit(moffset + 1)
                        * (s2 * RawImpedance(roffset + 9)
                            - sc * RawImpedance(roffset + 13))
                    + Misfit(moffset + 2)
                        * (s2 * RawImpedance(roffset + 10)
                            - sc * RawImpedance(roffset + 14))
                    + Misfit(moffset + 3)
                        * (s2 * RawImpedance(roffset + 11)
                            - sc * RawImpedance(roffset + 15))
                    + Misfit(moffset + 4)
                        * (-sc * RawImpedance(roffset + 8)
                            + c2 * RawImpedance(roffset + 12))
                    + Misfit(moffset + 5)
                        * (-sc * RawImpedance(roffset + 9)
                            + c2 * RawImpedance(roffset + 13))
                    + Misfit(moffset + 6)
                        * (-sc * RawImpedance(roffset + 10)
                            + c2 * RawImpedance(roffset + 14))
                    + Misfit(moffset + 7)
                        * (-sc * RawImpedance(roffset + 11)
                            + c2 * RawImpedance(roffset + 15));

              }
          }
        return result;
      }

    cmat CalcEExt(const rvec &Misfit, const std::vector<double> &C,
        const size_t startindex, const size_t freq_start_index,
        const std::complex<double> &Hx1, const std::complex<double> &Hx2,
        const std::complex<double> &Hy1, const std::complex<double> &Hy2)
      {
        const size_t siteindex = freq_start_index + startindex * 8;
        cmat AH(2, 2);
        const std::complex<double> magdet = 1. / (Hx1 * Hy2 - Hx2 * Hy1);
        const std::complex<double> A00(Misfit(siteindex), Misfit(siteindex + 1));
        const std::complex<double> A01(Misfit(siteindex + 2), Misfit(siteindex + 3));
        const std::complex<double> A10(Misfit(siteindex + 4), Misfit(siteindex + 5));
        const std::complex<double> A11(Misfit(siteindex + 6), Misfit(siteindex + 7));
        //x3d uses a different convention for the complex exponentials
        //so we have to use the complex conjugate for the source
        AH(0, 0) = magdet * (conj(A00) * Hy2 - conj(A01) * Hx2);
        AH(0, 1) = magdet * (-conj(A00) * Hy1 + conj(A01) * Hx1);
        AH(1, 0) = magdet * (conj(A10) * Hy2 - conj(A11) * Hx2);
        AH(1, 1) = magdet * (-conj(A10) * Hy1 + conj(A11) * Hx1);

        return AH;
      }

    void CalcU(const std::string &RootName,
        const std::vector<std::complex<double> > &XPolMoments,
        const std::vector<std::complex<double> > &YPolMoments,
        const std::vector<std::complex<double> > &ZPolMoments,
        std::vector<std::complex<double> > &Ux, std::vector<std::complex<double> > &Uy,
        std::vector<std::complex<double> > &Uz, const std::vector<size_t> &XSourceXIndex,
        const std::vector<size_t> &XSourceYIndex,
        const std::vector<double> &XSourceDepths,
        const std::vector<size_t> &YSourceXIndex,
        const std::vector<size_t> &YSourceYIndex,
        const std::vector<double> &YSourceDepths,
        const std::vector<size_t> &ZSourceXIndex,
        const std::vector<size_t> &ZSourceYIndex,
        const std::vector<double> &ZSourceDepths,
        const jif3D::ThreeDModelBase::t3DModelDim &ZCellBoundaries,
        const jif3D::ThreeDModelBase::t3DModelDim &ZCellSizes, const size_t ncellsx,
        const size_t ncellsy, const size_t ncellsz)
      {
        auto CompFunc = [](std::complex<double> a,std::complex<double> b ) -> bool
          { return std::norm(a) < std::norm(b);};
        auto XPolMax = std::max_element(XPolMoments.begin(), XPolMoments.end(), CompFunc);
        auto YPolMax = std::max_element(YPolMoments.begin(), YPolMoments.end(), CompFunc);
        auto ZPolMax = std::max_element(ZPolMoments.begin(), ZPolMoments.end(), CompFunc);
        bool HaveXPol = std::norm(*XPolMax) > 1e-30;
        bool HaveYPol = std::norm(*YPolMax) > 1e-30;
        bool HaveZPol = std::norm(*ZPolMax) > 1e-30;
        if (HaveXPol || HaveYPol || HaveZPol)
          {
            std::string DirName = RootName + dirext + "/";
#pragma omp critical(calcU_writesource)
              {
                WriteSourceFile(DirName + sourceafilename, XSourceXIndex, XSourceYIndex,
                    XSourceDepths, YSourceXIndex, YSourceYIndex, YSourceDepths,
                    ZSourceXIndex, ZSourceYIndex, ZSourceDepths, XPolMoments, YPolMoments,
                    ZPolMoments, ZCellBoundaries, ZCellSizes, ncellsx, ncellsy);
              }
            RunX3D(RootName);
#pragma omp critical(calcU_readema)
              {
                ReadEMA(DirName + emaname, Ux, Uy, Uz, ncellsx, ncellsy, ncellsz);
              }
          }
        else
          {

            Ux.resize(ncellsx * ncellsy * ncellsz, 0.0);
            Uy.resize(ncellsx * ncellsy * ncellsz, 0.0);
            Uz.resize(ncellsx * ncellsy * ncellsz, 0.0);
          }
        //boost::filesystem::remove_all(emaname);
      }

  }
