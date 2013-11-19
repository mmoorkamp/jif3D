//============================================================================
// Name        : CalcFreq.cpp
// Author      : 23 Feb 2013
// Version     : 
// Copyright   : 2013, mm489
//============================================================================

#include "CalcFreq.h"
#include <boost/filesystem/operations.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <map>
#include <fstream>
#include <boost/assign/list_of.hpp> // for 'map_list_of()'
#include "../Global/convert.h"
#include "../Global/FatalException.h"
#include "../MT/ReadWriteX3D.h"
#include "../MT/MTEquations.h"
#include "../ModelBase/CellBoundaries.h"


namespace fs = boost::filesystem;

//define some names that are always the same
//either because x3d uses this convention
//or because we use them in their own directory
//and want to keep them simple to make sure x3d can handle them
const std::string modelfilename("x3d.model");
const std::string resultfilename("x3d.result");
const std::string emaname = resultfilename + "0.ema";
const std::string sourceafilename = modelfilename + "0a.source";
const std::string sourcebfilename = modelfilename + "0b.source";
const std::string emoAname = resultfilename + "0a.emo";
const std::string emoBname = resultfilename + "0b.emo";
const std::string emaAname = resultfilename + "0a.ema";
const std::string emaBname = resultfilename + "0b.ema";
const std::string runext = "_run";
const std::string dirext = "_dir";
//associate a type of calculation with a name string
const std::map<jif3D::X3DModel::ProblemType, std::string> Extension =
    boost::assign::map_list_of(jif3D::X3DModel::MT, "MT")(jif3D::X3DModel::EDIP, "EDIP")(
        jif3D::X3DModel::MDIP, "MDIP");



inline void CheckField(const std::vector<std::complex<double> > &Field, size_t nelem)
  {
    if (Field.size() != nelem)
      throw jif3D::FatalException(
          "Number of read in elements in Field: " + jif3D::stringify(Field.size())
              + " does not match expected: " + jif3D::stringify(nelem));
  }

//check that the .hnk file for x3d are in a certain directory
bool CheckHNK(const fs::path &TargetDir)
  {
    return fs::exists(TargetDir / "ndec15.hnk") && fs::exists(TargetDir / "ndec20.hnk")
        && fs::exists(TargetDir / "ndec30.hnk") && fs::exists(TargetDir / "ndec40.hnk");
  }

//copy the .hnk files for x3d from SourceDir to TargetDir
void CopyHNK(const fs::path &SourceDir, const fs::path &TargetDir)
  {
    //copy file fails with an exception if the target exists
    //so we check before we do the actual copy
    if (!CheckHNK(TargetDir))
      {
        fs::copy(SourceDir / "ndec15.hnk", TargetDir / "ndec15.hnk");
        fs::copy(SourceDir / "ndec20.hnk", TargetDir / "ndec20.hnk");
        fs::copy(SourceDir / "ndec30.hnk", TargetDir / "ndec30.hnk");
        fs::copy(SourceDir / "ndec40.hnk", TargetDir / "ndec40.hnk");
      }
  }

//create a unique ID that we can use to name things and still
//perform parallel calculations
std::string ObjectID()
  {
    //a unique ID created on construction
    boost::uuids::uuid tag = boost::uuids::random_generator()();
    //make a unique filename for the sensitivity file created by this object
    //we use boost uuid to generate a unique identifier tag
    //and translate it to a string to generate the filename
    return "mt" + jif3D::stringify(getpid()) + jif3D::stringify(tag);
  }

std::string MakeUniqueName(jif3D::X3DModel::ProblemType Type, const size_t FreqIndex)
  {
    //we assemble the name from the id of the process
    //and the address of the current object
    std::string result(ObjectID());
    //the type of calculation
    result += Extension.find(Type)->second;
    //and the frequency index
    result += jif3D::stringify(FreqIndex);
    return result;
  }
//create a script that changes to the correct directory
//and executes x3d in that directory
void MakeRunFile(const std::string &NameRoot, const std::string DirName)
  {
    std::string RunFileName = NameRoot + runext;
    fs::create_directory(DirName);
    std::ofstream runfile;
    runfile.open(RunFileName.c_str());
    runfile << "#!/bin/bash\n";
    runfile << "cd " << DirName << "\n";
    runfile << "x3d > /dev/null\n";
    runfile << "cd ..\n";
    runfile.close();
    //we also copy the necessary *.hnk files
    //from the current directory to the work directory
    CopyHNK(fs::current_path(), DirName);
  }

//execute the script that runs x3d
void RunX3D(const std::string &NameRoot)
  {
    std::cout << "Running x3D from hpx " << std::endl;
	//instead of making the script executable
    //we run a bash with the scriptname as an argument
    //this turns out to be more robust

    const std::string runname = "bash " + NameRoot + runext;
    //it is important to include the std:: namespace specification
    //for the system call, otherwise the GNU compiler picks up
    //a version from the c library that gives trouble in threaded environments
    if (std::system(runname.c_str()))
      throw jif3D::FatalException("Cannot execute run script: " + runname);
  }

jif3D::rvec HPXCalculateFrequency(const jif3D::X3DModel &Model, size_t freqindex,
    std::string TempDirName)
  {

	const size_t nmeas = Model.GetMeasPosX().size();
	        const size_t nmodx = Model.GetXCoordinates().size();
	        const size_t nmody = Model.GetYCoordinates().size();
	        jif3D::rvec result(nmeas * 8);
	        fs::path RootName =  MakeUniqueName(jif3D::X3DModel::MT, freqindex);
	        fs::path DirName = RootName.string() + dirext;
	        std::vector<double> CurrFreq(1, Model.GetFrequencies()[freqindex]);
	        std::vector<double> ShiftDepth;
	        std::vector<size_t> MeasDepthIndices;
	        //construct a vector of indices of unique station depths
	        size_t nlevels = ConstructDepthIndices(MeasDepthIndices, ShiftDepth, Model);
	        //writing out files causes problems in parallel
	        // so we make sure it is done one at a time
	#pragma omp critical(forward_write_files)
	          {
	        	MakeRunFile(RootName.string(), DirName.string());
	        	jif3D::WriteProjectFile(DirName.string(), CurrFreq, jif3D::X3DModel::MT, resultfilename,
	                modelfilename);
	        	jif3D::Write3DModelForX3D((DirName / modelfilename).string(), Model.GetXCellSizes(),
	                Model.GetYCellSizes(), Model.GetZCellSizes(), ShiftDepth,
	                Model.GetConductivities(), Model.GetBackgroundConductivities(),
	                Model.GetBackgroundThicknesses());
	          }
	        //run x3d in parallel
	        RunX3D(RootName.string());
	        std::vector<std::complex<double> > Ex1, Ex2, Ey1, Ey2, Hx1, Hx2, Hy1, Hy2;
	        std::complex<double> Zxx, Zxy, Zyx, Zyy;
	        //read in the electric and magnetic field at the observe sites
	#pragma omp critical(forward_read_emo)
	          {
	        	jif3D::ReadEMO((DirName / emoAname).string(), Ex1, Ey1, Hx1, Hy1);
	        	jif3D::ReadEMO((DirName / emoBname).string(), Ex2, Ey2, Hx2, Hy2);
	          }
	        const size_t nval = (nmodx * nmody * nlevels);
	       CheckField(Ex1, nval);
	        CheckField(Ex2, nval);
	        CheckField(Ey1, nval);
	        CheckField(Ey2, nval);
	        CheckField(Hx1, nval);
	        CheckField(Hx2, nval);
	        CheckField(Hy1, nval);
	        CheckField(Hy2, nval);
	        //calculate impedances from the field spectra for all measurement sites
	        for (size_t j = 0; j < nmeas; ++j)
	          {
	            boost::array<jif3D::ThreeDModelBase::t3DModelData::index, 3> StationIndex =
	                Model.FindAssociatedIndices(Model.GetMeasPosX()[j],
	                    Model.GetMeasPosY()[j], Model.GetMeasPosZ()[j]);
	// with the current equations we cannot use interpolation as the position
	        	// of the source in the adjoint calculation is always smeared
	        	//across the whole cell, this works best for a cell in the centre
	        	//of the cell and any other position deteriorates convergence
	//            std::complex<double> Ex1Inter = InterpolateField(Ex1, Model, j,
	//                MeasDepthIndices);
	//            std::complex<double> Ex2Inter = InterpolateField(Ex2, Model, j,
	//                MeasDepthIndices);
	//            std::complex<double> Ey1Inter = InterpolateField(Ey1, Model, j,
	//                MeasDepthIndices);
	//            std::complex<double> Ey2Inter = InterpolateField(Ey2, Model, j,
	//                MeasDepthIndices);
	//            std::complex<double> Hx1Inter = InterpolateField(Hx1, Model, j,
	//                MeasDepthIndices);
	//            std::complex<double> Hx2Inter = InterpolateField(Hx2, Model, j,
	//                MeasDepthIndices);
	//            std::complex<double> Hy1Inter = InterpolateField(Hy1, Model, j,
	//                MeasDepthIndices);
	//            std::complex<double> Hy2Inter = InterpolateField(Hy2, Model, j,
	//                MeasDepthIndices);
	//            FieldsToImpedance(Ex1Inter, Ex2Inter, Ey1Inter, Ey2Inter, Hx1Inter, Hx2Inter,
	//                Hy1Inter, Hy2Inter, Zxx, Zxy, Zyx, Zyy);
	            const size_t offset = (nmodx * nmody) * MeasDepthIndices[j] + StationIndex[0] * nmody
	                + StationIndex[1];
	            jif3D::FieldsToImpedance(Ex1[offset], Ex2[offset], Ey1[offset], Ey2[offset],
	                Hx1[offset], Hx2[offset], Hy1[offset], Hy2[offset], Zxx, Zxy, Zyx, Zyy);
	            //result is a local array for this frequency
	            //so we can directly use it even in a threaded environment
	            const size_t meas_index = j * 8;
	            result(meas_index) = Zxx.real();
	            result(meas_index + 1) = Zxx.imag();
	            result(meas_index + 2) = Zxy.real();
	            result(meas_index + 3) = Zxy.imag();
	            result(meas_index + 4) = Zyx.real();
	            result(meas_index + 5) = Zyx.imag();
	            result(meas_index + 6) = Zyy.real();
	            result(meas_index + 7) = Zyy.imag();
	          }
	        return result;
  }

HPX_REGISTER_PLAIN_ACTION(CalculateFrequency_action)
