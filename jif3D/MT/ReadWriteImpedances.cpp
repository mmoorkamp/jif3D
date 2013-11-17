//============================================================================
// Name        : ReadWriteImpedance.cpp
// Author      : Jul 13, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================

#include <fstream>
#include <iomanip>
#include <algorithm>
#include <iostream>
#include <boost/cast.hpp>
#include "../Global/FatalException.h"
#include "../Global/NumUtil.h"
#include "../Global/convert.h"
#include "../ModelBase/NetCDFModelTools.h"
#include "../Global/NetCDFTools.h"
#include "ReadWriteImpedances.h"
#include "MTEquations.h"

namespace jif3D
{
using namespace std;
const std::string FreqDimName = "Frequency";

//write one compoment of the impedance tensor to a netcdf file
//this is an internal helper function
void WriteImpedanceComp(NcFile &NetCDFFile, NcDim *StatNumDim, NcDim *FreqDim,
		const jif3D::rvec &Impedances, const std::string &CompName,
		const size_t compindex)
{
	NcVar *CompVar = NetCDFFile.add_var(CompName.c_str(), ncDouble, FreqDim,
			StatNumDim);
	jif3D::rvec Component(FreqDim->size() * StatNumDim->size());
	for (size_t i = 0; i < Component.size(); ++i)
	{
		Component(i) = Impedances(i * 8 + compindex);
	}
	CompVar->put(&Component[0], FreqDim->size(), StatNumDim->size());
	//we write the impedance in units of Ohm
	CompVar->add_att("units", "Ohm");
}

//read one component of the impedance tensor from a netcdf file
//this is an internal helper function
void ReadImpedanceComp(NcFile &NetCDFFile, jif3D::rvec &Impedances,
		const std::string &CompName, const size_t compindex,
		const bool MustExist = true)
{
	NcError MyError(
			MustExist ? NcError::verbose_fatal : NcError::silent_nonfatal);
	NcVar *SizeVar;
	if ((SizeVar = NetCDFFile.get_var(CompName.c_str())))
	{
		const size_t nvalues = Impedances.size() / 8;
		assert(
				nvalues == boost::numeric_cast<size_t>( SizeVar->edges()[0] * SizeVar->edges()[1]));
		jif3D::rvec Temp(nvalues);
		SizeVar->get(&Temp[0], SizeVar->edges()[0], SizeVar->edges()[1]);

		for (size_t i = 0; i < nvalues; ++i)
		{
			Impedances(i * 8 + compindex) = Temp(i);
		}
	}
}

void WriteImpedancesToNetCDF(const std::string &filename,
		const std::vector<double> &Frequencies,
		const std::vector<double> &StatXCoord,
		const std::vector<double> &StatYCoord,
		const std::vector<double> &StatZCoord, const jif3D::rvec &Impedances,
		const jif3D::rvec &Errors)
{
	const size_t nstats = StatXCoord.size();
	const size_t nfreqs = Frequencies.size();
	const size_t nimp = nstats * nfreqs * 8;
	assert(nstats == StatYCoord.size());
	assert(nstats == StatYCoord.size());
	assert(Impedances.size() == nimp);
	//create a netcdf file
	NcFile DataFile(filename.c_str(), NcFile::Replace);
	//each station gets an index number that we use as a dimension
	//in the netcdf file
	NcDim *StatNumDim = DataFile.add_dim(StationNumberName.c_str(), nstats);
	std::vector<int> StationNumber;
	std::generate_n(back_inserter(StationNumber), nstats, IntSequence(0));
	NcVar *StatNumVar = DataFile.add_var(StationNumberName.c_str(), ncInt,
			StatNumDim);
	StatNumVar->put(&StationNumber[0], nstats);
	//write out the measurement coordinates
	WriteVec(DataFile, MeasPosXName, StatXCoord, StatNumDim, "m");
	WriteVec(DataFile, MeasPosYName, StatYCoord, StatNumDim, "m");
	WriteVec(DataFile, MeasPosZName, StatZCoord, StatNumDim, "m");
	//write out the frequencies that we store
	NcDim *FreqDim = DataFile.add_dim(FreqDimName.c_str(), Frequencies.size());
	NcVar *FreqVar = DataFile.add_var(FreqDimName.c_str(), ncDouble, FreqDim);
	FreqVar->put(&Frequencies[0], nfreqs);
	//and now we can write all the impedance components
	WriteImpedanceComp(DataFile, StatNumDim, FreqDim, Impedances, "Zxx_re", 0);
	WriteImpedanceComp(DataFile, StatNumDim, FreqDim, Impedances, "Zxx_im", 1);
	WriteImpedanceComp(DataFile, StatNumDim, FreqDim, Impedances, "Zxy_re", 2);
	WriteImpedanceComp(DataFile, StatNumDim, FreqDim, Impedances, "Zxy_im", 3);
	WriteImpedanceComp(DataFile, StatNumDim, FreqDim, Impedances, "Zyx_re", 4);
	WriteImpedanceComp(DataFile, StatNumDim, FreqDim, Impedances, "Zyx_im", 5);
	WriteImpedanceComp(DataFile, StatNumDim, FreqDim, Impedances, "Zyy_re", 6);
	WriteImpedanceComp(DataFile, StatNumDim, FreqDim, Impedances, "Zyy_im", 7);
	//now we deal with the errors, if no parameter has been explicitly passed
	//Errors is empty, so we just fill the vector with zeros
	jif3D::rvec ZErr(Errors);
	if (ZErr.empty())
	{
		ZErr.resize(nimp);
		std::fill_n(ZErr.begin(), nimp, 0.0);
	}
	//in any case we write some error information even if it is all zero
	//as the error is identical for real and imaginary part
	//we only have to write the elements corresponding to the real part
	WriteImpedanceComp(DataFile, StatNumDim, FreqDim, ZErr, "dZxx", 0);
	WriteImpedanceComp(DataFile, StatNumDim, FreqDim, ZErr, "dZxy", 2);
	WriteImpedanceComp(DataFile, StatNumDim, FreqDim, ZErr, "dZyx", 4);
	WriteImpedanceComp(DataFile, StatNumDim, FreqDim, ZErr, "dZyy", 6);
}

void ReadImpedancesFromNetCDF(const std::string &filename,
		std::vector<double> &Frequencies, std::vector<double> &StatXCoord,
		std::vector<double> &StatYCoord, std::vector<double> &StatZCoord,
		jif3D::rvec &Impedances, jif3D::rvec &ImpError)
{
	//open the netcdf file readonly
	NcFile DataFile(filename.c_str(), NcFile::ReadOnly);
	//read in the station coordinates in the three directions
	ReadVec(DataFile, MeasPosXName, StatXCoord);
	ReadVec(DataFile, MeasPosYName, StatYCoord);
	ReadVec(DataFile, MeasPosZName, StatZCoord);
	//read in the frequencies
	ReadVec(DataFile, FreqDimName, Frequencies);
	//for each frequency and each station we have 8 impedances
	//currently it is not possible to have a different number of frequencies
	//at each site
	const size_t nimp = Frequencies.size() * StatXCoord.size() * 8;
	Impedances.resize(nimp);
	ImpError.resize(nimp);
	//read the impedances
	ReadImpedanceComp(DataFile, Impedances, "Zxx_re", 0);
	ReadImpedanceComp(DataFile, Impedances, "Zxx_im", 1);
	ReadImpedanceComp(DataFile, Impedances, "Zxy_re", 2);
	ReadImpedanceComp(DataFile, Impedances, "Zxy_im", 3);
	ReadImpedanceComp(DataFile, Impedances, "Zyx_re", 4);
	ReadImpedanceComp(DataFile, Impedances, "Zyx_im", 5);
	ReadImpedanceComp(DataFile, Impedances, "Zyy_re", 6);
	ReadImpedanceComp(DataFile, Impedances, "Zyy_im", 7);
	//now read the errors, we make their presence in the file optional
	ReadImpedanceComp(DataFile, ImpError, "dZxx", 0, false);
	ReadImpedanceComp(DataFile, ImpError, "dZxy", 2, false);
	ReadImpedanceComp(DataFile, ImpError, "dZyx", 4, false);
	ReadImpedanceComp(DataFile, ImpError, "dZyy", 6, false);
	//we only read in the errors for the real part
	//and then assign the same error to the imaginary part
	for (size_t i = 0; i < nimp - 1; i += 2)
	{
		ImpError(i + 1) = ImpError(i);
	}
}

void ReadImpedancesFromMTT(const std::string &filename,
		std::vector<double> &Frequencies, jif3D::rvec &Impedances,
		jif3D::rvec &Errors)
{
	std::ifstream infile;
	double currentreal, currentimag;
	int nentries = 0;
	infile.open(filename.c_str());
	while (infile.good())
	{
		infile >> currentreal;
		++nentries;
	}
	infile.close();
	if (((nentries - 1) % 23) != 0)
		throw FatalException(
				"Number of records does not match expected: " + filename);
	const int nrecords = (nentries - 1) / 23;
	Impedances.resize(nrecords * 8);
	Errors.resize(nrecords * 8);
	Frequencies.resize(nrecords);
	infile.open(filename.c_str());
	int currentrecord = 0;
	if (infile)
	{
		while (infile.good()) // read in inputfile
		{
			infile >> Frequencies[currentrecord / 8];
			if (infile.good())
			{
				//read number of degrees of freedom in file an throw away
				infile >> currentreal;

				infile >> Impedances(currentrecord)
						>> Impedances(currentrecord + 1)
						>> Impedances(currentrecord + 2)
						>> Impedances(currentrecord + 3)
						>> Impedances(currentrecord + 4)
						>> Impedances(currentrecord + 5)
						>> Impedances(currentrecord + 6)
						>> Impedances(currentrecord + 7);
				// read in the impedance errors
				infile >> Errors(currentrecord) >> Errors(currentrecord + 2)
						>> Errors(currentrecord + 4)
						>> Errors(currentrecord + 6);
				//fpr the moment we ignore these values in the .mtt file
				//Tx
				infile >> currentreal >> currentimag;
				//Ty
				infile >> currentreal >> currentimag;
				//dTx
				infile >> currentreal;
				//dTy
				infile >> currentreal;
				//Coherence Rx
				infile >> currentreal;
				//Coherence Ry
				infile >> currentreal;
				//Coherence Rz
				infile >> currentreal;

				currentrecord += 8;
			}
		}
		infile.close();

		//we read the error information into the vector components corresponding to the real part
		//we have to copy that error to the components corresponding to the imaginary part
		for (size_t i = 0; i < Errors.size() - 1; i += 2)
		{
			Errors(i + 1) = Errors(i);
		}
		//convert the units in the .mtt file (km/s) into S.I units (Ohm)
		const double convfactor = 4.0 * 1e-4 * acos(-1.0);
		Impedances *= convfactor;
		Errors *= convfactor;
	}
	else
	{
		throw jif3D::FatalException("File not found: " + filename);
	}
}

void WriteImpedancesToMtt(const std::string &filenamebase,
		const std::vector<double> &Frequencies, const jif3D::rvec &Imp,
		const jif3D::rvec &Err)
{
	const double convfactor = 4.0 * 1e-4 * acos(-1.0);
	jif3D::rvec Impedances = 1.0 / convfactor * Imp;
	jif3D::rvec Errors = 1.0 / convfactor * Err;
	const size_t nfreq = Frequencies.size();
	const size_t nimp = Impedances.size();
	const size_t ndatapersite = nfreq * 8;
	const size_t nsites = nimp / ndatapersite;
	assert(nimp % ndatapersite == 0);
	for (size_t i = 0; i < nsites; ++i)
	{

		ofstream outfile;
		std::string currfilename = filenamebase + jif3D::stringify(i) + ".mtt";
		outfile.open(currfilename.c_str());

		for (unsigned int j = 0; j < nfreq; ++j) //write mtt-file
		{
			const size_t startindex = (j * nsites + i) * 8;
			outfile << Frequencies.at(j);
			outfile << "   1 \n";

			outfile << Impedances(startindex) << " ";
			outfile << Impedances(startindex + 1) << " ";
			outfile << Impedances(startindex + 2) << " ";
			outfile << Impedances(startindex + 3) << " ";
			outfile << Impedances(startindex + 4) << " ";
			outfile << Impedances(startindex + 5) << " ";
			outfile << Impedances(startindex + 6) << " ";
			outfile << Impedances(startindex + 7) << " ";
			outfile << "\n";
			outfile << Errors(startindex) << " ";
			outfile << Errors(startindex + 2) << " ";
			outfile << Errors(startindex + 4) << " ";
			outfile << Errors(startindex + 6) << " ";
			outfile << 0.0 << " ";
			outfile << 0.0 << " ";
			outfile << 0.0 << " ";
			outfile << 0.0 << " ";
			outfile << "\n";
			outfile << 0.0 << " ";
			outfile << 0.0 << " ";
			outfile << 0.0 << " ";
			outfile << 0.0 << " ";
			outfile << 0.0 << " ";
			outfile << "\n";
			//write out mtt file entries

		}
		outfile.close();
	}
}

inline double rad(double d)
{
	return d / 180.0 * boost::math::constants::pi<double>();
}

void ReadAppResFromAscii(const std::string &filename,
		std::vector<double> &Frequencies, std::vector<double> &StatXCoord,
		std::vector<double> &StatYCoord, std::vector<double> &StatZCoord,
		jif3D::rvec &Imp, jif3D::rvec &Err)
{
	std::ifstream infile;
	infile.open(filename.c_str());
	char dummy[5000];
	//skip header line
	infile.getline(dummy, 5000);
	std::string major, link, site;
	double x, y, elevation, fre, rxy, ryx, pxy, pyx, erxy, eryx, epxy, epyx;

	std::vector<double> ImpTemp, ErrTemp;
	while (infile.good())
	{
		infile >> major >> link >> site >> x >> y >> elevation >> fre >> rxy
				>> ryx >> pxy >> pyx >> erxy >> eryx >> epxy >> epyx;
		if (infile.good())
		{
			//std::cout << fre << " " << rxy << " " << pxy << std::endl;
			ImpTemp.push_back(0.0);
			ImpTemp.push_back(0.0);
			fre = std::pow(10, fre);
			const double absZxy = sqrt(twopimu * fre * rxy);
			const double cpxy = cos(rad(pxy));
			const double spxy = sin(rad(pxy));
			const double cpyx = cos(rad(pyx));
			const double spyx = sin(rad(pyx));
			ImpTemp.push_back(absZxy * cpxy);
			ImpTemp.push_back(absZxy * spxy);
			//std::cout << "Fre: " << fre << " Rho: " << AppRes(std::complex<double>(absZxy * cpxy,absZxy * spxy),fre) << std::endl;
			//std::cout << "Zxy_re: " << absZxy * cpxy << " Zxy_im " << absZxy * spxy << std::endl;
			const double absZyx = sqrt(twopimu * fre * ryx);
			ImpTemp.push_back(absZyx * cpyx);
			ImpTemp.push_back(absZyx * spyx);
			ImpTemp.push_back(0.0);
			ImpTemp.push_back(0.0);
			ErrTemp.push_back(0.0);
			ErrTemp.push_back(0.0);
			double dR = sqrt(
					twopimu * fre / (4 * rxy) * std::pow(erxy * cpxy, 2)
							+ std::pow(absZxy * spxy * rad(epxy), 2));
			double dI = sqrt(
					twopimu * fre / (4 * rxy) * std::pow(erxy * spxy, 2)
							+ std::pow(absZxy * cpxy * rad(epxy), 2));
			ErrTemp.push_back(std::max(dR, dI));
			ErrTemp.push_back(std::max(dR, dI));
			dR = sqrt(
					twopimu * fre / (4 * ryx) * std::pow(eryx * cpxy, 2)
							+ std::pow(absZyx * spyx * rad(epyx), 2));
			dI = sqrt(
					twopimu * fre / (4 * ryx) * std::pow(eryx * spxy, 2)
							+ std::pow(absZyx * cpyx * rad(epyx), 2));
			ErrTemp.push_back(std::max(dR, dI));
			ErrTemp.push_back(std::max(dR, dI));
			ErrTemp.push_back(0.0);
			ErrTemp.push_back(0.0);
			Frequencies.push_back(fre);
			StatXCoord.push_back(x);
			StatYCoord.push_back(y);
			StatZCoord.push_back(elevation);
		}

	}

	Imp.resize(ImpTemp.size());
	Err.resize(ErrTemp.size());
	//std::copy(ImpTemp.begin(), ImpTemp.end(), Imp.begin());
	//std::copy(ErrTemp.begin(), ErrTemp.end(), Err.begin());
	std::sort(Frequencies.begin(), Frequencies.end(), std::greater<double>());
	Frequencies.erase(
			std::unique(Frequencies.begin(), Frequencies.end(),
					[](double a, double b)
					{	return std::abs((a-b)/std::max(a,b)) < 0.001;}),
			Frequencies.end());
	StatXCoord.erase(std::unique(StatXCoord.begin(), StatXCoord.end()),
			StatXCoord.end());
	StatYCoord.erase(std::unique(StatYCoord.begin(), StatYCoord.end()),
			StatYCoord.end());
	StatZCoord.erase(std::unique(StatZCoord.begin(), StatZCoord.end()),
			StatZCoord.end());
	const size_t nfreq = Frequencies.size();
	const size_t nstat = StatXCoord.size();
	for (size_t i = 0; i < nstat; ++i)
	{
		for (size_t j = 0; j < nfreq; ++j)
		{
			size_t Impindex = (j * nstat + i) * 8;
			size_t TempIndex = (i * nfreq + j) * 8;
			std::copy(ImpTemp.begin() + TempIndex,
					ImpTemp.begin() + TempIndex + 8, Imp.begin() + Impindex);
			std::copy(ErrTemp.begin() + TempIndex,
					ErrTemp.begin() + TempIndex + 8, Err.begin() + Impindex);
		}
	}
	std::cout << "Stations: " << StatXCoord.size() << std::endl;
	std::cout << "Frequencies: " << Frequencies.size() << std::endl;
	std::copy(Frequencies.begin(), Frequencies.end(),
			std::ostream_iterator<double>(std::cout, " "));
	std::cout << std::endl;
	std::cout << "Impendances: " << Imp.size() << std::endl;
}

}
