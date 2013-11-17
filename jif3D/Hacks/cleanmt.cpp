#include "../Global/FileUtil.h"
#include <iostream>
#include "../MT/ReadWriteImpedances.h"
#include "../Global/convert.h"
#include "../ModelBase/VTKTools.h"


int main()
{
	std::string filename = jif3D::AskFilename("MT data:");
	std::vector<double> Freq, StatX, StatY, StatZ;
	jif3D::rvec Imp, Err;
	jif3D::ReadImpedancesFromNetCDF(filename, Freq, StatX, StatY, StatZ, Imp,
			Err);
	const size_t nstat = StatX.size();
	const size_t nfreq = Freq.size();
	double lowquant, highquant;
	std::cout << "Lower quantile: ";
	std::cin >> lowquant;
	std::cout << "Upper quantile: ";
		std::cin >> highquant;
	for (size_t i = 0; i < nfreq; ++i)
	{
		for (size_t k = 2; k < 6; ++k)
		{
			std::vector<double> CurrValues;
			for (size_t j = i * nstat * 8; j < (i + 1) * nstat * 8; j += 8)
			{
				CurrValues.push_back(Imp(j + k));
			}
			std::sort(CurrValues.begin(), CurrValues.end());
			size_t lowindex = std::max(0, int(lowquant * CurrValues.size()-1));
			size_t highindex = highquant * CurrValues.size()-1;
			size_t medindex = CurrValues.size()/2;
			double minthresh = CurrValues.at(lowindex);
			double maxthresh = CurrValues.at(highindex);
			double median = CurrValues.at(medindex);
			std::cout << "Frequency: " << Freq.at(i) << " " << minthresh << " "
					<< median << " " << maxthresh << std::endl;
			std::string filename = "Curr_" + jif3D::stringify(Freq.at(i)) + "_" + jif3D::stringify(k) + ".out";
			std::ofstream currfile;
			currfile.open(filename.c_str());
			std::copy(CurrValues.begin(), CurrValues.end(),
					std::ostream_iterator<double>(currfile, "\n"));
			currfile << "\n  \n 0 " << minthresh << "\n  " << CurrValues.size()
					<< "  " << minthresh << "\n";
			currfile << "\n  \n 0 " << maxthresh << "\n  " << CurrValues.size()
					<< "  " << maxthresh << "\n";
			currfile << "\n  \n 0 " << median << "\n  " << CurrValues.size()
					<< "  " << median << "\n";
			currfile.close();
//			double minselect, maxselect;
//			std::cout << "Minselect: ";
//			std::cin >> minselect;
//			std::cout << "Maxselect: ";
//			std::cin >> maxselect;
			jif3D::rvec CurrImp(StatX.size()), CurrErr(StatX.size());
			size_t currindex = 0;
			for (size_t j = i * nstat * 8; j < (i + 1) * nstat * 8; j += 8)
			{
				CurrImp(currindex) = Imp(j+k);
				CurrErr(currindex) = Err(j+k);
				if (Imp(j + k) < minthresh || Imp(j + k) > maxthresh)
				{
					Err(j + k) = std::max(std::abs(2*Imp(j + k)),Err(j + k));
				}
				currindex++;
			}
			std::string impfilename = "Imp_" + jif3D::stringify(Freq.at(i)) + "_" + jif3D::stringify(k) + ".vtk";
			std::string errfilename = "Err_" + jif3D::stringify(Freq.at(i)) + "_" + jif3D::stringify(k) + ".vtk";
			jif3D::Write3DDataToVTK(impfilename,"Z",CurrImp,StatX,StatY,StatZ);
			jif3D::Write3DDataToVTK(errfilename,"dZ",CurrErr,StatX,StatY,StatZ);
		}
		std::cout << std::endl;
	}
	jif3D::WriteImpedancesToNetCDF(filename + ".clean.nc", Freq, StatX, StatY,
			StatZ, Imp, Err);
}
