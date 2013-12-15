#include "../Global/FileUtil.h"
#include "../MT/ReadWriteImpedances.h"

double AvgComp(size_t compindex, size_t nstat, size_t freqindex,
		const jif3D::rvec &Imp, const jif3D::rvec &Err)
{
	double errthresh = 0.5;
	double avg = 0.0;
	size_t n = 0;
	for (size_t j = 0; j < nstat; ++j)
	{
		double currimp = Imp((freqindex * nstat + j) * 8 + compindex);
		double currerr = Err((freqindex * nstat + j) * 8 + compindex);
		if (currerr < errthresh * std::abs(currimp))
		{
			avg += currimp;
			++n;
		}
	}
	std::cout << "F: " << freqindex << " Comp: " << compindex << " N: " << n
			<< std::endl;
	if (n > 0)
	{
		return avg / n;
	}
	return 0.0;
}

int main()
{
	std::string filename = jif3D::AskFilename("Impedance file: ");

	std::vector<double> Frequencies, XCoord, YCoord, ZCoord,C;
	jif3D::rvec Imp, Err;
	jif3D::ReadImpedancesFromNetCDF(filename, Frequencies, XCoord, YCoord,
			ZCoord, Imp, Err, C);

	size_t nstat = XCoord.size();
	size_t nfreq = Frequencies.size();
	jif3D::rvec AvgImp(8 * nfreq, 0.0);

	for (size_t i = 0; i < nfreq; ++i)
	{
		for (size_t j = 0; j < 8; ++j)
		{
			AvgImp(i * 8 + j) = AvgComp(j, nstat, i, Imp, Err);
		}
	}
	jif3D::rvec AvgErr(8 * nfreq, 0.0);
	std::vector<double> AvgXCoord(1, 0.0), AvgYCoord(1, 0.0), AvgZCoord(1, 0.0);
	jif3D::WriteImpedancesToNetCDF(filename + ".avg.nc", Frequencies, AvgXCoord,
			AvgYCoord, AvgZCoord, AvgImp, AvgErr);
	jif3D::WriteImpedancesToMtt(filename + ".avg.nc", Frequencies, AvgImp,
			AvgErr);

}
