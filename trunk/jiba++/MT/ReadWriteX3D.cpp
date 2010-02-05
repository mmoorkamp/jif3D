//============================================================================
// Name        : ReadWriteX3D.cpp
// Author      : Feb 17, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#include <cassert>
#include <fstream>
#include <iomanip>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include "../Global/FileUtil.h"
#include "../Global/FatalException.h"
#include "../Global/convert.h"
#include "ReadWriteX3D.h"

namespace jiba
  {

    std::vector<std::complex<double> > ResortFields(const std::vector<
        std::complex<double> > &InField, const size_t nx, const size_t ny,
        const size_t nz);

    void Read3DModelFromX3D(const std::string &filename,
        ThreeDModelBase::t3DModelDim &XCellSizes,
        ThreeDModelBase::t3DModelDim &YCellSizes,
        ThreeDModelBase::t3DModelDim &ZCellSizes,
        ThreeDModelBase::t3DModelData &Data,
        std::vector<double> &bg_conductivities,
        std::vector<double> &bg_thicknesses)
      {
        std::ifstream infile(filename.c_str());
        std::string line = FindToken(infile, "Dx");
        double dx = 2.1, dy = 3.4;
        infile >> dx >> dy;

        line = FindToken(infile, "Thickness");
        double currvalue;
        //swallow the line with 0 thickness and conductivity
        char dummy[1024];
        infile.getline(dummy, 1024);
        while (!infile.fail())
          {
            infile >> currvalue;
            if (!infile.fail())
              bg_thicknesses.push_back(currvalue);
            infile >> currvalue;
            if (!infile.fail())
              {
                char dummy[1024];
                bg_conductivities.push_back(currvalue);
                infile.getline(dummy, 1024);
              }
          }
        infile.clear();
        std::vector<double> Zthick;
        std::vector<double> Values;
        bool havelayer = true;

        int startx, starty, endx, endy;
        while (havelayer)
          {
            try
              {
                line = FindToken(infile, "dzA(m)");
              } catch (FatalException &e)
              {
                havelayer = false;
              }
            if (havelayer)
              {
                infile >> currvalue;
                Zthick.push_back(currvalue);

                line = FindToken(infile, "cells_in_X-direction");
                infile >> startx >> endx;
                FindToken(infile, "cells_in_Y-direction");
                infile >> starty >> endy;
                FindToken(infile, "ARRAY");
                const unsigned int nelements = (endx - startx + 1) * (endy
                    - starty + 1);
                for (size_t i = 0; i < nelements; ++i)
                  {
                    infile >> currvalue;
                    Values.push_back(currvalue);
                  }
              }
          }//end of while
        const size_t nx = (endx - startx + 1);
        const size_t ny = (endy - starty + 1);
        const size_t nz = Zthick.size();
        assert(nx * ny * nz == Values.size());
        XCellSizes.resize(boost::extents[nx]);
        YCellSizes.resize(boost::extents[ny]);
        ZCellSizes.resize(boost::extents[nz]);
        std::copy(Zthick.begin(), Zthick.end(), ZCellSizes.begin());
        std::fill(XCellSizes.begin(), XCellSizes.end(), dx);
        std::fill(YCellSizes.begin(), YCellSizes.end(), dy);
        Data.resize(boost::extents[nx][ny][nz]);
        for (size_t i = 0; i < nx; ++i)
          for (size_t j = 0; j < ny; ++j)
            for (size_t k = 0; k < nz; ++k)
              Data[i][j][k] = Values.at((nx * ny) * k + nx * j + i);

      }

    void Write3DModelForX3D(const std::string &filename,
        const ThreeDModelBase::t3DModelDim &XCellSizes,
        const ThreeDModelBase::t3DModelDim &YCellSizes,
        const ThreeDModelBase::t3DModelDim &ZCellSizes,
        const double ObservationDepth,
        const ThreeDModelBase::t3DModelData &Data,
        const std::vector<double> &bg_conductivities,
        const std::vector<double> &bg_thicknesses)
      {
        assert(bg_conductivities.size() == bg_thicknesses.size());
        std::ofstream outfile(filename.c_str());
        outfile << "Version_of_X3D code (yyyy-mm-dd)\n";
        outfile << "2006-06-06\n\n";
        outfile << "Dx (m)       Dy (m)\n";
        outfile << XCellSizes[0] << "   " << YCellSizes[0] << "\n\n";
        outfile << "Background \n Thickness(m) \n";
        const double sigma_imag = 0.0;
        const double rel_eps = 1.0;
        const double rel_mu = 1.0;
        //we need a top layer with 0 thickness and 0 conductivity
        outfile << std::setw(10) << std::setprecision(5) << 0.0
            << std::setw(10) << std::setprecision(5) << 0.0 << std::setw(10)
            << std::setprecision(5) << 0.0 << std::setw(10)
            << std::setprecision(5) << 1.0 << std::setw(10)
            << std::setprecision(5) << 1.0 << "\n";
        for (size_t i = 0; i < bg_thicknesses.size(); ++i)
          {
            outfile << std::setw(10) << std::setprecision(5)
                << bg_thicknesses[i] << std::setw(15) << std::setprecision(5)
                << bg_conductivities[i] << std::setw(15)
                << std::setprecision(5) << sigma_imag << std::setw(15)
                << std::setprecision(5) << rel_eps << std::setw(15)
                << std::setprecision(5) << rel_mu << "\n";
          }
        outfile << "$\n imaginary_part_of_the_anomalous_conductivity \n N \n\n";
        outfile << " anomalous_dielectric_permittivity\n N \n\n";
        outfile
            << "First and last cells in X-direction  (nxAf_Common, nxAl_Common)\n";
        outfile << " 1 " << XCellSizes.size() << "\n\n";
        outfile
            << "First and last cells in Y-direction  (nyAf_Common, nyAl_Common) \n";
        outfile << " 1 " << YCellSizes.size() << "\n\n";

        double currdepth = 0.0;
        const size_t nzlayers = ZCellSizes.size();
        const size_t valuesperline = std::min(static_cast<size_t> (35),
            static_cast<size_t> (XCellSizes.size() + 1));
        const size_t valuewidth = 27;
        const size_t valueprec = 18;
        for (size_t i = 0; i < nzlayers; ++i)
          {
            outfile << "zA(m)  ( 'Depth_to_the_anomaly_layer' ) \n "
                << std::resetiosflags(std::ios::scientific)
                << std::setprecision(5) << currdepth << " \n\n";
            outfile << "dzA(m) ( 'Thickness_of_the_anomaly_layer' ) \n";
            outfile << std::resetiosflags(std::ios::scientific)
                << std::setprecision(5) << ZCellSizes[i] << "\n\n";
            outfile << "Thicknesses of sublayers (m) \n";
            outfile << std::resetiosflags(std::ios::scientific)
                << std::setprecision(5) << ZCellSizes[i] << "\n\n";
            currdepth += ZCellSizes[i];
            outfile
                << "Scale  ( the ARRAY will be multiplied by 'this_Scale' ) \n 1.0 \n\n";
            outfile << "First and last cells_in_X-direction \n";
            outfile << " 1 " << XCellSizes.size() << "\n\n";
            outfile << "First and last cells_in_Y-direction \n";
            outfile << " 1 " << YCellSizes.size() << "\n\n";

            outfile << "'FORMAT'\n (" << valuesperline << "E" << valuewidth
                << "." << valueprec << ") \n";
            outfile << " ARRAY ( 'scaled_and_compressed' )     \n";
            for (size_t j = 0; j < YCellSizes.size(); ++j)
              {
                for (size_t k = 0; k < XCellSizes.size(); ++k)
                  {
                    outfile << std::scientific << std::setw(valuewidth)
                        << std::setprecision(valueprec) << Data[k][j][i];
                    if (k > 0 && ((k + 1) % valuesperline) == 0)
                      outfile << "\n";
                  }
                outfile << "\n";
              }
            outfile << "\n\n";
          }//end of loop through all layers
        outfile << " First and last cells in X-direction (nxOf, nxOl)   \n";
        outfile << " 1 " << XCellSizes.size() << "\n\n";
        outfile << "   First and last cells in Y-direction (nyOf, nyOl)    \n";
        outfile << " 1 " << YCellSizes.size() << "\n\n";
        //at the moment we always set the z component of the site location to zero
        //this should be changed
        outfile << "   zO(m)  \n" << ObservationDepth << "\n\n";
        outfile
            << "Binding_cell_in_X-direction     X-coordinate of centre of Binding cell (m)  \n";
        outfile << " 1                              0.0\n";
        outfile
            << "Binding_cell_in_Y-direction     Y-coordinate of centre of Binding cell (m) \n";
        outfile << " 1                              0.0\n";

        if (outfile.bad())
          {
            throw jiba::FatalException("Problem writing model file.");
          }
      }

    void WriteProjectFile(const boost::filesystem::path &RootDir,
        const std::vector<double> &Frequencies, X3DModel::ProblemType Type,
        const std::string &ResultFilename, const std::string &ModelFilename)
      {
        const size_t nfreq = Frequencies.size();
        //the filename is always a.project
        std::ofstream outfile((RootDir / "a.project").file_string().c_str());
        //write out some information that x3d expects
        outfile << "  Version_of_X3D code (yyyy-mm-dd)\n";
        outfile << "  2006-06-06\n\n";
        outfile << "  Type_of_problem (0 - MT, 1 - CSMT, 2 - EDIP, 3 - MDIP)\n";
        outfile << " " << Type << "\n";
        //we can specify a model and a result for each frequency
        //we use the same model, but create a unique name for the frequencies
        //by appending the index at the end
        outfile
            << "  Frequency (Hz)    File_with_results     File with 3D formation      File with 1st source      File with 2nd source\n";
        for (size_t i = 0; i < nfreq; ++i)
          {
            outfile << "  " << Frequencies.at(i) << "  " << ResultFilename
                << jiba::stringify(i) << " " << ModelFilename
                << "                          " << ModelFilename
                << jiba::stringify(i) << "a.source             "
                << ModelFilename << jiba::stringify(i) << ".bsource\n";
          }
        //in principle some other parameters for the forward calculation can be changed
        //at the moment we leave them at their default values
        outfile << "$\n";
        outfile
            << "  Threshold_of relative residual norm  (default value = 0.003)\n";
        outfile << "  0.003\n\n";
        outfile
            << "  Maximum_number (Nmax) of iterates (default value = 500)\n";
        outfile << "  1000\n";
        outfile << "  Name_for_QAA (default value = value of Generic_name)\n";
        outfile << "  csmtQaa\n\n";
        outfile
            << "  Far_Zone (Y or N; default value = N; if Y it defines slow,but accurate calculus of EM fields in far zone from a dipole source)\n";
        outfile << "  Y\n\n";
        outfile
            << "  Format_of_Output (0 - mfo, 1 - ASCII, 2 - mfo+ASCII, 3 - mfo+ASCII+mfa; default value = 0)\n";
        outfile << "  3\n";
        outfile << "What_Green_function_at_step_1\n";
        outfile << "HST\n";
        outfile << "What_Green_function_at_step_4\n";
        outfile << "HST\n";

        if (outfile.bad())
          {
            throw jiba::FatalException("Problem writing project file.");
          }
      }

    void ReadEMO(const std::string &filename,
        std::vector<std::complex<double> > &Ex, std::vector<
            std::complex<double> > &Ey, std::vector<std::complex<double> > &Hx,
        std::vector<std::complex<double> > &Hy)
      {
        std::ifstream infile(filename.c_str());
        //find the description line for the electric fields
        FindToken(infile, "#        x (m)");
        //we have a few values in the file that we do not care about right now
        double dummy, real, imag;
        const std::complex<double> I(0.0, 1.0);
        char restline[2048];
        //read in numbers as long as we can, this will be the E-field
        while (infile.good())
          {
            infile >> dummy >> dummy >> dummy;
            if (infile.good())
              {
                infile >> real >> imag;
                Ex.push_back(std::complex<double>(real, -imag));
                infile >> real >> imag;
                Ey.push_back(std::complex<double>(real, -imag));
                infile.getline(restline, 2048);
              }
          }
        //swallow up the header line for the magnetic fields
        infile.clear();
        infile.getline(restline, 2048);

        //read in the magnetic fields
        while (infile.good())
          {
            infile >> dummy >> dummy >> dummy;
            if (infile.good())
              {
                infile >> real >> imag;
                Hx.push_back(std::complex<double>(real, -imag));
                infile >> real >> imag;
                Hy.push_back(std::complex<double>(real, -imag));
                infile.getline(restline, 2048);
              }
          }
        //make sure all fields have the same size
        assert(Ex.size()==Ey.size());
        assert(Ex.size()==Hx.size());
        assert(Ex.size()==Hy.size());
      }

    void ReadEMA(const std::string &filename,
        std::vector<std::complex<double> > &Ex, std::vector<
            std::complex<double> > &Ey, std::vector<std::complex<double> > &Ez,
        const size_t ncellsx, const size_t ncellsy, const size_t ncellsz)
      {
        std::ifstream infile(filename.c_str());
        //find the description line for the electric fields
        FindToken(infile, "#        x (m)");
        //we have a few values in the file that we do not care about right now
        double dummy, real, imaginary;
        //read in numbers as long as we can, this will be the E-field
        while (infile.good())
          {
            infile >> dummy >> dummy >> dummy;
            if (infile.good())
              {
                infile >> real >> imaginary;
                Ex.push_back(std::complex<double>(real, -imaginary));
                infile >> real >> imaginary;
                Ey.push_back(std::complex<double>(real, -imaginary));
                infile >> real >> imaginary;
                Ez.push_back(std::complex<double>(real, -imaginary));
              }
          }
        //make sure all fields have the same size
        assert(Ex.size() == ncellsx*ncellsy*ncellsz);
        assert(Ex.size()==Ey.size());
        assert(Ex.size()==Ez.size());
        //we use a different storage ordering then in the .ema files
        Ex = ResortFields(Ex, ncellsx, ncellsy, ncellsz);
        Ey = ResortFields(Ey, ncellsx, ncellsy, ncellsz);
        Ez = ResortFields(Ez, ncellsx, ncellsy, ncellsz);
      }

    void WriteSourceComp(std::ofstream &outfile, const boost::multi_array<
        std::complex<double>, 2> &Moments, const boost::function<
        const double &(const std::complex<double> &)> &CompFunc)
      {
        const size_t nx = Moments.shape()[0];
        const size_t ny = Moments.shape()[1];
        const size_t valuesperline = std::min(static_cast<size_t> (35), nx + 1);
        const size_t valuewidth = 27;
        const size_t valueprec = 18;

        outfile << "FORMAT\n (" << valuesperline << "E" << valuewidth << "."
            << valueprec << ") \n";
        outfile << "ARRAY\n";
        for (size_t j = 0; j < ny; ++j)
          {
            for (size_t k = 0; k < nx; ++k)
              {
                outfile << std::scientific << std::setw(valuewidth)
                    << std::setprecision(valueprec) << CompFunc(Moments[k][j]);
                if (k > 0 && ((k + 1) % valuesperline) == 0)
                  outfile << "\n";
              }
            outfile << "\n";
          }
        outfile << "$\n";
        //finished writing out all moments
      }

    void WriteGeometryInfo(std::ofstream &outfile, const size_t endx,
        const size_t endy)
      {
        outfile
            << " Scale  ( the ARRAY will be multiplied by this Scale ) \n 1.0 \n\n";
        outfile << "First and last cells in X-direction \n";
        outfile << " 1  " << endx << "\n";
        outfile << "First and last cells in Y-direction \n";
        outfile << " 1  " << endy << "\n\n";
      }

    void WriteEmptyArray(std::ofstream &outfile, const size_t XSize,
        const size_t YSize)
      {
        outfile << "ARRAY\n";
        outfile << YSize << "Lines: " << XSize << "*0.\n";
      }

    void WriteSourceFile(const std::string &filename, const double SourceDepth,
        const boost::multi_array<std::complex<double>, 2> &XPolMoments,
        const boost::multi_array<std::complex<double>, 2> &YPolMoments)
      {
        const size_t nx = XPolMoments.shape()[0];
        const size_t ny = XPolMoments.shape()[1];
        assert(nx == YPolMoments.shape()[0]);
        assert(ny == YPolMoments.shape()[1]);
        std::ofstream outfile(filename.c_str());
        outfile << "  Version_of_X3D code (yyyy-mm-dd)\n";
        outfile << "  2006-06-06\n\n";
        outfile
            << " is imaginary_part_of_the_sources assigned in this file?(Y / N; Default N)\n";
        outfile << " Y\n";
        outfile << "First and last cells in X-direction(nxSf, nxSl)\n";
        outfile << " 1  " << XPolMoments.shape()[0] << "\n";
        outfile << "First and last cells in Y-direction(nySf, nySl)\n";
        outfile << " 1  " << XPolMoments.shape()[1] << "\n";
        outfile << "zS(m)  (Depth to the source level)\n";
        outfile << SourceDepth << "\n $\n";

        //write real part of x-component
        WriteGeometryInfo(outfile, nx, ny);
        WriteSourceComp(outfile, XPolMoments, boost::bind<const double&>(
            &std::complex<double>::real, _1));

        //write imaginary part of x-component
        WriteGeometryInfo(outfile, nx, ny);
        WriteSourceComp(outfile, XPolMoments, boost::bind<const double&>(
            &std::complex<double>::imag, _1));

        //write real part of y-component
        WriteGeometryInfo(outfile, nx, ny);
        WriteSourceComp(outfile, YPolMoments, boost::bind<const double&>(
            &std::complex<double>::real, _1));

        //write imaginary part of y-component
        WriteGeometryInfo(outfile, nx, ny);
        WriteSourceComp(outfile, YPolMoments, boost::bind<const double&>(
            &std::complex<double>::imag, _1));

        //write real part of z-component, we assume it is 0
        WriteGeometryInfo(outfile, nx, ny);
        WriteEmptyArray(outfile, nx, ny);

        //write imaginary part of z-component, we assume it is 0
        WriteGeometryInfo(outfile, nx, ny);
        WriteEmptyArray(outfile, nx, ny);
        if (outfile.bad())
          {
            throw jiba::FatalException("Problem writing source file.");
          }
      }

    std::vector<std::complex<double> > ResortFields(const std::vector<
        std::complex<double> > &InField, const size_t nx, const size_t ny,
        const size_t nz)
      {
        const size_t nelements = nx * ny * nz;
        assert(nelements == InField.size());
        std::vector<std::complex<double> > result;
        result.reserve(nelements);
        for (size_t i = 0; i < nx; ++i)
          for (size_t j = 0; j < ny; ++j)
            for (size_t k = 0; k < nz; ++k)
              {
                result.push_back(InField[j + ny * i + (ny * nx) * k]);
              }
        return result;
      }

  }//end of namespace jiba
