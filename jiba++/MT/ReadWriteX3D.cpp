//============================================================================
// Name        : ReadWriteX3D.cpp
// Author      : Feb 17, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#include "../Global/FileUtil.h"
#include "../Global/FatalException.h"
#include "../Global/convert.h"
#include "ReadWriteX3D.h"
#include <fstream>
#include <iomanip>

namespace jiba
  {

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
        for (size_t i = 0; i < bg_thicknesses.size(); ++i)
          {
            outfile << std::setw(15) << std::setprecision(5)
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
        const size_t valuesperline = 90;
        const size_t valuewidth = 11;
        const size_t valueprec = 3;
        for (size_t i = 0; i < nzlayers; ++i)
          {
            outfile << "zA(m)  ( 'Depth_to_the_anomaly_layer' ) \n "
                << currdepth << " \n\n";
            outfile << "dzA(m) ( 'Thickness_of_the_anomaly_layer' ) \n";
            outfile << ZCellSizes[i] << "\n\n";
            outfile << "Thicknesses of sublayers (m) \n";
            outfile << ZCellSizes[i] << "\n\n";
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
        outfile << "   zO(m)  \n 0.0 \n\n";
        outfile
            << "Binding_cell_in_X-direction     X-coordinate of centre of Binding cell (m)  \n";
        outfile << " 1                              0.0\n";
        outfile
            << "Binding_cell_in_Y-direction     Y-coordinate of centre of Binding cell (m) \n";
        outfile << " 1                              0.0\n";
      }

    void WriteProjectFile(const std::vector<double> &Frequencies,
        X3DModel::ProblemType Type, const std::string &ResultFilename,
        const std::string &ModelFilename)
      {
        const size_t nfreq = Frequencies.size();
        std::ofstream outfile("a.project");
        outfile << "  Version_of_X3D code (yyyy-mm-dd)\n";
        outfile << "  2006-06-06\n\n";
        outfile << "  Type_of_problem (0 - MT, 1 - CSMT, 2 - EDIP, 3 - MDIP)\n";
        outfile << " " << Type << "\n";
        outfile
            << "  Frequency (Hz)    File_with_results     File with 3D formation      File with 1st source      File with 2nd source\n";
        for (size_t i = 0; i < nfreq; ++i)
          {
            outfile << "  " << Frequencies.at(i) << "  " << ResultFilename
                << jiba::stringify(i) << " " << ModelFilename
                << "                             csmt01_1.source                csmt01_2.source\n";
          }
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
        double dummy;
        const std::complex<double> I(0.0, 1.0);
        char restline[2048];
        //read in numbers as long as we can, this will be the E-field
        while (infile.good())
          {
            infile >> dummy >> dummy >> dummy;
            if (infile.good())
              {
                infile >> dummy;
                Ex.push_back(dummy);
                infile >> dummy;
                Ex.back() -= dummy * I;
                infile >> dummy;
                Ey.push_back(dummy);
                infile >> dummy;
                Ey.back() -= dummy * I;
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
                infile >> dummy;
                Hx.push_back(dummy);
                infile >> dummy;
                Hx.back() -= dummy * I;
                infile >> dummy;
                Hy.push_back(dummy);
                infile >> dummy;
                Hy.back() -= dummy * I;
                infile.getline(restline, 2048);
              }
          }
        //read in both fields
      }

    void ReadEMA(const std::string &filename,
        std::vector<std::complex<double> > &Ex, std::vector<
            std::complex<double> > &Ey)
      {
        std::ifstream infile(filename.c_str());
        //find the description line for the electric fields
        FindToken(infile, "#        x (m)");
        //we have a few values in the file that we do not care about right now
        double dummy;
        const std::complex<double> I(0.0, 1.0);
        char restline[2048];
        //read in numbers as long as we can, this will be the E-field
        while (infile.good())
          {
            infile >> dummy >> dummy >> dummy;
            if (infile.good())
              {
                infile >> dummy;
                Ex.push_back(dummy);
                infile >> dummy;
                Ex.back() += dummy * I;
                infile >> dummy;
                Ey.push_back(dummy);
                infile >> dummy;
                Ey.back() += dummy * I;
                infile.getline(restline, 2048);
              }
          }
      }

  }//end of namespace jiba

