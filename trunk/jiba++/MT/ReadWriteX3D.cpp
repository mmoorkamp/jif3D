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
    static const size_t maxlength = 2048;

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
        //find the line in the file that describes the cell size in the horizontal directions
        //all cells have the same size
        std::string line = FindToken(infile, "Dx");
        double dx = 2.1, dy = 3.4;
        infile >> dx >> dy;

        //find the information about the background
        line = FindToken(infile, "Thickness");
        double currvalue;
        //swallow the line with 0 thickness and conductivity
        char dummy[maxlength];
        infile.getline(dummy, maxlength);
        //read in background layer thickness and conductivity
        //as long as we can read in two consecutive numbers
        while (!infile.fail())
          {
            infile >> currvalue;
            //we read in and if the conversion went OK
            //we assign it to the background
            if (!infile.fail())
              bg_thicknesses.push_back(currvalue);
            infile >> currvalue;
            if (!infile.fail())
              {
                char dummy[maxlength];
                bg_conductivities.push_back(currvalue);
                infile.getline(dummy, maxlength);
              }
          }
        //clear the failbit
        infile.clear();
        std::vector<double> Zthick;
        std::vector<double> Values;
        bool havelayer = true;

        int startx, starty, endx, endy;
        //we do not know how many layers we have beforehand
        //so we search for the next occurence of dzA(m)
        //if it is there we read in another layer
        //if not there are no more layers
        while (havelayer)
          {
            //FindToken throws an exception if it cannot find the token
            //so we catch that and just signal that there are no more layers
            try
              {
                line = FindToken(infile, "dzA(m)");
              } catch (FatalException &e)
              {
                havelayer = false;
              }
            //if there is another layer in the file read in the information
            if (havelayer)
              {
                //first the layer thickness
                infile >> currvalue;
                Zthick.push_back(currvalue);
                //we read in the amount of cells in x-direction and y-direction
                //but we assume that these numbers are the same for each layer
                //In the file these numbers can vary, but we do not perform
                //any checking as this case is not interesting for us
                line = FindToken(infile, "cells_in_X-direction");
                infile >> startx >> endx;
                FindToken(infile, "cells_in_Y-direction");
                infile >> starty >> endy;
                FindToken(infile, "ARRAY");
                //calculate how many cells in a layer
                const unsigned int nelements = (endx - startx + 1) * (endy
                    - starty + 1);
                //and read them in the order they were written
                for (size_t i = 0; i < nelements; ++i)
                  {
                    infile >> currvalue;
                    Values.push_back(currvalue);
                  }
              }
          }//end of while
        //now allocate appropriate memory to store the model information
        //again note that we assume that all layers have the same size
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
        //we use a different storage ordering than x3d
        //so we resort the cell values while copying
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
        //write out some header information
        outfile << "Version_of_X3D code (yyyy-mm-dd)\n";
        outfile << "2006-06-06\n\n";
        //all cells have the same siyes in each horizontal direction
        //so we write only the first value, the access function
        //guarantees that all values in the size arrays are identical
        outfile << "Dx (m)       Dy (m)\n";
        outfile << XCellSizes[0] << "   " << YCellSizes[0] << "\n\n";
        //hten we write information about the surrounding backgroud
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
        //and then we write the background layers we specified in the model
        for (size_t i = 0; i < bg_thicknesses.size(); ++i)
          {
            outfile << std::setw(10) << std::setprecision(5)
                << bg_thicknesses[i] << std::setw(15) << std::setprecision(5)
                << bg_conductivities[i] << std::setw(15)
                << std::setprecision(5) << sigma_imag << std::setw(15)
                << std::setprecision(5) << rel_eps << std::setw(15)
                << std::setprecision(5) << rel_mu << "\n";
          }
        //there are some options in x3d that we do not support
        //so we always answer no
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
        //there can only be 1024 character per line in the ascii file
        //so we format the entries for high precision and adjust the number of values per line
        const size_t valuesperline = std::min(static_cast<size_t> (35),
            static_cast<size_t> (XCellSizes.size() + 1));
        const size_t valuewidth = 27;
        const size_t valueprec = 18;
        //write out the model grid layer by layer in z-direction
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
            //we do not apply any scaling to the conductivity values, but write them out as is
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
        //write out the depth at which the sites are located
        //we assume that all sites are at the same depth
        outfile << "   zO(m)  \n" << ObservationDepth << "\n\n";
        //the first cell always has the coordinate 0
        outfile
            << "Binding_cell_in_X-direction     X-coordinate of centre of Binding cell (m)  \n";
        outfile << " 1                             " << XCellSizes[0]/2.0 << " \n";
        outfile
            << "Binding_cell_in_Y-direction     Y-coordinate of centre of Binding cell (m) \n";
        outfile << " 1                             " << YCellSizes[0]/2.0 << " \n";
        //if something went wrong during writing we throw an exception
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
        char restline[maxlength];
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
                infile.getline(restline, maxlength);
              }
          }
        //swallow up the header line for the magnetic fields
        infile.clear();
        infile.getline(restline, maxlength);

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
                infile.getline(restline, maxlength);
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
        if (Ex.size() != ncellsx * ncellsy * ncellsz)
          {
            throw jiba::FatalException(
                "In ReadEma, number of electric field values does not match grid size");
          }
        if (Ex.size() != Ey.size())
          {
            throw jiba::FatalException("In ReadEma, size of Ex != Ey");
          }
        if (Ex.size() != Ez.size())
          {
            throw jiba::FatalException("In ReadEma, size of Ex != Ez");
          }
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
        //write out the required header for the source file
        //x3d is very picky about this
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

