/*
 * DiskGravityCalculator.cpp
 *
 *  Created on: Jun 3, 2009
 *      Author: mmoorkamp
 */

#include <unistd.h>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/filesystem.hpp>
#include "../Global/convert.h"
#include "../Global/VecMat.h"
#include "../Global/FatalException.h"
#include "DiskGravMagCalculator.h"
#include <fstream>

namespace jif3D
  {

    std::string DiskGravMagCalculator::MakeFilename()
      {
        //a unique ID created on construction
        boost::uuids::uuid tag = boost::uuids::random_generator()();
        //make a unique filename for the sensitivity file created by this object
        //we use boost uuid to generate a unique identifier tag
        //and translate it to a string to generate the filename
        return "grav" + jif3D::stringify(getpid()) + jif3D::stringify(tag);
      }

    rvec DiskGravMagCalculator::CalculateNewModel(
        const ThreeDGravityModel &Model)
      {
        //when we have to calculate a new model
        //delete the file with the old sensitivities
        boost::filesystem::remove_all(FullPath);
        //then forward the call to the implementation object
        return Imp.get()->Calculate(Model, *this);
      }

    void DiskGravMagCalculator::HandleSensitivities(const size_t measindex)
      {
        //whenever we have a new row of the sensitivity matrix
        //we append to the existing file
        std::fstream outfile(FullPath.string().c_str(), std::ios::out | std::ios::binary
            | std::ios::app);
        //depending on whether we have FTG or scalar data
        //the current segment of the sensitivity matrix can have several rows
        const size_t nrows = SetCurrentSensitivities().size1();
        const size_t ncolumns = SetCurrentSensitivities().size2();
        for (size_t i = 0; i < nrows; ++i)
          {
            //we have to copy the current row in a vector
            //because the matrix is stored in column major order
            jif3D::rvec Row(ublas::matrix_row<jif3D::rmat>(
                SetCurrentSensitivities(), i));
            //write the current row to a file
            outfile.write(reinterpret_cast<char *> (&Row[0]), ncolumns
                * sizeof(double));
          }
      }

    rvec DiskGravMagCalculator::CalculateRawData(
        const ThreeDGravityModel &Model)
      {
        //open the file where the sensitivities are stored
        std::fstream infile(FullPath.string().c_str(), std::ios::in | std::ios::binary);

        const size_t nmeas = Model.GetMeasPosX().size()
            * Imp.get()->RawDataPerMeasurement();
        const size_t ngrid = Model.GetDensities().num_elements();
        const size_t nmod = ngrid + Model.GetBackgroundThicknesses().size();

        rvec DensVector(nmod);
        //copy the 3D model structure and the background into a vector of densities
        std::copy(Model.GetDensities().origin(), Model.GetDensities().origin()
            + ngrid, DensVector.begin());
        std::copy(Model.GetBackgroundDensities().begin(),
            Model.GetBackgroundDensities().end(), DensVector.begin() + ngrid);
        rvec result(nmeas);
        rvec CurrSens(nmod);
        //for each measurement read the current row from the binary file
        //and perform a scalar product between the sensitivities and the densities
        for (size_t i = 0; i < nmeas; ++i)
          {
            infile.read(reinterpret_cast<char *> (&CurrSens.data()[0]), nmod
                * sizeof(double));
            if (infile.good())
              {
                result(i) = boost::numeric::ublas::inner_prod(CurrSens,
                    DensVector);
              }
            else
              {
                std::string error =
                    "In CalculateRawData, cannot read sensitivities from binary file: "
                        + FullPath.string();
                throw FatalException(error);
              }
          }
        return result;

      }

    rvec DiskGravMagCalculator::CalculateRawLQDerivative(
        const ThreeDGravityModel &Model, const rvec &Misfit)
      {
        //when we are in this routine we read the sensitivities
        //from a previously created binary file
        std::fstream infile(FullPath.string().c_str(), std::ios::in | std::ios::binary);

        const size_t nmeas = Model.GetMeasPosX().size()
            * Imp.get()->RawDataPerMeasurement();
        const size_t ngrid = Model.GetDensities().num_elements();
        const size_t nmod = ngrid + Model.GetBackgroundThicknesses().size();
        //we read the sensitivities row by row and multiply
        //by the corresponding misfit to calculate the gradient
        rvec Gradient(nmod, 0.0);
        rvec CurrSens(nmod, 0.0);
        //for each measurement
        for (size_t i = 0; i < nmeas; ++i)
          {
            //read in one row of the sensitivity matrix from the binary file
            infile.read(reinterpret_cast<char *> (&CurrSens.data()[0]), nmod
                * sizeof(double));
            if (infile.good())
              {
                //if reading was successful, we can compute the gradient
                //this is equivalent to J^T * delta d
                Gradient += CurrSens * Misfit(i);
              }
            else
              {
                std::string error =
                    "In CalculateRawLQDerivative, cannot read sensitivities from binary file: "
                        + FullPath.string();
                throw FatalException(error);
              }
          }
        return Gradient;
      }

    DiskGravMagCalculator::DiskGravMagCalculator(boost::shared_ptr<
        ThreeDGravMagImplementation> TheImp, boost::filesystem::path TDir) :
      FullSensitivityGravMagCalculator(TheImp)
      {
        if (!boost::filesystem::is_directory(TDir))
          throw FatalException("TDir is not a directory: " + TDir.string());
         FullPath = TDir / MakeFilename();
      }

    //! We need to define the copy constructor to make sure that filename stays unique among all created objects
    DiskGravMagCalculator::DiskGravMagCalculator(
        const DiskGravMagCalculator &Old) :
      FullSensitivityGravMagCalculator(Old)
      {
        //we keep the original directory, as this can be a temp directory
        //but generate a new filename for the new object to avoid overwriting
        //the file from the original object
        FullPath = Old.FullPath.parent_path() / MakeFilename() ;
      }

    //! We have to define a copy operator to make sure filename stays unique
    DiskGravMagCalculator& DiskGravMagCalculator::operator=(
        const DiskGravMagCalculator& source)
      {
        if (this == &source)
          return *this;
        FullSensitivityGravMagCalculator::operator=(source);
        //we only have to copy the base class information
        //DiskGravityCalculator only adds the field filename
        //this should be unqiue for each object, so we do
        //not copy it, but keep the name generated at construction
        return *this;
      }

    DiskGravMagCalculator::~DiskGravMagCalculator()
      {
        //make sure we clean up the file when the object disappears
        boost::filesystem::remove_all(FullPath);
      }

  }
