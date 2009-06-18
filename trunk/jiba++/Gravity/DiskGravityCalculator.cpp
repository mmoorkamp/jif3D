/*
 * DiskGravityCalculator.cpp
 *
 *  Created on: Jun 3, 2009
 *      Author: mmoorkamp
 */

#include <unistd.h>
#include <boost/filesystem.hpp>
#include "../Global/convert.h"
#include "../Global/VecMat.h"
#include "../Global/FatalException.h"
#include "DiskGravityCalculator.h"
#include <fstream>

namespace jiba
  {

    rvec DiskGravityCalculator::CalculateNewModel(
        const ThreeDGravityModel &Model)
      {
        //when we have to calculate a new model
        //delete the file with the old sensitivities
        boost::filesystem::remove_all(filename);
        //then forward the call to the implementation object
        return Imp.get()->Calculate(Model, *this);
      }

    void DiskGravityCalculator::HandleSensitivities(const size_t measindex)
      {
        //whenever we have a new row of the sensitivity matrix
        //we append to the existing file
        std::fstream outfile(filename.c_str(), std::ios::out | std::ios::binary
            | std::ios::app);
        //depending on whether we have FTG or scalar data
        //the current segment of the sensitivity matrix can have several rows
        const size_t nrows = SetCurrentSensitivities().size1();
        const size_t ncolumns = SetCurrentSensitivities().size2();
        for (size_t i = 0; i < nrows; ++i)
          {
            //we have to copy the current row in a vector
            //because the matrix is stored in column major order
            jiba::rvec Row(ublas::matrix_row<jiba::rmat>(
                SetCurrentSensitivities(), i));
            //write the current row to a file
            outfile.write(reinterpret_cast<char *> (&Row[0]), ncolumns
                * sizeof(double));
          }
      }

    rvec DiskGravityCalculator::CalculateRawData(
        const ThreeDGravityModel &Model)
      {
        //open the file where the sensitivities are stored
        std::fstream infile(filename.c_str(), std::ios::in | std::ios::binary);

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
                throw FatalException(
                    "Cannot read sensitivities from binary file !");
              }
          }
        return result;

      }

    rvec DiskGravityCalculator::CalculateRawLQDerivative(
        const ThreeDGravityModel &Model, const rvec &Misfit)
      {
        std::fstream infile(filename.c_str(), std::ios::in | std::ios::binary);

        const size_t nmeas = Model.GetMeasPosX().size()
            * Imp.get()->RawDataPerMeasurement();
        const size_t ngrid = Model.GetDensities().num_elements();
        const size_t nmod = ngrid + Model.GetBackgroundThicknesses().size();
        rvec result(nmod);
        result.clear();
        rvec CurrSens(nmod);
        for (size_t i = 0; i < nmeas; ++i)
          {
            infile.read(reinterpret_cast<char *> (&CurrSens.data()[0]), nmod
                * sizeof(double));
            result += CurrSens * Misfit(i);
          }
        return result;
      }

    DiskGravityCalculator::DiskGravityCalculator(boost::shared_ptr<
        ThreeDGravityImplementation> TheImp) :
      FullSensitivityGravityCalculator(TheImp)
      {
        //make a unique filename for the sensitivity file created by this object
        filename = "grav" + jiba::stringify(getpid()) + jiba::stringify(this);
      }

    DiskGravityCalculator::~DiskGravityCalculator()
      {
        //make sure we clean up the file when the object disappears
        boost::filesystem::remove_all(filename);
      }

  }
