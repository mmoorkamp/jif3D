//============================================================================
// Name        : WaveletCompressedGravityCalculator.cpp
// Author      : Dec 12, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#include "WaveletCompressedGravityCalculator.h"
#include "../Global/Wavelet.h"
#include "../Global/NumUtil.h"
#include "../Global/convert.h"
#include "../ModelBase/VTKTools.h"
#include "../ModelBase/NetCDFTools.h"

namespace jiba
  {

    WaveletCompressedGravityCalculator::WaveletCompressedGravityCalculator(
        boost::shared_ptr<ThreeDGravityImplementation> TheImp) :
      CachedGravityCalculator(TheImp), nmeas(), ngrid(), nbglayers(), nmod(),
          xsize(), ysize(), zsize(), accuracy(0.01), CurrRow(),
          transformsize(), SparseSens(), WhiteningVector()
      {

      }

    WaveletCompressedGravityCalculator::~WaveletCompressedGravityCalculator()
      {

      }

    rvec WaveletCompressedGravityCalculator::CalculateNewModel(
        const ThreeDGravityModel &Model)
      {
        nmeas = Model.GetMeasPosX().size();
        ngrid = Model.GetDensities().num_elements();
        nbglayers = Model.GetBackgroundDensities().size();
        nmod = ngrid + nbglayers;

        xsize = Model.GetDensities().shape()[0];
        ysize = Model.GetDensities().shape()[1];
        zsize = Model.GetDensities().shape()[2];

        boost::multi_array_types::size_type denssize[3];

        denssize[0] = xsize; // + 1;
        denssize[1] = ysize;
        denssize[2] = zsize;
        std::fill_n(transformsize, 3, 1.0);
        //we go through the 3 spatial directions
        //and determine the smallest power of 2 greates than the number of cells
        //we need this as we can only do a wavelet transform on data
        //with a power of 2 length
        for (size_t i = 0; i < 3; ++i)
          {
            while (transformsize[i] < denssize[i])
              {
                transformsize[i] *= 2;
              }
          }
        //allocate memory for the sensitivities for the extended domain
        CurrRow.resize(
            boost::extents[transformsize[0]][transformsize[1]][transformsize[2]]);
        SparseSens.resize(nmeas * Imp.get()->GetDataPerMeasurement(),
            CurrRow.num_elements(), false);
        SparseSens.clear();

        return Imp->Calculate(Model, *this);
      }

    rvec WaveletCompressedGravityCalculator::CalculateCachedResult(
        const ThreeDGravityModel &Model)
      {
        //if the WhiteningVector does not have the right length
        //set all elements to 1 so that it does not do anything
        //in the calculations below
        if (WhiteningVector.size() != ngrid)
          {
            WhiteningVector.resize(ngrid);
            std::fill(WhiteningVector.begin(), WhiteningVector.end(), 1.0);
          }
        std::fill_n(CurrRow.origin(), CurrRow.num_elements(), 0.0);
        //CurrRow and the Model have different sizes, but we have to
        //consider the spatial structure of the model, so we have
        //to use loops for copying and NOT std::copy
        for (size_t j = 0; j < xsize; ++j)
          for (size_t k = 0; k < ysize; ++k)
            {
              for (size_t l = 0; l < zsize; ++l)
                CurrRow[j][k][l] = Model.GetDensities()[j][k][l]
                    / WhiteningVector(j * (ysize * zsize) + k * zsize + l);
            }
        //we also need to copy the background densities to consider
        //them in the compressed sensitivity matrix
        for (size_t l = 0; l < nbglayers; ++l)
          CurrRow[xsize + 1][0][l] = Model.GetBackgroundDensities().at(l);
        jiba::WaveletTransform(CurrRow);

        jiba::rvec TransDens(CurrRow.num_elements());
        std::copy(CurrRow.origin(), CurrRow.origin() + CurrRow.num_elements(),
            TransDens.begin());

        jiba::rvec SparseResult(boost::numeric::ublas::prec_prod(SparseSens,
            TransDens));

        return SparseResult;

      }

    rvec WaveletCompressedGravityCalculator::CachedLQDerivative(
        const ThreeDGravityModel &Model, const rvec &Misfit)
      {
        const size_t ngrid = Model.GetDensities().num_elements();
        const size_t nmod = ngrid + Model.GetBackgroundThicknesses().size();
        rvec result(nmod);
        //when calculating the gradient from cached results we only
        //have to multiply the misfit vector with the transposed sensitivities
        return boost::numeric::ublas::prec_prod(trans(SparseSens), Misfit);
      }

    void WaveletCompressedGravityCalculator::HandleSensitivities(
        const size_t measindex)
      {
        if (WhiteningVector.size() != ngrid)
          {
            WhiteningVector.resize(ngrid);
            std::fill(WhiteningVector.begin(), WhiteningVector.end(), 1.0);
          }

        const size_t sparserowoffset = measindex * Imp->GetDataPerMeasurement();
        for (size_t i = 0; i < Imp->GetDataPerMeasurement(); ++i)
          {
            boost::numeric::ublas::matrix_row<jiba::rmat> SensRow(
                SetCurrentSensitivities(), i);

            //SensRow and CurrRow have different sizes, but we have to
            //consider the spatial structure of the underlying model, so we have
            //to use loops for copying and NOT std::copy
            std::fill_n(CurrRow.origin(), CurrRow.num_elements(), 0.0);
            for (size_t j = 0; j < xsize; ++j)
              for (size_t k = 0; k < ysize; ++k)
                {
                  for (size_t l = 0; l < zsize; ++l)
                    CurrRow[j][k][l] = SensRow(j * (ysize * zsize) + k * zsize
                        + l) * WhiteningVector(j * (ysize * zsize) + k * zsize
                        + l);
                }

            for (size_t l = 0; l < nbglayers; ++l)
              CurrRow[xsize + 1][0][l] = SetCurrentSensitivities()(i, xsize
                  * ysize * zsize + l);
            jiba::WaveletTransform(CurrRow);

            double normall = std::inner_product(CurrRow.origin(),
                CurrRow.origin() + CurrRow.num_elements(), CurrRow.origin(),
                0.0);
            double normdiscarded = normall;
            const double maximum = fabs(*std::max_element(CurrRow.origin(),
                CurrRow.origin() + CurrRow.num_elements(), absLess<double,
                    double> ()));
            //this is our initial estimate of the threshold
            double absthresh = 0.5 * maximum; //maximum * 2.0 * accuracy;

            while (normdiscarded / normall > accuracy)
              {
                absthresh /= 2.0;
                normdiscarded = 0.0;
                for (size_t j = 0; j < CurrRow.num_elements(); ++j)
                  {
                    if (fabs(*(CurrRow.origin() + j)) <= absthresh)

                      {
                        normdiscarded += pow(*(CurrRow.origin() + j), 2);
                      }
                  }
              }
            //now that we have determined the threshold we can assign the
            //elements to the sparse matrix
            for (size_t j = 0; j < CurrRow.num_elements(); ++j)
              {
                const double currelement = *(CurrRow.origin() + j);
                if (fabs(currelement) > absthresh)
                  {
                    SparseSens(sparserowoffset + i, j) = currelement;
                  }
              }
            //now we have thresholded this row of the sensitivity matrix
            //and we can do the next component (if we work with FTG data)
          }
      }
  }
