//============================================================================
// Name        : WaveletCompressedGravityCalculator.cpp
// Author      : Dec 12, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#include "WaveletCompressedGravityCalculator.h"
#include "../Global/Wavelet.h"

namespace jiba
  {

    WaveletCompressedGravityCalculator::WaveletCompressedGravityCalculator(
        boost::shared_ptr<ThreeDGravityImplementation> TheImp) :
      CachedGravityCalculator(TheImp)
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

        denssize[0] = xsize + 1;
        denssize[1] = ysize;
        denssize[2] = zsize;
        std::fill_n(transformsize, 3, 1.0);
        for (size_t i = 0; i < 3; ++i)
          {
            while (transformsize[i] < denssize[i])
              {
                transformsize[i] *= 2;
              }
          }

        CurrRow.resize(
            boost::extents[transformsize[0]][transformsize[1]][transformsize[2]]);
        SparseSens.resize(nmeas * Imp.get()->GetDataPerMeasurement(),
            CurrRow.num_elements(), false);

        return Imp->Calculate(Model, *this);
      }

    rvec WaveletCompressedGravityCalculator::CalculateCachedResult(
        const ThreeDGravityModel &Model)
      {

        std::fill_n(CurrRow.origin(), CurrRow.num_elements(), 0.0);
        for (size_t j = 0; j < xsize; ++j)
          for (size_t k = 0; k < ysize; ++k)
            {
              for (size_t l = 0; l < zsize; ++l)
                CurrRow[j][k][l] = Model.GetDensities()[j][k][l];
            }

        for (size_t l = 0; l < nbglayers; ++l)
          CurrRow[xsize + 1][0][l] = Model.GetBackgroundDensities().at(l);
        jiba::WaveletTransform(CurrRow);

        jiba::rvec TransDens(CurrRow.num_elements());
        std::fill_n(TransDens.begin(), CurrRow.num_elements(), 0.0);
        std::copy(CurrRow.origin(), CurrRow.origin() + CurrRow.num_elements(),
            TransDens.begin());

        jiba::rvec SparseResult(boost::numeric::ublas::prec_prod(SparseSens,
            TransDens));
        return SparseResult;

      }

    void WaveletCompressedGravityCalculator::HandleSensitivities(
        const size_t measindex)
      {
        for (size_t i = 0; i < Imp->GetDataPerMeasurement(); ++i)
          {
            std::fill_n(CurrRow.origin(), CurrRow.num_elements(), 0.0);
            for (size_t j = 0; j < xsize; ++j)
              for (size_t k = 0; k < ysize; ++k)
                {
                  for (size_t l = 0; l < zsize; ++l)
                    CurrRow[j][k][l] = SetCurrentSensitivities()(i, j * (ysize
                        * zsize) + k * zsize + l);

                }

            for (size_t l = 0; l < nbglayers; ++l)
              CurrRow[xsize + 1][0][l] = SetCurrentSensitivities()(i, xsize
                  * ysize * zsize + l);
            jiba::WaveletTransform(CurrRow);
            size_t currn = 0;
            double relthresh = 1e-1;
            double normall = std::inner_product(CurrRow.origin(),
                CurrRow.origin() + CurrRow.num_elements(), CurrRow.origin(),
                0.0);
            double normdiscarded = normall;
            const double maximum = *std::max_element(CurrRow.origin(),
                CurrRow.origin() + CurrRow.num_elements());
            while (normdiscarded / normall > 1e-2)
              {
                currn = 0;
                normdiscarded = 0.0;
                const double absthresh = maximum * relthresh;
                for (size_t j = 0; j < CurrRow.num_elements(); ++j)
                  {
                    if (fabs(*(CurrRow.origin() + j)) > absthresh)
                      {
                        SparseSens(measindex + i * Imp.get()->GetDataPerMeasurement(),
                            j) = *(CurrRow.origin() + j);
                        currn++;
                      }
                    else
                      {
                        normdiscarded += pow(*(CurrRow.origin() + j), 2);
                      }
                  }
                relthresh /= 2.0;
              }
          }
      }
  }
