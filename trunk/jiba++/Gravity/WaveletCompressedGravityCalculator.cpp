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
      CachedGravityCalculator(TheImp), accuracy(0.01)
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
        SparseSens.clear();

        return Imp->Calculate(Model, *this);
      }

    rvec WaveletCompressedGravityCalculator::CalculateCachedResult(
        const ThreeDGravityModel &Model)
      {
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

        for (size_t l = 0; l < nbglayers; ++l)
          CurrRow[xsize + 1][0][l] = Model.GetBackgroundDensities().at(l);
        jiba::WaveletTransform(CurrRow);

        //debug code start
        jiba::ThreeDGravityModel::t3DModelDim FakeSizeX(
            boost::extents[CurrRow.shape()[0]]), FakeSizeY(
            boost::extents[CurrRow.shape()[1]]), FakeSizeZ(
            boost::extents[CurrRow.shape()[2]]);
        std::fill_n(FakeSizeX.origin(), CurrRow.shape()[0], 1.0);
        std::fill_n(FakeSizeY.origin(), CurrRow.shape()[1], 1.0);
        std::fill_n(FakeSizeZ.origin(), CurrRow.shape()[2], 1.0);
        jiba::Write3DModelToVTK("waveletmod.vtk", "WaveletMod", FakeSizeX,
            FakeSizeY, FakeSizeZ, CurrRow);
        //debug code end

        jiba::rvec TransDens(CurrRow.num_elements());
        std::fill_n(TransDens.begin(), CurrRow.num_elements(), 0.0);
        std::copy(CurrRow.origin(), CurrRow.origin() + CurrRow.num_elements(),
            TransDens.begin());

        jiba::rvec SparseResult(boost::numeric::ublas::prec_prod(SparseSens,
            TransDens));

        //bad model identification start
        const size_t nelements = SparseSens.size2();
        jiba::rvec BadModel(nelements);
        for (size_t j = nelements/2; j < nelements; ++j)
          {
            if (SparseSens(0, j) != 0.0)
              {
                BadModel( j) = 0.0;
              }
            else
              {
                BadModel( j) = 1.0;
              }
          }
        std::copy(BadModel.begin(),BadModel.end(),CurrRow.origin());
        jiba::InvWaveletTransform(CurrRow);
        jiba::Write3DModelToVTK("badmod.vtk", "BadMod", Model.GetXCellSizes(),
            Model.GetYCellSizes(),  Model.GetZCellSizes(), CurrRow);
        jiba::Write3DModelToNetCDF("badmod.nc", "density","g/cm3", Model.GetXCellSizes(),
            Model.GetYCellSizes(),  Model.GetZCellSizes(), CurrRow);
        //bad model identification end

        return SparseResult;

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

            //debug code start
            jiba::ThreeDGravityModel::t3DModelDim FakeSizeX(
                boost::extents[CurrRow.shape()[0]]), FakeSizeY(
                boost::extents[CurrRow.shape()[1]]), FakeSizeZ(
                boost::extents[CurrRow.shape()[2]]);
            std::fill_n(FakeSizeX.origin(), CurrRow.shape()[0], 1.0);
            std::fill_n(FakeSizeY.origin(), CurrRow.shape()[1], 1.0);
            std::fill_n(FakeSizeZ.origin(), CurrRow.shape()[2], 1.0);
            jiba::Write3DModelToVTK(
                "waveletsens" + jiba::stringify(i) + ".vtk", "WaveletSens",
                FakeSizeX, FakeSizeY, FakeSizeZ, CurrRow);
            //debug code end

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
                std::cout << "Normall: " << normall << " Normdiscarded: "
                    << normdiscarded << std::endl;
              }
            std::cout << "Abs thresh: " << absthresh << " " << maximum
                << std::endl;
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
