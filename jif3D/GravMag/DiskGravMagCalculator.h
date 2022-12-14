/*
 * DiskGravityCalculator.h
 *
 *  Created on: Jun 3, 2009
 *      Author: mmoorkamp
 */

#ifndef DISKGRAVITYCALCULATOR_H_
#define DISKGRAVITYCALCULATOR_H_

#include "../Global/Serialization.h"
#include "../Global/Jif3DGlobal.h"
#include "../Global/convert.h"
#include "../Global/VecMat.h"
#include "../Global/FatalException.h"
#include "../Global/Jif3DPlatformHelper.h"
#include <fstream>
#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/filesystem.hpp>
#include "FullSensitivityGravMagCalculator.h"

namespace jif3D
  {
    //! This class works similar to FullSensitivityCalculator, only that it stores the sensitivities on the disk, not in memory
    /*! When we do not have enough RAM to store the sensitivities we can store them on the disk if
     * we do not want to recalculate them each time. On construction we generate a filename that consists
     * of the ID of the program and the memory location of the object and therefore should be unique.
     */
    template<class PotentialDataType>
    class J3DEXPORT DiskGravMagCalculator: public jif3D::FullSensitivityGravMagCalculator<
        PotentialDataType>
      {
    private:
      typedef typename PotentialDataType::ModelType ThreeDModelType;

      //! This routine creates a unique filename for each object, used upon construction
      std::string MakeFilename();
      //! The directory and name of the file to write the file with the sensitivities
      boost::filesystem::path FullPath;
      //calculate the raw data without any transformation from the sensitivities stored in the file
      virtual rvec CalculateRawData(const ThreeDModelType &Model,
          const PotentialDataType &Data) override;
      //when we calculate a completely new model we have to delete the file with the old sensitivities
      virtual rvec CalculateNewModel(const ThreeDModelType &Model,
          const PotentialDataType &Data) override;
      //calculate the raw derivative without any transformation from the sensitivities stored in the file
      virtual rvec CalculateRawLQDerivative(const ThreeDModelType &Model,
          const PotentialDataType &Data, const rvec &Misfit) override;
      friend class access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & base_object<FullSensitivityGravMagCalculator<PotentialDataType> >(*this);
          ar & FullPath;
        }
    public:
      //! Handle the current row of the sensitivity matrix by storing it in a file
      virtual void HandleSensitivities(const size_t measindex) override;
      //! The constructor needs a shared pointer to an implementation object, usually handled by CreateGravityCalculator
      explicit DiskGravMagCalculator(
          boost::shared_ptr<ThreeDGravMagImplementation<PotentialDataType> > TheImp,
          boost::filesystem::path TDir = boost::filesystem::current_path());
      //! We need to define the copy constructor to make sure that filename stays unique among all created objects
      explicit DiskGravMagCalculator(const DiskGravMagCalculator &Old);
      //! We have to define a copy operator to make sure filename stays unique
      DiskGravMagCalculator& operator=(const DiskGravMagCalculator& source);
      virtual ~DiskGravMagCalculator();
      };

    template<class PotentialDataType>
    std::string DiskGravMagCalculator<PotentialDataType>::MakeFilename()
      {
        //a unique ID created on construction
        boost::uuids::uuid tag = boost::uuids::random_generator()();
        //make a unique filename for the sensitivity file created by this object
        //we use boost uuid to generate a unique identifier tag
        //and translate it to a string to generate the filename
        return "grav" + jif3D::stringify(jif3D::platform::get_process_id())
            + jif3D::stringify(tag);
      }

    template<class PotentialDataType>
    rvec DiskGravMagCalculator<PotentialDataType>::CalculateNewModel(
        const ThreeDModelType &Model, const PotentialDataType &Data)
      {
        //when we have to calculate a new model
        //delete the file with the old sensitivities
        boost::filesystem::remove_all(FullPath);
        //then forward the call to the implementation object
        return ThreeDGravMagCalculator<PotentialDataType>::Imp.get()->Calculate(Model,
            Data, *this);
      }

    template<class PotentialDataType>
    void DiskGravMagCalculator<PotentialDataType>::HandleSensitivities(
        const size_t measindex)
      {
        //whenever we have a new row of the sensitivity matrix
        //we append to the existing file
        std::fstream outfile(FullPath.string().c_str(),
            std::ios::out | std::ios::binary | std::ios::app);
        //depending on whether we have FTG or scalar data
        //the current segment of the sensitivity matrix can have several rows
        const size_t nrows =
            ThreeDGravMagCalculator<PotentialDataType>::SetCurrentSensitivities().size1();
        const size_t ncolumns =
            ThreeDGravMagCalculator<PotentialDataType>::SetCurrentSensitivities().size2();
        for (size_t i = 0; i < nrows; ++i)
          {
            //we have to copy the current row in a vector
            //because the matrix is stored in column major order
            jif3D::rvec Row(
                ublas::matrix_row<jif3D::rmat>(
                    ThreeDGravMagCalculator<PotentialDataType>::SetCurrentSensitivities(),
                    i));
            //write the current row to a file
            outfile.write(reinterpret_cast<char *>(&Row[0]), ncolumns * sizeof(double));
          }
      }

    template<class PotentialDataType>
    rvec DiskGravMagCalculator<PotentialDataType>::CalculateRawData(
        const ThreeDModelType &Model, const PotentialDataType &Data)
      {
        //open the file where the sensitivities are stored
        std::fstream infile(FullPath.string().c_str(), std::ios::in | std::ios::binary);

        const size_t nmeas =
            Data.GetMeasPosX().size()
                * ThreeDGravMagCalculator<PotentialDataType>::Imp.get()->RawDataPerMeasurement();
        const size_t nmod = Model.GetNModelParm();

        rvec Vector = Model.GetModelParameters();

        rvec result(nmeas);
        rvec CurrSens(nmod);
        //for each measurement read the current row from the binary file
        //and perform a scalar product between the sensitivities and the densities
        for (size_t i = 0; i < nmeas; ++i)
          {
            infile.read(reinterpret_cast<char *>(&CurrSens.data()[0]),
                nmod * sizeof(double));
            if (infile.good())
              {
                result(i) = boost::numeric::ublas::inner_prod(CurrSens, Vector);
              }
            else
              {
                std::string error =
                    "In CalculateRawData, cannot read sensitivities from binary file: "
                        + FullPath.string();
                throw FatalException(error, __FILE__, __LINE__);
              }
          }
        return result;

      }
    template<class PotentialDataType>
    rvec DiskGravMagCalculator<PotentialDataType>::CalculateRawLQDerivative(
        const ThreeDModelType &Model, const PotentialDataType &Data, const rvec &Misfit)
      {
        //when we are in this routine we read the sensitivities
        //from a previously created binary file
        std::fstream infile(FullPath.string().c_str(), std::ios::in | std::ios::binary);

        const size_t nmeas =
            Data.GetMeasPosX().size()
                * ThreeDGravMagCalculator<PotentialDataType>::Imp.get()->RawDataPerMeasurement();
        const size_t nmod = Model.GetNModelParm();
        //we read the sensitivities row by row and multiply
        //by the corresponding misfit to calculate the gradient
        rvec Gradient(nmod, 0.0);
        rvec CurrSens(nmod, 0.0);
        //for each measurement
        for (size_t i = 0; i < nmeas; ++i)
          {
            //read in one row of the sensitivity matrix from the binary file
            infile.read(reinterpret_cast<char *>(&CurrSens.data()[0]),
                nmod * sizeof(double));
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
                throw FatalException(error, __FILE__, __LINE__);
              }
          }
        return Gradient;
      }

    template<class PotentialDataType>
    DiskGravMagCalculator<PotentialDataType>::DiskGravMagCalculator(
        boost::shared_ptr<ThreeDGravMagImplementation<PotentialDataType> > TheImp,
        boost::filesystem::path TDir) :
        FullSensitivityGravMagCalculator<PotentialDataType>(TheImp)
      {
        std::cout << "Storing sensitivities on disk " << std::endl;
        if (!boost::filesystem::is_directory(TDir))
          throw FatalException("TDir is not a directory: " + TDir.string(), __FILE__,
          __LINE__);
        FullPath = TDir / MakeFilename();
      }

    //! We need to define the copy constructor to make sure that filename stays unique among all created objects
    template<class PotentialDataType>
    DiskGravMagCalculator<PotentialDataType>::DiskGravMagCalculator(
        const DiskGravMagCalculator<PotentialDataType> &Old) :
        FullSensitivityGravMagCalculator<PotentialDataType>(Old)
      {
        //we keep the original directory, as this can be a temp directory
        //but generate a new filename for the new object to avoid overwriting
        //the file from the original object
        FullPath = Old.FullPath.parent_path() / MakeFilename();
      }

    //! We have to define a copy operator to make sure filename stays unique
    template<class PotentialDataType>
    DiskGravMagCalculator<PotentialDataType>& DiskGravMagCalculator<PotentialDataType>::operator=(
        const DiskGravMagCalculator<PotentialDataType>& source)
      {
        if (this == &source)
          return *this;
        FullSensitivityGravMagCalculator<PotentialDataType>::operator=(source);
        //we only have to copy the base class information
        //DiskGravityCalculator only adds the field filename
        //this should be unqiue for each object, so we do
        //not copy it, but keep the name generated at construction
        return *this;
      }

    template<class PotentialDataType>
    DiskGravMagCalculator<PotentialDataType>::~DiskGravMagCalculator()
      {
        //make sure we clean up the file when the object disappears
        boost::filesystem::remove_all(FullPath);
      }
  }

#endif /* DISKGRAVITYCALCULATOR_H_ */
