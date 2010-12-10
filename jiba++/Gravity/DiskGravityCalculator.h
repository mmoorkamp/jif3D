/*
 * DiskGravityCalculator.h
 *
 *  Created on: Jun 3, 2009
 *      Author: mmoorkamp
 */

#ifndef DISKGRAVITYCALCULATOR_H_
#define DISKGRAVITYCALCULATOR_H_

#include "FullSensitivityGravityCalculator.h"

namespace jiba
  {
    //! This class works similar to FullSensitivityCalculator, only that it stores the sensitivities on the disk, not in memory
    /*! When we do not have enough RAM to store the sensitivities we can store them on the disk if
     * we do not want to recalculate them each time. On construction we generate a filename that consists
     * of the ID of the program and the memory location of the object and therefore should be unique.
     */
    class DiskGravityCalculator: public jiba::FullSensitivityGravityCalculator
      {
    private:
      //! This routine creates a unique filename for each object, used upon construction
      std::string MakeFilename();
      //the filename where we store the sensitivities, generated by the constructor
      std::string filename;
      //calculate the raw data without any transformation from the sensitivities stored in the file
      virtual rvec CalculateRawData(const ThreeDGravityModel &Model);
      //when we calculate a completely new model we have to delete the file with the old sensitivities
      virtual rvec CalculateNewModel(const ThreeDGravityModel &Model);
      //calculate the raw derivative without any transformation from the sensitivities stored in the file
      virtual rvec CalculateRawLQDerivative(const ThreeDGravityModel &Model,
          const rvec &Misfit);
    public:
      //! Handle the current row of the sensitivity matrix by storing it in a file
      virtual void HandleSensitivities(const size_t measindex);
      //! The constructor needs a shared pointer to an implementation object, usually handled by CreateGravityCalculator
      explicit DiskGravityCalculator(
          boost::shared_ptr<ThreeDGravityImplementation> TheImp);
      //! We need to define the copy constructor to make sure that filename stays unique among all created objects
      explicit DiskGravityCalculator(const DiskGravityCalculator &Old);
      //! We have to define a copy operator to make sure filename stays unique
      DiskGravityCalculator& operator=(const DiskGravityCalculator& source);
      virtual ~DiskGravityCalculator();
      };

  }

#endif /* DISKGRAVITYCALCULATOR_H_ */
