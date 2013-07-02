//============================================================================
// Name        : InterpolateField.cpp
// Author      : 21 Jun 2013
// Version     : 
// Copyright   : 2013, mm489
//============================================================================

#include "../Global/FatalException.h"
#include "InterpolateField.h"

namespace jif3D
  {
    //! Interpolate magnetic and electric fields horizontally to allow for arbitrary horizontal site positions
    std::complex<double> InterpolateField(std::vector<std::complex<double> > &Field,
        const X3DModel &Model, size_t MeasIndex,
        const std::vector<size_t> &MeasDepthIndices)
      {
        const size_t nmodx = Model.GetXCoordinates().size();
        const size_t nmody = Model.GetYCoordinates().size();
        //get the coordinates of the measurement site of interest
        double MeasPosX = Model.GetMeasPosX()[MeasIndex];
        double MeasPosY = Model.GetMeasPosY()[MeasIndex];
        //where are we in the model grid, we need the indices of the cell
        //in which the site is located
        boost::array<ThreeDModelBase::t3DModelData::index, 3> StationIndex =
            Model.FindAssociatedIndices(MeasPosX, MeasPosY,
                Model.GetMeasPosZ()[MeasIndex]);
        //get the position of the cell center
        double CellCenterX = Model.GetXCoordinates()[StationIndex[0]]
            + Model.GetXCellSizes()[StationIndex[0]] / 2.0;
        double CellCenterY = Model.GetYCoordinates()[StationIndex[1]]
            + Model.GetYCellSizes()[StationIndex[1]] / 2.0;
        //depending where we are with respect to the center
        //we need to use different neighbor cells
        int NextX = 0, NextY = 0;
        MeasPosX > CellCenterX ?
            NextX = StationIndex[0] + 1 : NextX = StationIndex[0] - 1;
        MeasPosY > CellCenterY ?
            NextY = StationIndex[1] + 1 : NextY = StationIndex[1] - 1;
        //if one of the indices is out of range we cannot interpolate
        //this means that the sites have to be at least gridspacing/2
        //away from the boundaries
        if (NextX < 0 || NextY < 0 || NextX >= nmodx || NextY >= nmody)
          throw FatalException("Station outside interpolation range. ");
        //get the coordinates of the centers of the adjacent cells
        //as our grid is regular we only need 2 additional coordinates
        //we can construct the coordinates of the 4 cells from these
        //and the original coordinates
        double NextCellCenterX = Model.GetXCoordinates()[NextX]
            + Model.GetXCellSizes()[NextX] / 2.0;
        double NextCellCenterY = Model.GetYCoordinates()[NextY]
            + Model.GetYCellSizes()[NextY] / 2.0;

        //calculate the offset in memory for all the field values
        const size_t Offset11 = (nmodx * nmody) * MeasDepthIndices[MeasIndex]
            + StationIndex[0] * nmody + StationIndex[1];
        const size_t Offset12 = (nmodx * nmody) * MeasDepthIndices[MeasIndex]
            + StationIndex[0] * nmody + NextY;
        const size_t Offset21 = (nmodx * nmody) * MeasDepthIndices[MeasIndex]
            + NextX * nmody + StationIndex[1];
        const size_t Offset22 = (nmodx * nmody) * MeasDepthIndices[MeasIndex]
            + NextX * nmody + NextY;
        //implement bilinear interpolation equation on a regular grid
        std::complex<double> InterField = Field[Offset11] * (NextCellCenterX - MeasPosX)
            * (NextCellCenterY - MeasPosY);
        InterField += Field[Offset21] * (MeasPosX - CellCenterX)
            * (NextCellCenterY - MeasPosY);
        InterField += Field[Offset12] * (NextCellCenterX - MeasPosX)
            * (MeasPosY - CellCenterY);
        InterField += Field[Offset22] * (MeasPosX - CellCenterX)
            * (MeasPosY - CellCenterY);
        InterField /= (NextCellCenterX - CellCenterX) * (NextCellCenterY - CellCenterY);
        return InterField;
      }
  }
