//============================================================================
// Name        : ReadWriteX3D.cpp
// Author      : Feb 17, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#include "../Global/FileUtil.h"
#include "../Global/FatalException.h"
#include "ReadWriteX3D.h"
#include <fstream>
namespace jiba
  {

    void Read3DModelFromX3D(const std::string &filename,
        ThreeDModelBase::t3DModelDim &XCellSizes,
        ThreeDModelBase::t3DModelDim &YCellSizes,
        ThreeDModelBase::t3DModelDim &ZCellSizes,
        ThreeDModelBase::t3DModelData &Data)
      {
        std::ifstream infile(filename.c_str());
        std::string line = FindToken(infile, "Dx");
        double dx = 2.1 , dy = 3.4;
        infile >> dx >> dy;
        std::vector<double> Zthick;
        std::vector<double> Values;
        bool havelayer = true;
        double currvalue;
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
        std::cout << nx << " " << ny << " " << nz << std::endl;
        XCellSizes.resize(boost::extents[nx]);
        YCellSizes.resize(boost::extents[ny]);
        ZCellSizes.resize(boost::extents[nz]);
        std::copy(Zthick.begin(),Zthick.end(),ZCellSizes.begin());
        std::fill(XCellSizes.begin(),XCellSizes.end(),dx);
        std::fill(YCellSizes.begin(),YCellSizes.end(),dy);
        Data.resize(boost::extents[nx][ny][nz]);
        for (size_t i = 0; i < nx; ++i)
          for (size_t j = 0; j < ny; ++j)
            for (size_t k = 0; k < nz; ++k)
              Data[i][j][k] = Values.at((nx*ny)*k+ nx*j+i);

      }

  }//end of namespace jiba

