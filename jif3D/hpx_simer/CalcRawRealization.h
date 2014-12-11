//============================================================================
// Name        : CalcRawRealization.h
// Author      : 5 Dec 2014
// Version     : 
// Copyright   : 2014, mm489
//============================================================================

#ifndef CALCRAWREALIZATION_H_
#define CALCRAWREALIZATION_H_

#include <vector>
#include <string>
#include "../Global/VecMat.h"

jif3D::rvec CalcRawRealization(int nx, int ny, int nz, double delta, double topthick,
    double freq, double bgcond, std::vector<double> Conductivities, std::string tempdir,
    std::string x3dname);

#endif
