//============================================================================
// Name        : CalcRawRealization.h
// Author      : 5 Dec 2014
// Version     :
// Copyright   : 2014, mm489
//============================================================================

#ifndef CALCRAWREALIZATION_H_
#define CALCRAWREALIZATION_H_

#ifdef HAVEHPX
#include <hpx/config.hpp>
#endif
#include <boost/python.hpp>
#include <string>
#include "../Global/VecMat.h"
#include "../Global/Jif3DGlobal.h"

J3DEXPORT double CalcRawRealization(int nx, int ny, int nz, double delta, double topthick,
    double freq, double bgcond, boost::python::list& Conductivities, std::string tempdir,
    std::string x3dname);

#endif
