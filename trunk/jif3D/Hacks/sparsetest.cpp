//============================================================================
// Name        : sparsetest.cpp
// Author      : Sep 29, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================

#include "../Global/VecMat.h"
#include "../Inversion/MatrixTools.h"
int main()
  {
    const size_t n = 100;
    jif3D::rsparse SparseMat(n, n), Inverse(n, n);
    for (size_t i = 0; i < n; ++i)
      {
        SparseMat(i, i) =  4.0;
      }

    jif3D::InvertMatrix(SparseMat, Inverse);
  }

