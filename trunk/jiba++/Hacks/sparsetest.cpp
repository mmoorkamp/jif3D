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
    jiba::rsparse SparseMat(n, n), Inverse(n, n);
    for (size_t i = 0; i < n; ++i)
      {
        SparseMat.push_back(i, i, 4);
      }

    jiba::InvertMatrix(SparseMat, Inverse);
  }

