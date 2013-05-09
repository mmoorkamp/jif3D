//============================================================================
// Name        : ModelTransforms.h
// Author      : Apr 15, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#ifndef MODELTRANSFORMS_H_
#define MODELTRANSFORMS_H_

#include "../ModelTransforms/ChainedTransform.h"
#include "../ModelTransforms/ConductivityTransform.h"
#include "../ModelTransforms/DensityTransform.h"
#include "../ModelTransforms/LogTransform.h"
#include "../ModelTransforms/ModelCopyTransform.h"
#include "../ModelTransforms/MultiSectionTransform.h"
#include "../ModelTransforms/NormalizeTransforms.h"
#include "../ModelTransforms/TanhTransform.h"
#include "../ModelTransforms/WaveletModelTransform.h"
#include "ModelDistributor.h"

/*! \file ModelTransforms.h
 * This file contains includes header files for various classes to transform model parameters within an inversion.
 * These are used to either make the inversion more well behaved or to calculate
 * one physical quantity from another, like density from slowness in LogDensityTransform.
 */


#endif /* MODELTRANSFORMS_H_ */
