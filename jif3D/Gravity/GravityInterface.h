//============================================================================
// Name        : GravityInterface.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#ifndef GRAVITYINTERFACE_H_
#define GRAVITYINTERFACE_H_

/*! \file GravityInterface.h
 * Here we declare the interface to the 3D Gravity code for use with R. This has to be done
 * by declaring C-style functions that pass all parameters by reference. As we cannot determine
 * the size of vectorial quantities we cannot perform any consistency checking whatsoever. It is
 * completely up to the caller to ensure that all structures have been allocated and the dimensions match.
 *
 * Also we us a static object internally to perform these calculations. This means that these functions are not
 * thread safe !.
 */

extern "C"
  {

    /** \addtogroup gravity Gravity forward modeling, display and inversion  */
    /* @{ */

    //! This function allocates storage for the model and has to be called before any calculations
    /*! In order to allow to specify the storage of sensitivity matrices and the associated accelerations,
     * we allocate the model object dynamically through this function.
     * @param storescalar If > 0 we store the sensitivites for scalar calculations. Increases memory usage, but accelerates calculations
     * @param storetensor If > 0 we store the sensitivites for tensor calculations. Increases memory usage, but accelerates calculations
     */
    void AllocateModel(const int *storescalar, const int *storetensor);

    //! Perform a forward calculation for scalar gravity data
    /*! Calculate the scalar gravity response of the specified model
     * @param XSizes The sizes of the model cells in x-direction in m
     * @param nx The number of cells in x-direction
     * @param YSizes The sizes of the model cells in y-direction in m
     * @param ny The number of cells in y-direction
     * @param ZSizes The sizes of the model cells in z-direction in m
     * @param nz The number of cells in z-direction
     * @param Densities The densities of all model cells in g/cm^3 the size of this vector must be nx * ny * nz. The z-direction varies fastest.
     * @param BG_Thicknesses An array containing the thickness in m for each background layer
     * @param BG_Densities An array containing the density in g/cm^3 for each background layer
     * @param nbg_layers The number of background layers
     * @param XMeasPos The vector of x-coordinates of the measurement points in m
     * @param YMeasPos The vector of y-coordinates of the measurement points in m
     * @param ZMeasPos The vector of z-coordinates of the measurement points in m
     * @param nmeas The number of measurement points. Must match the size of all three vectors above
     * @param GravAcceleration The calculated gravitational acceleration in m/s^2. Must be allocated with nmeas elements
     */
    void CalcScalarForward(const double *XSizes, const unsigned int *nx,
        const double *YSizes, const unsigned int *ny, const double *ZSizes,
        const unsigned int *nz, const double *Densities,
        const double *BG_Thicknesses, const double *BG_Densities, const unsigned int *nbg_layers,
        const double *XMeasPos, const double *YMeasPos, const double *ZMeasPos,
        const unsigned int *nmeas, double *GravAcceleration);

    //! Perform a forward calculation for tensor gravity data
    /*! Calculate the tensor gravity response of the specified model
     * @param XSizes The sizes of the model cells in x-direction in m
     * @param nx The number of cells in x-direction
     * @param YSizes The sizes of the model cells in y-direction in m
     * @param ny The number of cells in y-direction
     * @param ZSizes The sizes of the model cells in z-direction in m
     * @param nz The number of cells in z-direction
     * @param Densities The densities of all model cells in g/cm^3 the size of this vector must be nx * ny * nz. The z-direction varies fastest.
     * @param BG_Thicknesses An array containing the thickness in m for each background layer
     * @param BG_Densities An array containing the density in g/cm^3 for each background layer
     * @param nbg_layers The number of background layers
     * @param XMeasPos The vector of x-coordinates of the measurement points in m
     * @param YMeasPos The vector of y-coordinates of the measurement points in m
     * @param ZMeasPos The vector of z-coordinates of the measurement points in m
     * @param nmeas The number of measurement points. Must match the size of all three vectors above
     * @param GravTensor The calculated FTG tensor in vectorial form. Must be allocated with 9*nmeas elements.
     * Each block of 9 represents one FTG tensor at a measurement point.
     */
    void CalcTensorForward(const double *XSizes, const unsigned int *nx,
        const double *YSizes, const unsigned int *ny, const double *ZSizes,
        const unsigned int *nz, const double *Densities,
        const double *BG_Thicknesses, const double *BG_Densities, const unsigned int *nbg_layers,
        const double *XMeasPos, const double *YMeasPos, const double *ZMeasPos,
        const unsigned int *nmeas, double *GravTensor);

  /* @} */
  }
#endif /* GRAVITYINTERFACE_H_ */
