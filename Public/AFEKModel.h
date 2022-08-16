/*
 *  AFEKMesh.h
 *  AdaptableFiniteElementKit
 *
 *  Copyright © 2022 MECV Software.  All rights reserved.
 */
#ifndef AFEKModel_h
#define AFEKModel_h

#ifndef __METAL_VERSION__

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#import <simd/simd.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 *  A protocol describing a model to be used when simulating a AFEKMesh.
 *
 *  For models in which the force and stiffness matrices can be expressed
 *  as constant coefficients multiplied with shape function values and their
 *  derivatives.  Constant here refers to coefficients which do not depend
 *  on the particular nodal entry/entries in the force or stiffness arrays.
 *
 *  This requires only a description of the number of degrees of freedom a
 *  model is using.  One of AFEKModelSource3D or AFEKModelSource2D must be used
 *  to complete the description of the physical model.
 */
@protocol AFEKModelSource <NSCopying, NSObject>

@required

-(NSUInteger) degreesOfFreedom;

@optional

/*
 *  For models in which not all shape terms appear in the virtual work
 *  expression, this method can report which values are non-zero.  Doing so
 *  can allow AFEK to skip calculating non-zero coefficient terms.  Otherwise
 *  all terms will be used.
 */
-(bool) isShapeCoefficientNonZeroAtIndex: (NSUInteger) index;

/*
 *  If incompressible, return an incompressibility factor K where p = K * (deltaV - 1).
 */
-(double) kappa;

/*
 *  Methods to support GPU execution.
 */
/*
 *  Return the amount of space needed, in bytes, to hold the state results.
 *
 *  Inputs:
 *
 *      numberOfPoints:     The number of points for which state values should
 *                          be computed.
 *
 *      device:             The device on which the state is to be used on.
 *                          Should be the same device used to initialize the mesh.
 */
-(NSUInteger) stateSizeForNumberOfPoints: (NSUInteger) numberOfPoints device: (id<MTLDevice> __nonnull) device;

/*
 *  Return the MTLFunction corresponding to the coefficient compute function.
 *  The result is autoreleased.
 */
-(nullable id<MTLFunction>) getCoefficientComputeFunction;

/*
 *  Return a MTLBuffer holding arguments which make up the state.
 */
-(nullable id<MTLBuffer>) getArgumentBufferForState: (nullable id) state;

@end    // AFEKModelSource

@protocol AFEKModelSource3D <AFEKModelSource>

@optional

/*
 *  Given vectors u, du_dz1, du_dz2, du_dz3 where u is a vector of displacements
 *  for each degree of freedom.  z1, z2, and z3 are element local coordinate
 *  axes, and du_dz1, du_dz2, du_dz3 are derivatives of the displacements with
 *  respect to each of those directions respectively.  The i'th generalized
 *  displacement corresponds to the variable whose virtual displacement
 *  appeared multiplied against the i'th coefficient in the virtual work
 *  expression.
 *
 *  This method should return force and stiffness coefficients.
 *
 *  Force coefficients are returned in a 3 dimensional array of size
 *  4 x numPoints x degreesOfFreedom.  The 4 numPoints x degreesOfFreedom slices
 *  are densley packed together and have a stride between rows given by
 *  forceCoefficientsStride.  The order of these 4 slices corresponds to how
 *  the final force result is composed with the shape function values.  Slice 0
 *  is combined with shape values, whereas slices 1, 2, and 3 correspond to
 *  derivatives of shape values with respect to element local coordinates z1
 *  ,z2, and z3 respectively.
 *
 *  Stiffness coefficients are returned in a 5 dimensional array of size
 *  4 x 4 x numPoints x degreesOfFreedom x degreesOfFreedom.  Each
 *  numPoints x degreesOfFreedom x degreesOfFreedom array is regarded as
 *  flattened into a numPoints x (degreesOfFreedom**2) array with a stride
 *  between rows of stiffnessCoefficientsStride.  Consider this result as
 *  a 2 dimensional 4 x 4 array where each entry represents the flattened
 *  numPoints x (degreesOfFreedom**2) array.  Each entry is thus combined with
 *  products of shape values such that the i,j'th coefficient corresponds to
 *  the product of dS_dzi and dS_dzj, where S is the shape function and zi and
 *  zj are the local surface coordinates.  dS_dz0 is taken to be S itself.
 *
 *  All arrays are sized to hold data for numTotalPoints.  startPoint and
 *  numPoints represent the subset of points at which to compute results.
 *
 *  If a non-null state is provided it has been constructed with information
 *  from numTotalPoints.
 *
 *  Inputs:
 *      displacements           The vector of displacements, u.  The order of
 *                              generalized displacements in this vector is
 *                              consistent with the order that coefficients are
 *                              written to the result coefficient arrays, i.e.
 *                              u[i] corresponds to the i'th column of any
 *                              slice of the force coefficients array.
 *
 *      localDiffDisplacements1 The vector of derivative of displacements with
 *                              respect to the first surface coordinate, du_dz1.
 *
 *      localDiffDisplacements2 The vector of derivative of displacements with
 *                              respect to the second surface coordinate, du_dz2.
 *
 *      localDiffDisplacements3 The vector of derivative of displacements with
 *                              respect to the third surface coordinate, du_dz3.
 *
 *      displacementsStride     The stride between rows of the displacement
 *                              vectors.
 *
 *      forceCoefficients       Pointer to the beginning of the array of force
 *                              coefficients.
 *
 *      forceCoefficientsStride The stride between rows of forceCoefficients.
 *
 *      stiffnessCoefficients   Pointer to the beginning of the array of
 *                              stiffness coefficients.
 *
 *      stiffnessCoefficientsStride The stride between rows of stiffnessCoefficients.
 *
 *      state                   If non-null a state containing some
 *                              pre-computed data corresponding to the same
 *                              points as the displacements.
 *
 *      numTotalPoints          The total number of points for which the
 *                              displacements are supplied.  Also the number of
 *                              points which can be held by the result arrays
 *                              and the number of points at which the state
 *                              has been prepared.
 *
 *      startPoint              The index of the first point at which to begin
 *                              reading data and at which to begin writing
 *                              results.
 *
 *      numPoints               The number of points whose results are to be
 *                              computed.
 *
 *  Output:
 *
 *      Upon return forceCoefficients and stiffnessCoefficients will contain
 *      the computed results beginning at startPoint.
 */
// TODO: Move to function pointer, see below.  As a fcn ptr, reduces risk
// of user referencing state not contained in the actual state.
-(void) computeCoefficientsFromDisplacements: (double const* __nonnull) displacements
                     localDiffDisplacements1: (double const* __nonnull) localDiffDisplacements1
                     localDiffDisplacements2: (double const* __nonnull) localDiffDisplacements2
                     localDiffDisplacements3: (double const* __nonnull) localDiffDisplacements3
                         displacementsStride: (NSUInteger) displacementsStride
                           forceCoefficients: (double* __nonnull) forceCoefficients
                     forceCoefficientsStride: (NSUInteger) forceCoefficientsStride
                       stiffnessCoefficients: (double* __nonnull) stiffnessCoefficients
                 stiffnessCoefficientsStride: (NSUInteger) stiffnessCoefficientsStride
                                       state: (nullable id) state
                              numTotalPoints: (NSUInteger) numTotalPoints
                                  startPoint: (NSUInteger) startPoint
                                   numPoints: (NSUInteger) numPoints;

/*
 *  Compute applied force values at a point with a specified coordinate and
 *  local derivatives of that coordinate.
 *
 *  The result is a vector of length degreesOfFreedom.
 *
 *  Inputs:
 *      coordinate              The coordinates of the point at which force is
 *                              applied.
 *
 *      localDiffCoordinate1    The local derivative of the element coordinates
 *                              with respect to the first coordinate at the
 *                              point where coordinate is specified.
 *
 *      localDiffCoordinate2    The local derivative of the element coordinates
 *                              with respect to the second coordinate at the
 *                              point where coordinate is specified.
 *
 *      localDiffCoordinate3    The local derivative of the element coordinates
 *                              with respect to the third coordinate at the
 *                              point where coordinate is specified.
 *
 *      appliedForce            The applied force.
 *
 *      result                  The result array to place results in.
 *
 *  Output:
 *
 *      Upon return result contains the computed force values.
 */
-(void) computeForceAtCoordinate: (vector_double3) coordinate
            localDiffCoordinate1: (vector_double3) localDiffCoordinate1
            localDiffCoordinate2: (vector_double3) localDiffCoordinate2
            localDiffCoordinate3: (vector_double3) localDiffCoordinate3
                    appliedForce: (vector_double3) appliedForce
                          result: (double* __nonnull) result;

/*
 *  Prepare precomputed state necessary for computing force and stiffness
 *  coefficients.  Inputs are coordinates and their derivatives with respect to
 *  element local coordinates.
 *
 *  Each array contains numberOfPoints points of data for which the state
 *  is to be prepared.
 *
 *  Inputs
 *      coordinates             A vector of coordinates.
 *
 *      localDiffCoordinates1   A vector of derivatives of coordinates with
 *                              respect to the first element coordinate.
 *
 *      localDiffCoordinates2   A vector of derivatives of coordinates with
 *                              respect to the second element coordinate.
 *
 *      localDiffCoordinates3   A vector of derivatives of coordinates with
 *                              respect to the second element coordinate.
 *
 *      numberOfPoints          The number of points for which input coordinates
 *                              and normals, etc, are specified.
 *
 *  Output:
 *
 *      Returns an autoreleased object containing any necessary state to be
 *      used with computeCoefficientsFromDisplacements:...
 */
-(nullable id) prepareStateForCoordinates: (vector_double3 const* __nonnull) coordinates
                    localDiffCoordinates1: (vector_double3 const* __nonnull) localDiffCoordinates1
                    localDiffCoordinates2: (vector_double3 const* __nonnull) localDiffCoordinates2
                    localDiffCoordinates3: (vector_double3 const* __nonnull) localDiffCoordinates3
                          numberOfPoints: (NSUInteger) numberOfPoints;

-(nullable id) prepareStateForCoordinates: (vector_double3 const* __nonnull) coordinates
                   localDiffCoordinates1: (vector_double3 const* __nonnull) localDiffCoordinates1
                   localDiffCoordinates2: (vector_double3 const* __nonnull) localDiffCoordinates2
                   localDiffCoordinates3: (vector_double3 const* __nonnull) localDiffCoordinates3
                         numberOfPoints: (NSUInteger) numberOfPoints
                                heap: (id<MTLHeap> __nullable) heap;

-(void) computeFloatCoefficientsFromDisplacements: (float const* __nonnull) displacements
                    localDiffDisplacements1: (float const* __nonnull) localDiffDisplacements1
                    localDiffDisplacements2: (float const* __nonnull) localDiffDisplacements2
                    localDiffDisplacements3: (float const* __nonnull) localDiffDisplacements3
                        displacementsStride: (NSUInteger) displacementsStride
                          forceCoefficients: (float* __nonnull) forceCoefficients
                    forceCoefficientsStride: (NSUInteger) forceCoefficientsStride
                      stiffnessCoefficients: (float* __nonnull) stiffnessCoefficients
                stiffnessCoefficientsStride: (NSUInteger) stiffnessCoefficientsStride
                                      state: (nullable id) state
                             numTotalPoints: (NSUInteger) numTotalPoints
                                 startPoint: (NSUInteger) startPoint
                                  numPoints: (NSUInteger) numPoints;
@end    // AFEKModelSource3D

/*
 *  AFEKModelSource2D describes a model of a shell.
 */
@protocol AFEKModelSource2D <AFEKModelSource>

@optional

/*
 *  Given vectors u, du_dz1, du_dz2 where u is a vector of displacements for
 *  each degree of freedom.  z1 and z2 are surface local curvilinear coordinate
 *  axes, and du_dz1 and du_dz2 are derivatives of the displacements with
 *  respect to each of those directions respectively.  The i'th generalized
 *  displacement corresponds to the variable whose virtual displacement
 *  appeared multiplied against the i'th coefficient in the virtual work
 *  expression.
 *
 *  This method should return force and stiffness coefficients.
 *
 *  Force coefficients are returned in a 3 dimensional array of size
 *  3 x numPoints x degreesOfFreedom.  The 3 numPoints x degreesOfFreedom slices
 *  are densley packed together and have a stride between rows given by
 *  forceCoefficientsStride.  The order of these 3 slices corresponds to how
 *  the final force result is composed with the shape function values.  Slice 0
 *  is combined with shape values, whereas slices 1 and 2 correspond to
 *  derivatives of shape values with respect to surface local coordinates z1
 *  and z2 respectively.
 *
 *  Stiffness coefficients are returned in a 5 dimensional array of size
 *  3 x 3 x numPoints x degreesOfFreedom x degreesOfFreedom.  Each
 *  numPoints x degreesOfFreedom x degreesOfFreedom array is regarded as
 *  flattened into a numPoints x (degreesOfFreedom**2) array with a stride
 *  between rows of stiffnessCoefficientsStride.  Consider this result as
 *  a 2 dimensional 3 x 3 array where each entry represents the flattened
 *  numPoints x (degreesOfFreedom**2) array.  Each entry is thus combined with
 *  products of shape values such that the i,j'th coefficient corresponds to
 *  the product of dS_dzi and dS_dzj, where S is the shape function and zi and
 *  zj are the local surface coordinates.  dS_dz0 is taken to be S itself.
 *
 *  All arrays are sized to hold data for numTotalPoints.  startPoint and
 *  numPoints represent the subset of points at which to compute results.
 *
 *  If a non-null state is provided it has been constructed with information
 *  from numTotalPoints.
 *
 *  Inputs:
 *      displacements           The vector of displacements, u.  The order of
 *                              generalized displacements in this vector is
 *                              consistent with the order that coefficients are
 *                              written to the result coefficient arrays, i.e.
 *                              u[i] corresponds to the i'th column of any
 *                              slice of the force coefficients array.
 *
 *      localDiffDisplacements1 The vector of derivative of displacements with
 *                              respect to the first surface coordinate, du_dz1.
 *
 *      localDiffDisplacements2 The vector of derivative of displacements with
 *                              respect to the second surface coordinate, du_dz2.
 *
 *      displacementsStride     The stride between rows of the displacement
 *                              vectors.
 *
 *      forceCoefficients       Pointer to the beginning of the array of force
 *                              coefficients.
 *
 *      forceCoefficientsStride The stride between rows of forceCoefficients.
 *
 *      stiffnessCoefficients   Pointer to the beginning of the array of
 *                              stiffness coefficients.
 *
 *      stiffnessCoefficientsStride The stride between rows of stiffnessCoefficients.
 *
 *      state                   If non-null a state containing some
 *                              pre-computed data corresponding to the same
 *                              points as the displacements.
 *
 *      numTotalPoints          The total number of points for which the
 *                              displacements are supplied.  Also the number of
 *                              points which can be held by the result arrays
 *                              and the number of points at which the state
 *                              has been prepared.
 *
 *      startPoint              The index of the first point at which to begin
 *                              reading data and at which to begin writing
 *                              results.
 *
 *      numPoints               The number of points whose results are to be
 *                              computed.
 *
 *  Output:
 *
 *      Upon return forceCoefficients and stiffnessCoefficients will contain
 *      the computed results beginning at startPoint.
 */
-(void) computeCoefficientsFromDisplacements: (double const* __nonnull) displacements
                     localDiffDisplacements1: (double const* __nonnull) localDiffDisplacements1
                     localDiffDisplacements2: (double const* __nonnull) localDiffDisplacements2
                         displacementsStride: (NSUInteger) displacementsStride
                           forceCoefficients: (double* __nonnull) forceCoefficients
                     forceCoefficientsStride: (NSUInteger) forceCoefficientsStride
                       stiffnessCoefficients: (double* __nonnull) stiffnessCoefficients
                 stiffnessCoefficientsStride: (NSUInteger) stiffnessCoefficientsStride
                                       state: (nullable id) state
                              numTotalPoints: (NSUInteger) numTotalPoints
                                  startPoint: (NSUInteger) startPoint
                                   numPoints: (NSUInteger) numPoints;

/*
 *  Compute applied force values at a point with a specified coordinate and
 *  local derivatives of that coordinate, as well as the surface normal
 *  at that point and its associated local derivatives.
 *
 *  The result is a vector of length degreesOfFreedom.
 *
 *  Inputs:
 *      coordinate              The coordinates of the point at which force is
 *                              applied.
 *
 *      localDiffCoordinate1    The local derivative of the surface coordinates
 *                              with respect to the first surface coordinate
 *                              at the point where coordinate is specified.
 *
 *      localDiffCoordinate2    The local derivative of the surface coordinates
 *                              with respect to the second surface coordinate
 *                              at the point where coordinate is specified.
 *
 *      normal                  The surface normal at the point at which force is
 *                              applied.
 *
 *      localDiffNormal1        The local derivative of the surface normal
 *                              with respect to the first surface coordinate
 *                              at the point where coordinate is specified.
 *
 *      localDiffNormal2        The local derivative of the surface normal
 *                              with respect to the second surface coordinate
 *                              at the point where coordinate is specified.
 *
 *      appliedForce            The applied force.
 *
 *      result                  The result array to place results in.
 *
 *  Output:
 *
 *      Upon return result contains the computed force values.
 */
-(void) computeForceAtCoordinate: (vector_double3) coordinate
            localDiffCoordinate1: (vector_double3) localDiffCoordinate1
            localDiffCoordinate2: (vector_double3) localDiffCoordinate2
                          normal: (vector_double3) normal
                localDiffNormal1: (vector_double3) localDiffNormal1
                localDiffNormal2: (vector_double3) localDiffNormal2
                    appliedForce: (vector_double3) appliedForce
                          result: (double* __nonnull) result;

/*
 *  Prepare precomputed state necessary for computing force and stiffness
 *  coefficients.  Inputs are coordinates, normals, and their derivatives
 *  with respect to the surface local curvilinear coordinates.
 *
 *  Each array contains numberOfPoints points of data for which the state
 *  is to be prepared.
 *
 *  Inputs
 *      coordinates             A vector of coordinates.
 *
 *      localDiffCoordinates1   A vector of derivatives of coordinates with
 *                              respect to the first surface coordinate.
 *
 *      localDiffCoordinates2   A vector of derivatives of coordinates with
 *                              respect to the second surface coordinate.
 *
 *      normals                 A vector of surface mid-plane normals.
 *
 *      localDiffNormals1       A vector of derivatives of normals with respect
 *                              to the first surface coordinate.
 *
 *      localDiffNormals2       A vector of derivatives of normals with respect
 *                              to the second surface coordinate.
 *
 *      numberOfPoints          The number of points for which input coordinates
 *                              and normals, etc, are specified.
 *
 *  Output:
 *
 *      Returns an autoreleased object containing any necessary state to be
 *      used with computeCoefficientsFromDisplacements:...
 */
-(nullable id) prepareStateForCoordinates: (vector_double3 const* __nonnull) coordinates
                    localDiffCoordinates1: (vector_double3 const* __nonnull) localDiffCoordinates1
                    localDiffCoordinates2: (vector_double3 const* __nonnull) localDiffCoordinates2
                                  normals: (vector_double3 const* __nonnull) normals
                        localDiffNormals1: (vector_double3 const* __nonnull) localDiffNormals1
                        localDiffNormals2: (vector_double3 const* __nonnull) localDiffNormals2
                           numberOfPoints: (NSUInteger) numberOfPoints;

-(nullable id) prepareStateForCoordinates: (vector_double3 const* __nonnull) coordinates
                   localDiffCoordinates1: (vector_double3 const* __nonnull) localDiffCoordinates1
                   localDiffCoordinates2: (vector_double3 const* __nonnull) localDiffCoordinates2
                                 normals: (vector_double3 const* __nonnull) normals
                       localDiffNormals1: (vector_double3 const* __nonnull) localDiffNormals1
                       localDiffNormals2: (vector_double3 const* __nonnull) localDiffNormals2
                          numberOfPoints: (NSUInteger) numberOfPoints
                                    heap: (id<MTLHeap> __nullable) heap;

// TODO: Temporary

-(void) computeFloatCoefficientsFromDisplacements: (float const* __nonnull) displacements
                    localDiffDisplacements1: (float const* __nonnull) localDiffDisplacements1
                    localDiffDisplacements2: (float const* __nonnull) localDiffDisplacements2
                        displacementsStride: (NSUInteger) displacementsStride
                          forceCoefficients: (float* __nonnull) forceCoefficients
                    forceCoefficientsStride: (NSUInteger) forceCoefficientsStride
                      stiffnessCoefficients: (float* __nonnull) stiffnessCoefficients
                stiffnessCoefficientsStride: (NSUInteger) stiffnessCoefficientsStride
                                      state: (nullable id) state
                             numTotalPoints: (NSUInteger) numTotalPoints
                                 startPoint: (NSUInteger) startPoint
                                  numPoints: (NSUInteger) numPoints;

@end    // AFEKModelSource2D

// TODO: make function pointer type declaration.
//-(void) AFEKModelSource3DComputeCoefficients(double const* __nonnull displacements,
//                                             double const* __nonnull localDiffDisplacements1,
//                                             double const* __nonnull localDiffDisplacements2,
//                                             double const* __nonnull localDiffDisplacements3,
//                                             NSUInteger displacementsStride,
//                                             double* __nonnull forceCoefficients,
//                                             NSUInteger forceCoefficientsStride,
//                                             double* __nonnull stiffnessCoefficients,
//                                             NSUInteger stiffnessCoefficientsStride,
//                                             nullable id state,
//                                             NSUInteger numTotalPoints,
//                                             NSUInteger startPoint,
//                                             NSUInteger numPoints);
typedef void (*AFEKModelSource3DComputeCoefficients)(double const* __nonnull displacements,
                                                     double const* __nonnull localDiffDisplacements1,
                                                     double const* __nonnull localDiffDisplacements2,
                                                     double const* __nonnull localDiffDisplacements3,
                                                     NSUInteger displacementsStride,
                                                     double* __nonnull forceCoefficients,
                                                     NSUInteger forceCoefficientsStride,
                                                     double* __nonnull stiffnessCoefficients,
                                                     NSUInteger stiffnessCoefficientsStride,
                                                     void* __nonnull state,
                                                     NSUInteger numTotalPoints,
                                                     NSUInteger startPoint,
                                                     NSUInteger numPoints);


//typedef void (*AFEKModelSource2DComputeCoefficients)(double const* __nonnull displacements,
//                                                     double const* __nonnull localDiffDisplacements1,
//                                                     double const* __nonnull localDiffDisplacements2,
//                                                     NSUInteger displacementsStride,
//                                                     double* __nonnull forceCoefficients,
//                                                     NSUInteger forceCoefficientsStride,
//                                                     double* __nonnull stiffnessCoefficients,
//                                                     NSUInteger stiffnessCoefficientsStride,
//                                                     void* __nonnull state,
//                                                     NSUInteger numTotalPoints,
//                                                     NSUInteger startPoint,
//                                                     NSUInteger numPoints);
//
//Also return coefficients for the mean-dilation incompressibility terms.
//typedef void (*AFEKModelSource3DComputeCoefficientsIncompressible)(double const* __nonnull displacements,
//                                                     double const* __nonnull localDiffDisplacements1,
//                                                     double const* __nonnull localDiffDisplacements2,
//                                                     double const* __nonnull localDiffDisplacements3,
//                                                     NSUInteger displacementsStride,
//                                                     double* __nonnull forceCoefficients,
//                                                     NSUInteger forceCoefficientsStride,
//                                                     double* __nonnull stiffnessCoefficients,
//                                                     NSUInteger stiffnessCoefficientsStride,
//                                                                   double* __nonnull incompressibilityCoefficients,
//                                                                   NSUInteger incompressibilityCoefficientsStride,
//                                                     void* __nonnull state,
//                                                     NSUInteger numTotalPoints,
//                                                     NSUInteger startPoint,
//                                                     NSUInteger numPoints);

#ifdef __cplusplus
}
#endif

#else   // __METAL_VERSION__
#include <AdaptableFiniteElementKit/AdaptableFiniteElementKit.h>

/*
 *  A Metal compute shader function to evaluate force and stiffness coefficients.
 */
typedef struct AFEKModelStateStruct device* AFEKModelState;

using AFEKModelSource3DComputeCoefficientsType = void(device float const*    displacements,
                                                      device float const*    localDiffDisplacements1,
                                                      device float const*    localDiffDisplacements2,
                                                      device float const*    localDiffDisplacements3,
                                                      uint                   displacementsStride,
                                                      device float*          forceCoefficients,
                                                      uint                   forceCoefficientsStride,
                                                      uint                   forceCoefficientsSliceStride,
                                                      device float*          stiffnessCoefficients,
                                                      uint                   stiffnessCoefficientsStride,
                                                      uint                   stiffnessCoefficientsSliceStride,
                                                      AFEKModelState         state,
                                                      uint                   numTotalPoints,
                                                      uint                   startPoint,
                                                      uint                   numPoints);

// TODO: Change so exact function name isn't required, i.e. try to get the name
// from the returned MTLFunction or similar...  Looks like use a visible function
// table.  As this stands a given program can only have this defined once.
[[visible]] void AFEKModelSource3DComputeCoefficients(device float const*    displacements,
                                                      device float const*    localDiffDisplacements1,
                                                      device float const*    localDiffDisplacements2,
                                                      device float const*    localDiffDisplacements3,
                                                      uint                   displacementsStride,
                                                      device float*          forceCoefficients,
                                                      uint                   forceCoefficientsStride,
                                                      uint                   forceCoefficientsSliceStride,
                                                      device float*          stiffnessCoefficients,
                                                      uint                   stiffnessCoefficientsStride,
                                                      uint                   stiffnessCoefficientsSliceStride,
                                                      AFEKModelState         state,
                                                      uint                   numTotalPoints,
                                                      uint                   startPoint,
                                                      uint                   numPoints);
#endif  // __METAL_VERSION__

#endif /* AFEKModel_h */
