/*
 *  AFEKElement.h
 *  AdaptableFiniteElementKit
 *
 *  Copyright © 2022 MECV Software.  All rights reserved.
 */
#ifndef AFEKElement_h
#define AFEKElement_h

#import <Foundation/Foundation.h>
#import <simd/simd.h>

#ifdef __cplusplus
extern "C" {
#endif

// Options for specifying the data type used in the simulation.
typedef NS_ENUM(NSUInteger, AFEKElementGeometry)
{
    AFEKElementGeometryTriangular CF_SWIFT_NAME(triangular) = 1 << 0,
    AFEKElementGeometryQuadrilateral CF_SWIFT_NAME(quadrilateral) = 1 << 1,
    AFEKElementGeometryTetrahedral CF_SWIFT_NAME(tetrahedral) = 1 << 2,
    AFEKElementGeometryHexahedral CF_SWIFT_NAME(hexahedral) = 1 << 3,
};

/*
 *  A protocol describing an element.  Note that elements provided to AFEKMesh
 *  must further conform to AFEKElementSource2D or AFEKElementSource3D.
 */
@protocol AFEKElementSource <NSCopying, NSObject>

@required

/*
 *  The number of nodes in the element.
 */
-(NSUInteger) numberOfNodes;

/*
 *  The order of the shape functions.
 */
-(NSUInteger) order;

/*
 *  The number of dimensions of this element.  
 */
-(NSUInteger) dimensions;

/*
 *  The element geometry type.
 */
-(AFEKElementGeometry) geometry;

@end

/*
 *  A protocol describing a 3D volume element.
 */
@protocol AFEKElementSource3D <AFEKElementSource>

@required

/*
 *  Evaluate shape functions for a provided array of points in element local
 *  coordinates.
 *
 *  Inputs
 *      points          Local coordinates of points at which to evaluate the
 *                      shape functions.
 *
 *      numberOfPoints  The number of points in the points array.
 *
 *      shapeResults    The location at which to write values of the shape
 *                      function results.
 *
 *      localDiffShape1Results  The location at which to write values of the
 *                              local derivative of the shape functions with
 *                              respect to the first element local coordinate.
 *
 *      localDiffShape2Results  The location at which to write values of the
 *                              local derivative of the shape functions with
 *                              respect to the second element local coordinate.
 *
 *      localDiffShape3Results  The location at which to write values of the
 *                              local derivative of the shape functions with
 *                              respect to the third element local coordinate.
 *
 *      shapeRowStride          The stride between rows of the result arrays.
 *
 *  Output
 *
 *      Upon return the results of shape function evaluation (and derivatives)
 *      are in the provided result arrays.  Each array is sized to hold
 *      numberOfPoints x numberOfNodes values.
 */
-(void) generateShapeValuesAtPoints: (vector_double3 const* __nonnull) points
                     numberOfPoints: (NSUInteger) numberOfPoints
                       shapeResults: (double* __nonnull) shapeResults
             localDiffShape1Results: (double* __nonnull) localDiffShape1Results
             localDiffShape2Results: (double* __nonnull) localDiffShape2Results
             localDiffShape3Results: (double* __nonnull) localDiffShape3Results
                     shapeRowStride: (NSUInteger) shapeRowStride;
@end    // AFEKElementSource3D

/*
 *  A protocol describing a 2D surface element.
 */
@protocol AFEKElementSource2D <AFEKElementSource>

@required

/*
 *  Evaluate shape functions for a provided array of points in element local
 *  coordinates.
 *
 *  Inputs
 *      points          Local coordinates of points at which to evaluate the
 *                      shape functions.
 *
 *      numberOfPoints  The number of points in the points array.
 *
 *      shapeResults    The location at which to write values of the shape
 *                      function results.
 *
 *      localDiffShape1Results  The location at which to write values of the
 *                              local derivative of the shape functions with
 *                              respect to the first surface local coordinate.
 *
 *      localDiffShape2Results  The location at which to write values of the
 *                              local derivative of the shape functions with
 *                              respect to the second surface local coordinate.
 *
 *      shapeRowStride          The stride between rows of the result arrays.
 *
 *  Output
 *
 *      Upon return the results of shape function evaluation (and derivatives)
 *      are in the provided result arrays.  Each array is sized to hold
 *      numberOfPoints x numberOfNodes values.
 */
-(void) generateShapeValuesAtPoints: (vector_double2 const* __nonnull) points
                     numberOfPoints: (NSUInteger) numberOfPoints
                       shapeResults: (double* __nonnull) shapeResults
             localDiffShape1Results: (double* __nonnull) localDiffShape1Results
             localDiffShape2Results: (double* __nonnull) localDiffShape2Results
                     shapeRowStride: (NSUInteger) shapeRowStride;
@end

#ifdef __cplusplus
}
#endif

#endif /* AFEKElement_h */
