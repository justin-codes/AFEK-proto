/*
 *  AFEKMesh.h
 *  AdaptableFiniteElementKit
 *
 *  Copyright © 2022 MECV Software.  All rights reserved.
 */

#ifndef AFEKMesh_h
#define AFEKMesh_h

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#import <AdaptableFiniteElementKit/AFEKElement.h>
#import <AdaptableFiniteElementKit/AFEKModel.h>
#import <simd/simd.h>

#ifdef __cplusplus
extern "C" {
#endif

// Options for specifying the data type used in the simulation.
typedef NS_ENUM(NSUInteger, AFEKDataType)
{
    AFEKDataTypeFloat64 CF_SWIFT_NAME(float64) = 1 << 0,
    AFEKDataTypeFloat32 CF_SWIFT_NAME(float32) = 1 << 1,
};

// Options for specifying the integration method used when assembling
// the stiffness terms.
typedef NS_ENUM(NSUInteger, AFEKIntegrationMethod)
{
    // Integrate using full Gaussian quadrature.
    AFEKIntegrationMethodGaussFull CF_SWIFT_NAME(gaussFull) = 1,
};

/*
 *  AFEKMesh defines the object mesh.
 */
@interface AFEKMesh : NSObject

// The integration method currently being used by this mesh.
@property (readonly, nonatomic) AFEKIntegrationMethod integrationMethod;

// The number of elements in this mesh.
@property (readonly, nonatomic) NSUInteger numberOfElements;

// The number of vertices in this mesh.
@property (readonly, nonatomic) NSUInteger numberOfVertices;

// The data type this mesh has been prepared for.
@property (readonly, nonatomic) AFEKDataType dataType;

/*
 *  Initialize a new surface mesh with the following parameters.
 *
 *  Inputs:
 *      elements        An array specifying which vertices make up each
 *                      element in the mesh.  This is an array of size M x N
 *                      where M is the number of elements in the mesh and N
 *                      is the number of vertices in each element.
 *
 *                      Note that the order in which vertices are specified
 *                      for a given element must match the order of the shape
 *                      values provided by the AFEKElementSource input object.
 *
 *      elementStride   The stride between consecutive rows of the 'elements'
 *                      array.
 *
 *      elementType     An object which conforms to AFEKElementSource2D that
 *                      defines the type of element to use for this mesh.
 *
 *      coordinates     An array of coordinate values for each vertex.  It is
 *                      a densely packed array of size 'numberOfVertices' x 3.
 *
 *      normals         An arry of normal vectors for each vertex.  It is a
 *                      densely packed array of size 'numberOfVertices' x 3.
 *
 *      numberOfElements    The number of elements in the mesh.
 *
 *      numberOfVertices    The number of vertices in the mesh.
 *
 *      model           An object which conforms to AFEKModelSource2D that
 *                      defines the physical model used for the simulation.
 *
 *      dataType        The data type to be used in the simulation.
 */
-(nonnull instancetype) initWithElements: (NSUInteger const* __nonnull) elements
                           elementStride: (NSUInteger) elementStride
                             elementType: (const id<AFEKElementSource2D> __nonnull) elementType
                             coordinates: (double const* __nonnull) coordinates
                                 normals: (double const* __nonnull) normals
                        numberOfElements: (NSUInteger) numberOfElements
                        numberOfVertices: (NSUInteger) numberOfVertices
                                   model: (const id<AFEKModelSource2D> __nonnull) model
                                dataType: (AFEKDataType) dataType;

/*
 *  Initialize a new surface mesh (same as above) for execution on a GPU.
 *
 *  Inputs:
 *      ...
 *
 *      device      The device on which to prepare data, execute iterations,
 *                  and produce results.
 */
-(nonnull instancetype) initWithElements: (NSUInteger const* __nonnull) elements
                          elementStride: (NSUInteger) elementStride
                            elementType: (const id<AFEKElementSource2D> __nonnull) elementType
                            coordinates: (double const* __nonnull) coordinates
                                normals: (double const* __nonnull) normals
                       numberOfElements: (NSUInteger) numberOfElements
                       numberOfVertices: (NSUInteger) numberOfVertices
                                  model: (const id<AFEKModelSource2D> __nonnull) model
                                dataType: (AFEKDataType) dataType
                                  device: (id<MTLDevice> __nonnull) device;

/*
 *  Initialize a new volumetric mesh with the following parameters.
 *
 *  Inputs:
 *      elements        An array specifying which vertices make up each
 *                      element in the mesh.  This is an array of size M x N
 *                      where M is the number of elements in the mesh and N
 *                      is the number of vertices in each element.
 *
 *                      Note that the order in which vertices are specified
 *                      for a given element must match the order of the shape
 *                      values provided by the AFEKElementSource input object.
 *
 *      elementStride   The stride between consecutive rows of the 'elements'
 *                      array.
 *
 *      elementType     An object which conforms to AFEKElementSource3D that
 *                      defines the type of element to use for this mesh.
 *
 *      coordinates     An array of coordinate values for each vertex.  It is
 *                      a densely packed array of size 'numberOfVertices' x 3.
 *
 *      normals         An arry of normal vectors for each vertex.  It is a
 *                      densely packed array of size 'numberOfVertices' x 3.
 *
 *      numberOfElements    The number of elements in the mesh.
 *
 *      numberOfVertices    The number of vertices in the mesh.
 *
 *      model           An object which conforms to AFEKModelSource3D that
 *                      defines the physical model used for the simulation.
 *
 *      dataType        The data type to be used in the simulation.
 */
-(nonnull instancetype) initWithElements: (NSUInteger const* __nonnull) elements
                           elementStride: (NSUInteger) elementStride
                             elementType: (const id<AFEKElementSource3D> __nonnull) elementType
                             coordinates: (double const* __nonnull) coordinates
                        numberOfElements: (NSUInteger) numberOfElements
                        numberOfVertices: (NSUInteger) numberOfVertices
                                   model: (const id<AFEKModelSource3D> __nonnull) model
                                dataType: (AFEKDataType) dataType;

-(nonnull instancetype) initWithElements: (NSUInteger const* __nonnull) elements
                           elementStride: (NSUInteger) elementStride
                             elementType: (const id<AFEKElementSource3D> __nonnull) elementType
                             coordinates: (double const* __nonnull) coordinates
                        numberOfElements: (NSUInteger) numberOfElements
                        numberOfVertices: (NSUInteger) numberOfVertices
                                   model: (const id<AFEKModelSource3D> __nonnull) model
           coefficientEvaluationFunction: (AFEKModelSource3DComputeCoefficients __nonnull) coefficientEvaluationFunction
                                dataType: (AFEKDataType) dataType;

/*
 *  Initialize a new volumetric mesh (same as above) for execution on a GPU.
 *
 *  Inputs:
 *      ...
 *
 *      device      The device on which to prepare data, execute iterations,
 *                  and produce results.
 */
-(nonnull instancetype) initWithElements: (NSUInteger const* __nonnull) elements
                          elementStride: (NSUInteger) elementStride
                            elementType: (const id<AFEKElementSource3D> __nonnull) elementType
                            coordinates: (double const* __nonnull) coordinates
                       numberOfElements: (NSUInteger) numberOfElements
                       numberOfVertices: (NSUInteger) numberOfVertices
                                  model: (const id<AFEKModelSource3D> __nonnull) model
                                dataType: (AFEKDataType) dataType
                                  device: (id<MTLDevice> __nonnull) device;

/*
 *  Completely fix a given node at its initialized location.  The displacement
 *  of each degree of freedom is fixed to zero.
 */
-(void) fixNode: (NSUInteger) node;

/*
 *  Apply a point force to a given location.  The location is specified by
 *  the barycentric coordinates within a given element.
 */
-(NSUInteger) applyPointForce: (vector_double3) force
                      element: (NSUInteger) element
                     location: (vector_double3) location;

-(NSUInteger) applySurfaceForce: (vector_double3) force
                        element: (NSUInteger) element
                           face: (NSUInteger) face;

-(NSUInteger) applyVolumeForce: (vector_double3) force
                       element: (NSUInteger) element;

/*
 *  Methods for executing iterations on the CPU.
 */
// Execute a single iteration of the simulation using the CPU.  Only valid
// when called from AFEKMesh objects initialized for CPU execution.
-(void) runIteration;

// Execute a single iteration asynchronously.  Upon completion the
// completionBlock, if non-null, will execute with a pointer to the solution
// data and the L2 norm error.  The callback may execute on a thread different
// than on which this method was called; the provided data will be guaranteed to
// be available but any other data may not be.  The user should take care to
// ensure that accesses are properly synchronized and ordered.  Only valid
// when called from AFEKMesh objects initialized for CPU execution.
-(void) runIterationWithCompletionHandler: (nullable void(^)(void const* __nonnull, void const* __nonnull)) completionBlock;

//-(void) runMaximumIterationsWithCompletionHandler: (nullable void(^)(double const* __nonnull, double)) completionBlock;

// Wait for pending iterations to complete.  Only valid when called from
// AFEKMesh objects initialized for CPU execution.
-(void) waitUntilCompleted;

-(void) setDisplacements: (float const* __nonnull) displacements;

/*
 *  Methods for executing iterations on the GPU.
 */

-(void) setErrorBuffer: (id<MTLBuffer> __nullable) errorBuffer
                offset: (NSUInteger) offset;

-(void) setDisplacementBuffer: (id<MTLBuffer> __nullable) displacementBuffer;

//-(NSUInteger) sizeForElements: (NSUInteger const* __nonnull) elements
//                elementStride: (NSUInteger) elementStride
//                  elementType: (const id<AFEKElementSource> __nonnull) elementType
//             numberOfElements: (NSUInteger) numberOfElements
//             numberOfVertices: (NSUInteger) numberOfVertices
//                     dataType: (AFEKDataType) dataType
//                       device: (id<MTLDevice> __nonnull) device;

// TODO: Use below method as not only is the predicate possibly determined
// via GPU feedback but the force as well.  Possibly location, etc, as well?
-(void) encodeApplyForce: (vector_double3) force
                 element: (NSUInteger) element
                location: (vector_double3) location
         predicateBuffer: (id<MTLBuffer> __nullable) predicateBuffer
         predicateOffset: (NSUInteger) predicateOffset
           commandBuffer: (id<MTLCommandBuffer> __nonnull) commandBuffer;



//-(void) encodeApplyPointForce: (id<MTLBuffer> __nonnull) force
//             forceOffset: (NSUInteger) forceOffset
//                 element: (NSUInteger) element
//                location: (vector_double3) location
//         predicateBuffer: (id<MTLBuffer> __nullable) predicateBuffer
//         predicateOffset: (NSUInteger) predicateOffset
//           commandBuffer: (id<MTLCommandBuffer> __nonnull) commandBuffer;

//-(void) encodeApplySurfaceForce: (id<MTLBuffer> __nonnull) force
//             forceOffset: (NSUInteger) forceOffset
//                 element: (NSUInteger) element
//                face: (NSUInteger) face
//         predicateBuffer: (id<MTLBuffer> __nullable) predicateBuffer
//         predicateOffset: (NSUInteger) predicateOffset
//           commandBuffer: (id<MTLCommandBuffer> __nonnull) commandBuffer;

-(void) encodeMaximumIterations: (NSUInteger) maximumIterations
                 errorThreshold: (id<MTLBuffer> __nullable) errorThreshold
                  commandBuffer: (id<MTLCommandBuffer> __nonnull) commandBuffer;

@end    // AFEKMesh

#ifdef __cplusplus
}
#endif

#endif /* AFEKMesh_h */
