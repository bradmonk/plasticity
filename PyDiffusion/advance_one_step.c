/*==========================================================
 *
 * NOTE: This example has been cribbed from the MATLAB MEX
 *       documentation (arrayProduct.c). For now, this is a shell
 *       of functionality to ensure that the Cython generated file
 *       can be compiled with MEX.
 *
 * This is a MEX-file for MATLAB.
 *
 * The calling syntax is:
 *
 *     outMatrix = arrayProduct(multiplier, inMatrix, xyz_loc, ...
 *                              face_indices, k, initial_point, ...
 *                              initial_face_index, all_vertices, ...
 *                              triangles, face_local_bases, ...
 *                              neighbor_faces)
 *
 *========================================================
 */

#include <Python.h>
#include "_cython_interface.h"
#include "mex.h"

extern void advance_one_step_c(
    int num_points, int num_vertices, int num_triangles, double *xyz_loc,
    long *face_indices, double k, double *initial_point,
    long initial_face_index, double *all_vertices, long *triangles,
    double *face_local_bases, long *neighbor_faces);

/* The computational routine */
void arrayProduct(double x, double *y, double *z, mwSize n)
{
    mwSize i;
    /* multiply each element y by x */
    for (i=0; i<n; i++) {
        z[i] = x * y[i];
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double multiplier;              /* input scalar */
    double *inMatrix;               /* 1xN input matrix */
    size_t ncols;                   /* size of matrix */
    double *outMatrix;              /* output matrix */

    size_t num_points;              /* size of xyz_loc */
    size_t num_vertices;            /* size of all_vertices */
    size_t num_triangles;           /* size of triangles */

    double *xyz_loc;                /* 2. Mx3 input matrix */
    long *face_indices;             /* 3. Mx1 input matrix */
    double k;                       /* 4. input scalar */
    double *initial_point;          /* 5. 1x3 input matrix */
    long initial_face_index;        /* 6. input scalar */
    double *all_vertices;           /* 7. Vx3 input matrix */
    long *triangles;                /* 8. Tx3 input matrix */
    double *face_local_bases;       /* 9. Tx6 input matrix */
    long *neighbor_faces;           /* 10. Tx3 input matrix */

    /* check for proper number of arguments */
    if(nrhs!=11) {
        mexErrMsgIdAndTxt("advance_one_step:nrhs","Eleven inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("advance_one_step:nlhs","One output required.");
    }

    /* make sure the first input argument is scalar */
    if( !mxIsDouble(prhs[0]) ||
         mxIsComplex(prhs[0]) ||
         mxGetNumberOfElements(prhs[0])!=1 ) {
        mexErrMsgIdAndTxt("advance_one_step:notScalar","Input multiplier must be a scalar.");
    }

    /* make sure the second input argument is type double */
    if( !mxIsDouble(prhs[1]) ||
         mxIsComplex(prhs[1])) {
        mexErrMsgIdAndTxt("advance_one_step:notDouble","Input matrix must be type double.");
    }

    /* check that number of rows in second input argument is 1 */
    if(mxGetM(prhs[1])!=1) {
        mexErrMsgIdAndTxt("advance_one_step:notRowVector","Input must be a row vector.");
    }

    /* make sure xyz_loc is double array (of correct shape) */
    if(mxGetN(prhs[2])!=3) {
        mexErrMsgIdAndTxt("advance_one_step:notXYZCoords",
                          "xyz_loc must have 3 columns.");
    }
    if( !mxIsDouble(prhs[2]) ||
         mxIsComplex(prhs[2])) {
        mexErrMsgIdAndTxt("advance_one_step:notDouble",
                          "xyz_loc must be type double.");
    }

    /* make sure face_indices is long array (of correct shape) */
    if(mxGetN(prhs[3])!=1) {
        mexErrMsgIdAndTxt("advance_one_step:notColVector",
                          "face_indices must have 1 column.");
    }
    /* NOTE: We don't check that rows(face_indices) == rows(xyz_loc). */
    if( !mxIsInt64(prhs[3]) ) {
        mexErrMsgIdAndTxt("advance_one_step:notLong",
                          "face_indices must be type long.");
    }

    /* make sure k is scalar double */
    if( !mxIsDouble(prhs[4]) ||
         mxIsComplex(prhs[4]) ||
         mxGetNumberOfElements(prhs[4])!=1 ) {
        mexErrMsgIdAndTxt("advance_one_step:notScalarDouble",
                          "k must be a scalar double.");
    }

    /* make sure initial_point is double 1x3 array */
    if(mxGetN(prhs[5])!=3 || mxGetM(prhs[5])!=1) {
        mexErrMsgIdAndTxt("advance_one_step:notXYZVec",
                          "initial_point must be 1x3.");
    }
    if( !mxIsDouble(prhs[5]) ||
         mxIsComplex(prhs[5])) {
        mexErrMsgIdAndTxt("advance_one_step:notDouble",
                          "initial_point must be type double.");
    }

    /* make sure initial_face_index is a scalar long */
    if( !mxIsInt64(prhs[6]) ||
         mxGetNumberOfElements(prhs[6])!=1 ) {
        mexErrMsgIdAndTxt("advance_one_step:notScalar",
                          "initial_face_index must be a scalar long.");
    }

    /* make sure all_vertices is double array (of correct shape) */
    if(mxGetN(prhs[7])!=3) {
        mexErrMsgIdAndTxt("advance_one_step:notXYZCoords",
                          "all_vertices must have 3 columns.");
    }
    if( !mxIsDouble(prhs[7]) ||
         mxIsComplex(prhs[7])) {
        mexErrMsgIdAndTxt("advance_one_step:notDouble",
                          "all_vertices must be type double.");
    }

    /* make sure triangles is long array (of correct shape) */
    if(mxGetN(prhs[8])!=3) {
        mexErrMsgIdAndTxt("advance_one_step:notVertexIndices",
                          "triangles must have 3 columns.");
    }
    if( !mxIsInt64(prhs[8]) ) {
        mexErrMsgIdAndTxt("advance_one_step:notLong",
                          "triangles must be type long.");
    }

    /* make sure face_local_bases is double array (of correct shape) */
    if(mxGetN(prhs[9])!=6) {
        mexErrMsgIdAndTxt("advance_one_step:notFaceBases",
                          "face_local_bases must have 6 columns.");
    }
    /* NOTE: We don't check that rows(face_local_bases) == rows(triangles). */
    if( !mxIsDouble(prhs[9]) ||
         mxIsComplex(prhs[9])) {
        mexErrMsgIdAndTxt("advance_one_step:notDouble",
                          "face_local_bases must be type double.");
    }
    /* make sure neighbor_faces is long array (of correct shape) */
    if(mxGetN(prhs[10])!=3) {
        mexErrMsgIdAndTxt("advance_one_step:notNeighborIndices",
                          "neighbor_faces must have 3 columns.");
    }
    /* NOTE: We don't check that rows(neighbor_faces) == rows(triangles). */
    if( !mxIsInt64(prhs[10]) ) {
        mexErrMsgIdAndTxt("advance_one_step:notLong",
                          "neighbor_faces must be type long.");
    }

    /* get the value of the scalar input  */
    multiplier = mxGetScalar(prhs[0]);

    /* create a pointer to the real data in the input matrix  */
    inMatrix = mxGetPr(prhs[1]);

    /* get dimensions of the input matrix */
    ncols = mxGetN(prhs[1]);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)ncols,mxREAL);

    /* get a pointer to the real data in the output matrix */
    outMatrix = mxGetPr(plhs[0]);

    /* call the computational routine */
    arrayProduct(multiplier,inMatrix,outMatrix,(mwSize)ncols);
}
