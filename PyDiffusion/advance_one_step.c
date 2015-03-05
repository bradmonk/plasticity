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
 *                              face_indices)
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

    /* check for proper number of arguments */
    if(nrhs!=4) {
        mexErrMsgIdAndTxt("advance_one_step:nrhs","Four inputs required.");
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
