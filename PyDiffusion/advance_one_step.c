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
 *     [xyz_loc, face_indices] = advance_one_step(...
 *         xyz_loc, face_indices, k, initial_point, ...
 *         initial_face_index, all_vertices, ...
 *         triangles, face_local_bases, neighbor_faces)
 *
 *========================================================
 */

#include <dlfcn.h>
#include <Python.h>
#include <stdio.h>
#include "_cython_interface.h"
#include "mex.h"

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    size_t num_points;              /* size of xyz_loc */
    size_t num_vertices;            /* size of all_vertices */
    size_t num_triangles;           /* size of triangles */

    double *xyz_loc;                /* 0. Mx3 input/output matrix */
    long *face_indices;             /* 1. Mx1 input/output matrix */
    double k;                       /* 2. input scalar */
    double *initial_point;          /* 3. 1x3 input matrix */
    long initial_face_index;        /* 4. input scalar */
    double *all_vertices;           /* 5. Vx3 input matrix */
    long *triangles;                /* 6. Tx3 input matrix */
    double *face_local_bases;       /* 7. Tx6 input matrix */
    long *neighbor_faces;           /* 8. Tx3 input matrix */

    /* check for proper number of arguments */
    if (nrhs != 9) {
        mexErrMsgIdAndTxt("advance_one_step:nrhs", "Nine inputs required.");
    }
    if (nlhs != 2) {
        mexErrMsgIdAndTxt("advance_one_step:nlhs", "Two outputs required.");
    }

    /* make sure xyz_loc is double array (of correct shape) */
    if (mxGetN(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("advance_one_step:notXYZCoords",
                          "xyz_loc (arg 1) must have 3 columns.");
    }
    if( !mxIsDouble(prhs[0]) ||
         mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("advance_one_step:notDouble",
                          "xyz_loc (arg 1) must be type double.");
    }

    /* make sure face_indices is long array (of correct shape) */
    if(mxGetN(prhs[1])!=1) {
        mexErrMsgIdAndTxt("advance_one_step:notColVector",
                          "face_indices (arg 2) must have 1 column.");
    }
    if( !mxIsInt64(prhs[1]) ) {
        mexErrMsgIdAndTxt("advance_one_step:notLong",
                          "face_indices (arg 2) must be type long.");
    }

    /* make sure k is scalar double */
    if( !mxIsDouble(prhs[2]) ||
         mxIsComplex(prhs[2]) ||
         mxGetNumberOfElements(prhs[2])!=1 ) {
        mexErrMsgIdAndTxt("advance_one_step:notScalarDouble",
                          "k (arg 3) must be a scalar double.");
    }

    /* make sure initial_point is double 1x3 array */
    if(mxGetN(prhs[3])!=3 || mxGetM(prhs[3])!=1) {
        mexErrMsgIdAndTxt("advance_one_step:notXYZVec",
                          "initial_point (arg 4) must be 1x3.");
    }
    if( !mxIsDouble(prhs[3]) ||
         mxIsComplex(prhs[3])) {
        mexErrMsgIdAndTxt("advance_one_step:notDouble",
                          "initial_point (arg 4) must be type double.");
    }

    /* make sure initial_face_index is a scalar long */
    if( !mxIsInt64(prhs[4]) ||
         mxGetNumberOfElements(prhs[4])!=1 ) {
        mexErrMsgIdAndTxt("advance_one_step:notScalar",
                          "initial_face_index (arg 5) must be a scalar long.");
    }

    /* make sure all_vertices is double array (of correct shape) */
    if(mxGetN(prhs[5])!=3) {
        mexErrMsgIdAndTxt("advance_one_step:notXYZCoords",
                          "all_vertices (arg 6) must have 3 columns.");
    }
    if( !mxIsDouble(prhs[5]) ||
         mxIsComplex(prhs[5])) {
        mexErrMsgIdAndTxt("advance_one_step:notDouble",
                          "all_vertices (arg 6) must be type double.");
    }

    /* make sure triangles is long array (of correct shape) */
    if(mxGetN(prhs[6])!=3) {
        mexErrMsgIdAndTxt("advance_one_step:notVertexIndices",
                          "triangles (arg 7) must have 3 columns.");
    }
    if( !mxIsInt64(prhs[6]) ) {
        mexErrMsgIdAndTxt("advance_one_step:notLong",
                          "triangles (arg 7) must be type long.");
    }

    /* make sure face_local_bases is double array (of correct shape) */
    if(mxGetN(prhs[7])!=6) {
        mexErrMsgIdAndTxt("advance_one_step:notFaceBases",
                          "face_local_bases (arg 8) must have 6 columns.");
    }
    if( !mxIsDouble(prhs[7]) ||
         mxIsComplex(prhs[7])) {
        mexErrMsgIdAndTxt("advance_one_step:notDouble",
                          "face_local_bases (arg 8) must be type double.");
    }

    /* make sure neighbor_faces is long array (of correct shape) */
    if(mxGetN(prhs[8])!=3) {
        mexErrMsgIdAndTxt("advance_one_step:notNeighborIndices",
                          "neighbor_faces (arg 9) must have 3 columns.");
    }
    if( !mxIsInt64(prhs[8]) ) {
        mexErrMsgIdAndTxt("advance_one_step:notLong",
                          "neighbor_faces (arg 9) must be type long.");
    }

    /* get dimensions of the inputs */
    num_points = mxGetM(prhs[0]);
    num_vertices = mxGetM(prhs[5]);
    num_triangles = mxGetM(prhs[6]);

    /* Check rows that should match:
     *     rows(face_indices) == rows(xyz_loc)
     *     rows(face_local_bases) == rows(triangles)
     *     rows(neighbor_faces) == rows(triangles)
     */
    if(mxGetM(prhs[1])!=num_points) {
        mexErrMsgIdAndTxt(
            "advance_one_step:mismatchRows",
            "face_indices (arg 2) and xyz_loc (arg 0) must have same # rows.");
    }
    if(mxGetM(prhs[7])!=num_triangles) {
        mexErrMsgIdAndTxt(
            "advance_one_step:mismatchRows",
            "face_local_bases (arg 8) and triangles (arg 7) must have same # rows.");
    }
    if(mxGetM(prhs[8])!=num_triangles) {
        mexErrMsgIdAndTxt(
            "advance_one_step:mismatchRows",
            "neighbor_faces (arg 9) and triangles (arg 7) must have same # rows.");
    }

    /* get scalar values */
    k = mxGetScalar(prhs[2]);
    initial_face_index = mxGetScalar(prhs[4]);

    /* get matrix values */
    initial_point = mxGetPr(prhs[3]);
    all_vertices = mxGetPr(prhs[5]);
    triangles = (long *)mxGetData(prhs[6]);
    face_local_bases = mxGetPr(prhs[7]);
    neighbor_faces = (long *)mxGetData(prhs[8]);

    /* Copy xyz_loc and face_indices before editing. (If
     * not copied, will cause segfault.)
     */
    plhs[0] = mxDuplicateArray(prhs[0]);
    plhs[1] = mxDuplicateArray(prhs[1]);

    xyz_loc = mxGetPr(plhs[0]);
    face_indices = (long *)mxGetData(plhs[1]);

    void* handle = dlopen("libpython2.7.so", RTLD_LAZY | RTLD_GLOBAL);
    Py_Initialize();
    init_cython_interface();
    advance_one_step_matlab((mwSize) num_points, (mwSize) num_vertices,
                            (mwSize) num_triangles,
                            xyz_loc, face_indices, k, initial_point,
                            initial_face_index, all_vertices, triangles,
                            face_local_bases, neighbor_faces);
    Py_Finalize();
    dlclose(handle);
}
