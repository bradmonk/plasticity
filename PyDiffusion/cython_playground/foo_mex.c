#include <stdio.h>
#include <Python.h>
#include "foo.h"
#include "mex.h"

void call_py(double *zap, size_t zap_size, size_t zap_rows)
{
    Py_Initialize();
    initfoo();
    foo_stream(stderr);
    foo_direct();
    from_bar();

    int i;
    for (i = 0; i < zap_size; i++) {
      printf("zap[%d] before: %f\n", i, zap[i]);
    }
    from_waldo(zap, zap_rows);
    for (i = 0; i < zap_size; i++) {
      printf("zap[%d] after: %f\n", i, zap[i]);
    }
    Py_Finalize();
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *zap;
    size_t zap_size;
    size_t zap_rows;

    if (nrhs != 2) {
        mexErrMsgIdAndTxt("foo_mex:nrhs", "Two inputs required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("foo_mex:nlhs", "One output required.");
    }

    if (mxGetM(prhs[0]) != 1) {
        mexErrMsgIdAndTxt("foo_mex:row_vec",
                          "First input expected to be row vector.");
    }
    zap_size = mxGetN(prhs[0]);
    /* TODO: Is there a better way to get a size_t scalar? */
    zap_rows = mxGetScalar(prhs[1]);
    zap = mxGetPr(prhs[0]);
    /* We are just going to return zap, but it will be modified. */
    plhs[0] = mxCreateDoubleMatrix(1, (mwSize)zap_size, mxREAL);
    mxSetPr(plhs[0], zap);

    call_py(zap, zap_size, zap_rows);
}
