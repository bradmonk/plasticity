#include <dlfcn.h>
#include <Python.h>
#include <stdio.h>
#include "foo.h"
#include "mex.h"

void call_py(double *zap, size_t zap_size, size_t zap_rows)
{
    void* handle = dlopen("libpython2.7.so", RTLD_LAZY | RTLD_GLOBAL);
    Py_Initialize();
    initfoo();
    from_waldo(zap, zap_rows);
    Py_Finalize();
    dlclose(handle);
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *out_matrix;
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

    /* Create a copy of input to be modified. */
    plhs[0] = mxDuplicateArray(prhs[0]);
    out_matrix = mxGetPr(plhs[0]);
    call_py(out_matrix, zap_size, zap_rows);
}
