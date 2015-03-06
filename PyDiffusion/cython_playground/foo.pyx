import numpy as np

from cython cimport view
from libc.stdio cimport fprintf
from libc.stdio cimport FILE
from libc.stdio cimport stderr

from bar import waldo
from bar import qux


cdef public void foo_stream(FILE *stream):
    fprintf(stderr, "STREAM: EMILIO %d\n", 0);


cdef public void foo_direct():
    fprintf(stderr, "DIRECT: EMILIO %d\n", 0);


cdef public void from_bar():
    cdef double x = qux()
    fprintf(stderr, "FRMBAR: EMILIO %f\n", x);


cdef public void from_waldo(double* zap, size_t zap_size):
    cdef view.array py_zap = <double[:zap_size]> zap
    waldo(np.asarray(py_zap))
