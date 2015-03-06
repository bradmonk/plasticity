import numpy as np

from cython cimport view
from libc.stdio cimport fprintf
from libc.stdio cimport FILE
from libc.stdio cimport stderr

from bar import waldo
from bar import qux


cdef public void foo_stream(FILE *stream):
    fprintf(stderr, "  STREAM: EMILIO %d\n", 0);


cdef public void foo_direct():
    fprintf(stderr, "  DIRECT: EMILIO %d\n", 0);


cdef public void from_bar():
    cdef double x = qux()
    fprintf(stderr, "FROM BAR: EMILIO %f\n", x);


cdef public void from_waldo(double* zap, size_t zap_rows):
    cdef view.array py_zap = <double[:zap_rows, :2]> zap
    waldo(np.asarray(py_zap))
