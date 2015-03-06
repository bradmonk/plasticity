import numpy as np

from cython cimport view

from bar import waldo


cdef public void from_waldo(double* zap, size_t zap_rows):
    cdef view.array py_zap = <double[:zap_rows, :2]> zap
    waldo(np.asarray(py_zap))
