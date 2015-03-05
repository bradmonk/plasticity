from libc.stdio cimport fprintf
from libc.stdio cimport FILE
from libc.stdio cimport stderr


cdef public void foo_stream(FILE *stream):
    fprintf(stderr, "STREAM: EMILIO %d\n", 0);


cdef public void foo_direct():
    fprintf(stderr, "DIRECT: EMILIO %d\n", 0);
