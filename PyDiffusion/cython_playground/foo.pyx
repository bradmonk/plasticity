from libc.stdio cimport fprintf
from libc.stdio cimport FILE
from libc.stdio cimport stderr

cdef public void foo_stream(FILE *stream):
    fprintf(stream, "HEY EMILIO %d\n", 0);
    fprintf(stream, "HEY EMILIO %d\n", 1);
    fprintf(stream, "HEY EMILIO %d\n", 100);


cdef public void foo_direct():
    fprintf(stderr, "HEY EMILIO %d\n", 0);
    fprintf(stderr, "HEY EMILIO %d\n", 1);
    fprintf(stderr, "HEY EMILIO %d\n", 100);
