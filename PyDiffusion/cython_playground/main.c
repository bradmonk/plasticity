#include <stdio.h>
#include <Python.h>
#include "foo.h"

int main(void)
{
    Py_Initialize();
    initfoo();
    foo_stream(stderr);
    foo_direct();
    from_bar();
    double zap[3] = {1.0, 2.0, 3.0};
    fprintf(stderr, "zap[0] before: %f\n", zap[0]);
    fprintf(stderr, "zap[1] before: %f\n", zap[1]);
    fprintf(stderr, "zap[2] before: %f\n", zap[2]);
    from_waldo(&zap[0], 3);
    fprintf(stderr, "zap[0] after: %f\n", zap[0]);
    fprintf(stderr, "zap[1] after: %f\n", zap[1]);
    fprintf(stderr, "zap[2] after: %f\n", zap[2]);
    Py_Finalize();
    return 0;
}
