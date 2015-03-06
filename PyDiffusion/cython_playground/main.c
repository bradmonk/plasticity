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
    double zap[5] = {1.0, 2.0, 3.0, 4.0, 5.0};
    size_t zap_size = 5;

    int i;
    for (i = 0; i < zap_size; i++) {
      printf("zap[%d] before: %f\n", i, zap[i]);
    }
    from_waldo(&zap[0], zap_size - 2);
    for (i = 0; i < zap_size; i++) {
      printf("zap[%d] after: %f\n", i, zap[i]);
    }
    Py_Finalize();
    return 0;
}
