#include <stdio.h>
#include <Python.h>
#include "foo.h"

int main(void)
{
    Py_Initialize();
    initfoo();
    double zap[7] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
    size_t zap_size = 7;
    size_t zap_rows = 3; // 3x2 is not quite 7

    int i;
    for (i = 0; i < zap_size; i++) {
        printf("zap[%d] before: %f\n", i, zap[i]);
    }
    from_waldo(&zap[0], zap_rows);
    puts("================================================================");
    for (i = 0; i < zap_size; i++) {
        printf("zap[%d] after: %f\n", i, zap[i]);
    }
    Py_Finalize();
    return 0;
}
