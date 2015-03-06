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
    Py_Finalize();
    return 0;
}
