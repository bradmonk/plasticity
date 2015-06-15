#ifndef __PYX_HAVE___cython_interface
#define __PYX_HAVE___cython_interface


#ifndef __PYX_HAVE_API___cython_interface

#ifndef __PYX_EXTERN_C
  #ifdef __cplusplus
    #define __PYX_EXTERN_C extern "C"
  #else
    #define __PYX_EXTERN_C extern
  #endif
#endif

#ifndef DL_IMPORT
  #define DL_IMPORT(_T) _T
#endif

__PYX_EXTERN_C DL_IMPORT(void) advance_one_step_c(size_t, size_t, size_t, double *, long *, double, double *, long, double *, long *, double *, long *);

#endif /* !__PYX_HAVE_API___cython_interface */

#if PY_MAJOR_VERSION < 3
PyMODINIT_FUNC init_cython_interface(void);
#else
PyMODINIT_FUNC PyInit__cython_interface(void);
#endif

#endif /* !__PYX_HAVE___cython_interface */
