# Directory Contents
----

This library contains mostly scripts, but has some utility modules
which can be re-used.

Some of the scripts generate data and some depend on data to run simulations
and/or generate images.

The generated data and images are stored in their own subdirectories:

```
DIFFUSION-RELATED UTILS FOUND IN:
    advance_one_step_mex/

MESH GENERATION UTILS FOUND IN:
    build_mesh/
```


#### Subdirectory Contents: advance_one_step_mex/ 

1. FILE: Makefile
    RUN: $make test_mex
    MAKES: advance_one_step.mexmaci64
   This file is necessary to generate the particle diffusion steps in...
1. FILE: generate_diffusion_paths.py
    Uses advance_one_step.mexmaci64 (no `dolfin` required) to generate paths
    (serialize mesh file can be used directly in python to generate paths)


#### Subdirectory Contents: build_mesh/
1. FILE: generate_mesh.py
    Creates mesh for points to travel on (`dolfin`) and save as XML file
1. FILE: serialize_mesh.py
    Converts full mesh to serialized mesh.
    Saves only the necessary surface geometry from Dolfin-generated mesh


----
# Reusable Modules
----

```
_cython_interface.pyx
dolfin_mesh_utils.py
particle_diffusion_on_mesh.py
particle_diffusion_plot_utils.py
```

In particular `_cython_interface.pyx` generates a Python importable shared
object file `_cython_interface.so` when `make cython_interface` is called.
This shared object file can / will also be used in a MEX file to call
the Python code from MATLAB.


#### Module Notes

The scripts can be divided into two purposes:

1.  Play around creating / refining meshes

    ```
    dendritic_shaft_mesh.py
    face_neighbors_and_top_components.py
    find_near_top_faces.py
    find_synapse_top_and_exterior_faces.py
    full_dendrite_mesh.py
    plot_trisurf_mesh.py
    refine_boundary_mesh.py
    ```

1.  Play around simulating particle movement on a given mesh

    ```
    sample_code.py
    ```
----
## Issues
----

The `dolfin` library depends on `ScientificPython` (a very old library that
built off of `numpy` **long ago**). This dependency is broken by
`numpy>=1.9`.

Danny [suggested][1] (and has [deployed][2]) a hack to work-around this
versioning issue.

**CGAL** is no longer supported in `dolfin` and many of the old scripts
stopped working. A [quote][3] (mailing list) from November 2014:

> Yes, CGAL support was disabled in favor of `mshr`:
>
> https://bitbucket.org/benjamik/mshr/wiki/Home

[1]: https://bitbucket.org/khinsen/scientificpython/issue/13/
[2]: https://gist.github.com/dhermes/38d8ff05267e861a4b01
[3]: http://fenicsproject.org/pipermail/fenics-support/2014-November/000961.html
