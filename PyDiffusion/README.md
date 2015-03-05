## Pipeline for Running Diffusion Simulation

1. Create a mesh for the points to travel on (`dolfin`)
   and save to an XML file
1. Convert the Dolfin-generated mesh into just the parts we
   need for the surface geometry (e.g. `sample_code.save_serialized_mesh()`
1. Use the serialized mesh (no `dolfin` required) in
   `particle_diffusion_on_mesh.run_simulation` (or via
   `Mesh.from_file` directly)

### Directory Contents

This library contains mostly scripts, but has some utility modules
which can be re-used.

Some of the scripts generate data and some depend on data to run simulations
and/or generate images.

The generated data and images are stored in their own subdirectories:

```
data/
plots/
```

The re-usable modules are

```
dolfin_mesh_utils.py
particle_diffusion_on_mesh.py
particle_diffusion_plot_utils.py
```

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

### Issues

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
