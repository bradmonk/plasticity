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
