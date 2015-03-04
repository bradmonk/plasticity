## Pipeline for Running Diffusion Simulation

1. Create a mesh for the points to travel on (`dolfin`)
   and save to an XML file
1. Convert the Dolfin-generated mesh into just the parts we
   need for the surface geometry (e.g. `sample_code.save_serialized_mesh()`
1. Use the serialized mesh (no `dolfin` required) in
   `particle_diffusion_on_mesh.run_simulation` (or via
   `Mesh.from_file` directly)
