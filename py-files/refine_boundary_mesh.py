import dolfin


def main():
  resolution = 96  # 32 * 3
  filename = 'mesh_res_%d_boundary.xml' % resolution

  mesh_3d = dolfin.Mesh(filename)
  mesh_3d = dolfin.mesh.refine(mesh_3d)

  print '=' * 60
  print 'Plotting in interactive mode'
  print '=' * 60
  dolfin.plot(mesh_3d, '3D mesh', interactive=True)


if __name__ == '__main__':
  main()
