import dolfin
import numpy as np
import os
import scipy.io

############################  SHELL ENVIRONMENT  ############################
# ENABLED:
# export PATH="$PATH:/Library/Frameworks/Python.framework/Versions/2.7/bin"
# source "/Applications/FEniCS.app/Contents/Resources/share/fenics/fenics.conf"
# export PATH="$PATH:/Applications/FEniCS.app/Contents/Resources/share/fenics"
# DISABLED:
# # source "/Users/bradleymonk/Library/Enthought/Canopy_64bit/User/bin/activate"
# # export PATH="$PATH:/Users/bradleymonk/Library/Enthought/Canopy_64bit/User/bin"
# # export PATH="$PATH:/Users/bradleymonk/anaconda/bin"
############################  SHELL ENVIRONMENT  ############################

if not dolfin.has_cgal():
  raise ImportError('Does not have CGAL.')
else:
  print('Done importing')
  print('================================')

def main():

  filename = 'data/surfmesh.xml'

  mesh_3d = dolfin.Mesh(filename)
  
  # Plot in either case.
  print('Plotting in interactive mode')
  fh = dolfin.plot(mesh_3d, '3D mesh', interactive=True,
    window_width=900, window_height=700)


if __name__ == '__main__':
  main()
