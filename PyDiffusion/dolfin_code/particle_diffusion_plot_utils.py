"""Plotting partner to particle_diffusion_on_mesh."""

import matplotlib
matplotlib.use('TKAgg')
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3

from particle_diffusion_on_mesh import Point


class PlotBoundary(object):
    """Container object to define a plot area for an animation.

    Example:
        >>> # 17.0 is just slightly bigger than 15
        >>> plot_boundary = PlotBoundary(-17.0, 17.0, -17.0, 17.0, 0.0, 50.0)
    """

    def __init__(self, min_x, max_x, min_y, max_y, min_z, max_z):
        self.min_x = min_x
        self.max_x = max_x
        self.min_y = min_y
        self.max_y = max_y
        self.min_z = min_z
        self.max_z = max_z


def default_color_function(point):
    return 'b'


# Mostly borrowed from particle_diffusion.py (without using dolfin).
def plot_simulation(num_points, mesh_wrapper, plot_boundary,
                    color_function=default_color_function, show_mesh=False,
                    num_frames=200, print_frequency=None,
                    interval=30, filename=None):
    points = [Point(mesh_wrapper) for _ in xrange(num_points)]

    # Attaching 3D axis to the figure
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    # Set default view.
    ax.view_init(elev=24, azim=-144)

    # Create lines with a single point as a scatter.
    all_points = [ax.plot([pt.x], [pt.y], [pt.z], c='b', marker='o')[0]
                  for pt in points]
    if show_mesh:
        # Add permanent (unchanged) features to plot.
        x = mesh_wrapper.all_vertices[:, 0]
        y = mesh_wrapper.all_vertices[:, 1]
        z = mesh_wrapper.all_vertices[:, 2]
        ax.plot_trisurf(x, y, z, triangles=mesh_wrapper.triangles,
                        color=(0, 0, 0, 0), edgecolor='Gray', linewidth=0.05)

    def update_plot(step_num):
        if print_frequency is not None and step_num % print_frequency == 0:
            print 'Step Number:', step_num

        for pt, point_container in zip(points, all_points):
            pt.move()
            point_container.set_color(color_function(pt))
            point_container.set_data([pt.x], [pt.y])
            point_container.set_3d_properties([pt.z])

        return all_points

    # Setting the axes properties
    ax.set_xlim3d([plot_boundary.min_x, plot_boundary.max_x])
    ax.set_xlabel('X')

    ax.set_ylim3d([plot_boundary.min_y, plot_boundary.max_y])
    ax.set_ylabel('Y')

    ax.set_zlim3d([plot_boundary.min_z, plot_boundary.max_z])
    ax.set_zlabel('Z')

    create_gif = filename is not None
    anim = animation.FuncAnimation(fig, update_plot,
                                   repeat=False, frames=num_frames,
                                   interval=interval, blit=create_gif)

    if create_gif:
        anim.save(filename, writer='imagemagick_file')
    else:
        # Could also use plt.draw() for interactive work, not for 3D though.
        plt.show(block=False)
