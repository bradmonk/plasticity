import numpy as np


class PointType(object):
  SPHERE = 1
  VERT_DENDRITE = 2
  HORIZ_DENDRITE = 3


class Point(object):

  SPHERE_Z_CENTER = 0.7
  SPHERE_RADIUS = 0.1
  SPHERE_RADIUS_SQUARED = 0.1**2
  SPHERE_CHANGE_POINT = 0.7 - np.sqrt(0.0075)

  VERT_DENDRITE_RADIUS = 0.05
  VERT_DENDRITE_RADIUS_SQUARED = 0.05**2
  DENDRITE_CHANGE_TOP = 0.5
  DENDRITE_CHANGE_BOTTOM = np.sqrt(0.2475)
  HORIZ_DENDRITE_RADIUS = 0.5
  HORIZ_DENDRITE_RADIUS_SQUARED = 0.5**2

  BEGIN_MIXED_HORIZ_THETA = np.arctan2(np.sqrt(0.2475), 0.05)
  END_MIXED_HORIZ_THETA = np.arctan2(np.sqrt(0.2475), -0.05)

  def __init__(self, x, y, z, k):
    self.x = x
    self.y = y
    self.z = z
    self.k = k

    self.classify()

  def classify(self):
    if self.z >= self.SPHERE_CHANGE_POINT:
      self.point_type = PointType.SPHERE
      self.verify_sphere()
    elif self.z >= self.DENDRITE_CHANGE_TOP:
      self.point_type = PointType.VERT_DENDRITE
      self.verify_vert_dendrite()
    elif self.z <= self.DENDRITE_CHANGE_BOTTOM:
      self.point_type = PointType.HORIZ_DENDRITE
      self.verify_horiz_dendrite()
    else:
      self.classify_between_dendrites()

  def verify_sphere(self):
    """Verifies that a point is on the sphere.

    Assumes the calling code has already checked the z-value.
    """
    effective_radius_squared = (self.x**2 + self.y**2 +
                                (self.z - self.SPHERE_Z_CENTER)**2)
    if not np.allclose(self.SPHERE_RADIUS_SQUARED,
                       effective_radius_squared):
      raise ValueError('Not on sphere.')

    # self.top_point = None
    # self.right_point = None

  def verify_vert_dendrite(self):
    """Verifies that a point is on the vertical dendrite.

    NOTE: This will be changed as we expand from one vertical dendrite
          to many and in that case will likely store more data on the object.

    Assumes the calling code has already checked the z-value.
    """
    effective_radius_squared = self.x**2 + self.y**2
    if not np.allclose(self.VERT_DENDRITE_RADIUS_SQUARED,
                       effective_radius_squared):
      raise ValueError('Not on vertical dendrite.')

    self.theta = np.arctan2(self.y, self.x)

  def verify_horiz_dendrite(self):
    """Verifies that a point is on the horizontal dendrite.

    Assumes the calling code has already checked the z-value.
    """
    effective_radius_squared = self.y**2 + self.z**2
    if not np.allclose(self.HORIZ_DENDRITE_RADIUS_SQUARED,
                       effective_radius_squared):
      raise ValueError('Not on horizontal dendrite.')

    self.theta = np.arctan2(self.z, self.y)

  def classify_between_dendrites(self):
    """Determines which dendrite a point lies on.

    Assumes the calling code has already checked the z-value.
    """
    if self.x**2 + self.y**2 > self.VERT_DENDRITE_RADIUS_SQUARED:
      self.point_type = PointType.HORIZ_DENDRITE
      self.verify_horiz_dendrite()
    else:
      self.point_type = PointType.VERT_DENDRITE
      self.verify_vert_dendrite()

  def __str__(self):
    type_str = 'Sphere'
    if self.point_type == PointType.VERT_DENDRITE:
      type_str = 'VerticalDendrite'
    elif self.point_type == PointType.HORIZ_DENDRITE:
      type_str = 'HorizontalDendrite'
    return '%s(x=%2.2f, y=%2.2f, z=%2.2f), k=%2.2f' % (
        type_str, self.x, self.y, self.z, self.k)

  def __repr__(self):
    return self.__str__()

  def move(self):
    """Moves the point.

    Assumes the point type is correctly set.
    """
    if self.point_type == PointType.SPHERE:
      self.move_on_sphere()
    elif self.point_type == PointType.VERT_DENDRITE:
      self.move_on_vert_dendrite()
    else:  # Assumes HORIZ_DENDRITE
      self.move_on_horiz_dendrite()

  def move_on_sphere(self, x_rand=None, y_rand=None):
    x_rand = x_rand or self.k * np.random.randn()
    y_rand = y_rand or self.k * np.random.randn()

  def move_on_vert_dendrite(self, x_rand=None, y_rand=None):
    x_rand = x_rand or self.k * np.random.randn()
    y_rand = y_rand or self.k * np.random.randn()

    new_z = self.z + y_rand
    if new_z > self.SPHERE_CHANGE_POINT:
      # self.move_from_vert_dendrite_to_sphere(x_rand, y_rand)
      print 'move_from_vert_dendrite_to_sphere'
      return
    elif new_z < self.DENDRITE_CHANGE_TOP:
      # self.move_from_vert_dendrite_to_horiz_dendrite(x_rand, y_rand)
      print 'move_from_vert_dendrite_to_horiz_dendrite'
      return

    # If our z-value doesn't exceed the boundaries, we can do simply math.
    self.z = new_z
    # (L / (2 pi R)) * (2 pi)
    theta_delta = x_rand / self.VERT_DENDRITE_RADIUS
    if np.abs(theta_delta) > np.pi:
      raise ValueError('Extremely large rotation. Consider adjusting k.')
    self.theta += theta_delta
    self.x = self.VERT_DENDRITE_RADIUS * np.cos(self.theta)
    self.y = self.VERT_DENDRITE_RADIUS * np.sin(self.theta)

  def move_on_horiz_dendrite(self, x_rand=None, y_rand=None):
    x_rand = x_rand or self.k * np.random.randn()
    y_rand = y_rand or self.k * np.random.randn()

    new_x = self.x + x_rand
    # If we started left and stayed left or started right and stayed
    # right, we are safe.
    if ((self.x <= -self.VERT_DENDRITE_RADIUS and
        new_x <= -self.VERT_DENDRITE_RADIUS) or
        (self.x >= self.VERT_DENDRITE_RADIUS and
         new_x >= self.VERT_DENDRITE_RADIUS)):
      self.x = new_x

      # (L / (2 pi R)) * (2 pi)
      theta_delta = y_rand / self.HORIZ_DENDRITE_RADIUS
      if np.abs(theta_delta) > np.pi:
        raise ValueError('Extremely large rotation. Consider adjusting k.')
      self.theta += theta_delta
      self.y = self.HORIZ_DENDRITE_RADIUS * np.cos(self.theta)
      self.z = self.HORIZ_DENDRITE_RADIUS * np.sin(self.theta)

      # We are finished updating, so we return.
      return

    # In this case, our x-values either started in the same range as the
    # vertical dendrite or have crossed those x-values.
    if self.x <= -self.VERT_DENDRITE_RADIUS:  # Started left
      print 'Started left'
    elif self.x >= self.VERT_DENDRITE_RADIUS:  # Started right
      print 'Started right'
    else:  # Started in range of cylinder.
      print 'Started in range of cylinder'


def test_init():
  dummy_k = 0.01
  sphere_north_pole = Point.SPHERE_Z_CENTER + Point.SPHERE_RADIUS
  p1 = Point(0, 0, sphere_north_pole, dummy_k)
  assert p1.point_type == PointType.SPHERE, 'p1'

  sphere_south_pole = Point.SPHERE_Z_CENTER - Point.SPHERE_RADIUS
  error_occurred = False
  try:
    p2 = Point(0, 0, sphere_south_pole, dummy_k)
  except ValueError:
    error_occurred = True
  assert error_occurred, 'p2'

  p3 = Point(0, Point.HORIZ_DENDRITE_RADIUS, 0, dummy_k)
  assert p3.point_type == PointType.HORIZ_DENDRITE, 'p3 type'
  assert np.allclose(p3.theta, 0), 'p3'
  p4 = Point(0, 0, -p3.y, dummy_k)
  assert p4.point_type == PointType.HORIZ_DENDRITE, 'p4 type'
  assert np.allclose(p4.theta, -np.pi / 2), 'p3'
  p5 = Point(0, -p3.y, 0, dummy_k)
  assert p5.point_type == PointType.HORIZ_DENDRITE, 'p5 type'
  assert np.allclose(p5.theta, np.pi), 'p5'
  p6 = Point(0, -np.sqrt(Point.HORIZ_DENDRITE_RADIUS_SQUARED - 0.0001**2),
             -0.0001, dummy_k)
  assert p6.point_type == PointType.HORIZ_DENDRITE, 'p6 type'
  assert 0 <= p6.theta + np.pi < 0.01, 'p6'

  p7 = Point(Point.VERT_DENDRITE_RADIUS, 0,
             Point.DENDRITE_CHANGE_TOP, dummy_k)
  assert p7.point_type == PointType.VERT_DENDRITE, 'p7'

p_v = Point(Point.VERT_DENDRITE_RADIUS, 0,
            Point.DENDRITE_CHANGE_TOP, 0.01)
p_h = Point(0, Point.HORIZ_DENDRITE_RADIUS, 0, 0.01)


mid_point_z = 0.5 * (Point.SPHERE_CHANGE_POINT + Point.DENDRITE_CHANGE_TOP)
points = [Point(Point.VERT_DENDRITE_RADIUS, 0, mid_point_z, 0.01)
          for i in xrange(50)]
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
xs = [pt.x for pt in points]
ys = [pt.y for pt in points]
zs = [pt.z for pt in points]
ax.scatter(xs, ys, zs, c='b', marker='o')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()
plt.clf()

[pt.move() for pt in points]
