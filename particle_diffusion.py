import matplotlib
matplotlib.use('TKAgg')
import math
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import mpmath
import numpy as np
# To catch operations that produce nans.
np.seterr(invalid='raise')


DEBUG = False


class PointType(object):
  SPHERE = 1
  VERT_DENDRITE = 2
  HORIZ_DENDRITE = 3

  COLOR_MAP = {
      SPHERE: 'b',
      VERT_DENDRITE: 'g',
      HORIZ_DENDRITE: 'r',
  }


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
    # Consider also storing the values on the object.
    self.move_counter = 0

    self.classify()

  def classify(self):
    if self.z >= self.SPHERE_CHANGE_POINT:
      self.verify_sphere()
    elif self.z >= self.DENDRITE_CHANGE_TOP:
      self.verify_vert_dendrite()
    elif self.z <= self.DENDRITE_CHANGE_BOTTOM:
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

    self.point_type = PointType.SPHERE

    first_two_coords_norm = np.sqrt(self.x**2 + self.y**2)
    self.first_two_coords_norm = first_two_coords_norm
    if not np.allclose(first_two_coords_norm, 0):
      self.centered_top = np.array([
          self.x * (self.SPHERE_Z_CENTER - self.z) / first_two_coords_norm,
          self.y * (self.SPHERE_Z_CENTER - self.z) / first_two_coords_norm,
          first_two_coords_norm,
      ])
      self.centered_right = np.array([
          - self.SPHERE_RADIUS * self.y / first_two_coords_norm,
          self.SPHERE_RADIUS * self.x / first_two_coords_norm,
          0,
      ])
    else:
      # This is the north pole. We pick an arbitrary direction -- (1, 0, 0)
      # -- for the "top" and then use the right hand rule to get the other
      # direction: (1, 0, 0) x (0, 0, 1) = (0, -1, 0).
      self.centered_top = np.array([self.SPHERE_RADIUS, 0, 0])
      self.centered_right = np.array([0, -self.SPHERE_RADIUS, 0])

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

    self.point_type = PointType.VERT_DENDRITE
    self.theta = np.arctan2(self.y, self.x)

  def verify_horiz_dendrite(self):
    """Verifies that a point is on the horizontal dendrite.

    Assumes the calling code has already checked the z-value.
    """
    effective_radius_squared = self.y**2 + self.z**2
    if not np.allclose(self.HORIZ_DENDRITE_RADIUS_SQUARED,
                       effective_radius_squared):
      raise ValueError('Not on horizontal dendrite.')

    self.point_type = PointType.HORIZ_DENDRITE
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

    self.move_counter += 1

  def _move_on_sphere(self, x_rand, y_rand):
    L = np.linalg.norm([x_rand, y_rand])
    theta = np.arctan2(y_rand, x_rand)

    theta_prime = L / self.SPHERE_RADIUS
    if np.abs(theta_prime) > np.pi:
      raise ValueError('Extremely large rotation. Consider adjusting k.')
    centered_point = np.array([self.x, self.y,
                               self.z - self.SPHERE_Z_CENTER])
    new_centered_point = (
        math.cos(theta_prime) * centered_point +
        math.sin(theta_prime) * (
            math.cos(theta) * self.centered_right +
            math.sin(theta) * self.centered_top
        )
    )
    self.x = new_centered_point[0]
    self.y = new_centered_point[1]
    self.z = new_centered_point[2] + self.SPHERE_Z_CENTER

  def move_on_sphere(self, x_rand=None, y_rand=None):
    x_rand = x_rand or self.k * np.random.randn()
    y_rand = y_rand or self.k * np.random.randn()

    L = np.linalg.norm([x_rand, y_rand])
    theta = np.arctan2(y_rand, x_rand)

    # BOTTOM HALF and moving downwards
    if self.z < self.SPHERE_Z_CENTER and theta <= 0:
      x_y_contrib = ((self.VERT_DENDRITE_RADIUS /
                      self.SPHERE_RADIUS_SQUARED) *
                     self.first_two_coords_norm)
      z_contrib = (self.SPHERE_Z_CENTER - self.z) * (
          (self.SPHERE_Z_CENTER - self.SPHERE_CHANGE_POINT) /
          self.SPHERE_RADIUS_SQUARED)
      d = self.SPHERE_RADIUS * np.arccos(x_y_contrib + z_contrib)

      theta1 = theta + np.pi/2
      H = d / math.cos(theta1)
      x_rot = H * math.sin(theta1)

      if L > H:
        self.move_from_sphere_to_vert_dendrite(H, x_rot, L, theta1)
        return

    self._move_on_sphere(x_rand, y_rand)
    self.verify_sphere()

  def move_from_sphere_to_vert_dendrite(self, H, x_rot, L, theta1):
    x_lip_cyl = (self.VERT_DENDRITE_RADIUS * self.x /
                 self.first_two_coords_norm)
    y_lip_cyl = (self.VERT_DENDRITE_RADIUS * self.y /
                 self.first_two_coords_norm)
    # The z is always SPHERE_CHANGE_POINT on the lip of the cylinder.
    theta_cyl = np.arctan2(y_lip_cyl, x_lip_cyl)

    # Update the z value down from the lip of the cylinder/sphere boundary.
    self.z = self.SPHERE_CHANGE_POINT - (L - H) * math.cos(theta1)

    x_displace = (L - H) * math.sin(theta1) + x_rot
    theta_delta = x_displace / self.VERT_DENDRITE_RADIUS
    if np.abs(theta_delta) > np.pi:
      raise ValueError('Extremely large rotation. Consider adjusting k.')

    self.theta = theta_cyl + theta_delta
    self.x = self.VERT_DENDRITE_RADIUS * math.cos(self.theta)
    self.y = self.VERT_DENDRITE_RADIUS * math.sin(self.theta)
    self.verify_vert_dendrite()

  def _move_on_vert_dendrite(self, x_rand, y_rand):
    self.z = self.z + y_rand
    # (L / (2 pi R)) * (2 pi)
    theta_delta = x_rand / self.VERT_DENDRITE_RADIUS
    if np.abs(theta_delta) > np.pi:
      raise ValueError('Extremely large rotation. Consider adjusting k.')
    self.theta += theta_delta
    self.x = self.VERT_DENDRITE_RADIUS * math.cos(self.theta)
    self.y = self.VERT_DENDRITE_RADIUS * math.sin(self.theta)

  def move_on_vert_dendrite(self, x_rand=None, y_rand=None):
    x_rand = x_rand or self.k * np.random.randn()
    y_rand = y_rand or self.k * np.random.randn()

    new_z = self.z + y_rand
    if new_z > self.SPHERE_CHANGE_POINT and y_rand > 0:
      # y_rand > 0 is redundant here since self.z <= SPHERE_CHANGE_POINT.
      self.move_from_vert_dendrite_to_sphere(x_rand, y_rand)
      return
    elif new_z < self.DENDRITE_CHANGE_TOP and y_rand < 0:
      # This is only a problem if we are moving down.
      self.move_from_vert_dendrite_to_horiz_dendrite(x_rand, y_rand)
      return

    # If our z-value doesn't exceed the boundaries, we can do simply math.
    self._move_on_vert_dendrite(x_rand, y_rand)
    # Since new_z stayed on the vertical dendrite, we don't need to change
    # the point_type.
    self.verify_vert_dendrite()

  def move_from_vert_dendrite_to_sphere(self, x_rand, y_rand):
    """Moves a point from the vertical dendrite to the sphere.

    Assumes new_z > SPHERE_CHANGE_POINT, this forces y_rand to be
    positive.
    """
    height_change_in_dendrite = self.SPHERE_CHANGE_POINT - self.z
    rot_displ_in_dendrite = (x_rand / y_rand) * height_change_in_dendrite
    self._move_on_vert_dendrite(rot_displ_in_dendrite,
                                height_change_in_dendrite)
    # We need to make sure the point is now on the sphere and set the
    # sphere-relevant data (e.g. the top and right points).
    self.verify_sphere()

    height_change_in_sphere = y_rand - height_change_in_dendrite
    rot_displ_in_sphere = x_rand - rot_displ_in_dendrite
    self._move_on_sphere(rot_displ_in_sphere,
                         height_change_in_sphere)
    self.verify_sphere()

  def move_from_vert_dendrite_to_horiz_dendrite(self, x_rand, y_rand):
    """Attempts to determine move at boundary of dendrites in vertical.

    May not actually leave the vertical dendrite, but has been identified
    as a potential to do so.

    Since we are moving down, we assume that y_rand < 0.
    """
    z0 = self.z
    theta0 = self.theta

    def curved_fn(theta):
      return mpmath.sqrt(
          self.HORIZ_DENDRITE_RADIUS_SQUARED -
          self.VERT_DENDRITE_RADIUS_SQUARED * mpmath.sin(theta)**2)

    change_regions = False
    intersection_x_rand = x_rand
    intersection_y_rand = y_rand
    m = None
    if x_rand == 0:
      z_intersect = float(curved_fn(theta0))
      # If the intersection point is above the new z, then we must change
      # regions.
      if z_intersect > self.z + y_rand:
        change_regions = True
        # We assume that this is negative, but only a portion of y_rand.
        intersection_y_rand = z_intersect - self.z
        # intersection_x_rand is already 0
    else:
      m = y_rand / x_rand

      def line_fn(theta):
        return z0 + m * (theta - theta0)

      def intersect_fn(theta):
        return curved_fn(theta) - line_fn(theta)

      intersection_theta = float(mpmath.findroot(intersect_fn, theta0))

      endpoints = sorted((theta0, theta0 + x_rand))
      if endpoints[0] <= intersection_theta <= endpoints[1]:
        change_regions = True
        z_intersect = float(curved_fn(intersection_theta))
        intersection_y_rand = z_intersect - self.z
        # We assume y_rand < 0, hence != 0
        intersection_fraction = intersection_y_rand / y_rand
        # This should be equal in absolute value to the distance from
        # theta0 and intersection_theta.
        intersection_x_rand = x_rand * intersection_fraction

    # We move up to the point of intersection (if it occurs) first since we
    # are still on the vertical dendrite.
    self._move_on_vert_dendrite(intersection_x_rand, intersection_y_rand)

    if not change_regions:
      # NOTE: This may cause problems due to ties.
      self.verify_vert_dendrite()
    else:
      # Verify before moving within the new region so all the correct
      # values are set.
      try:
        self.verify_horiz_dendrite()
      except ValueError:
        if DEBUG:
          print 'Not on horizontal dendrite in exact sense.'
          old_z = self.z
          print 'Approximating by updating z by',
        # Failed test to be on HorizontalDendrite, so force the
        # z-value to make this occur.
        z_sq = self.HORIZ_DENDRITE_RADIUS_SQUARED - self.y**2
        self.z = np.sqrt(z_sq)
        if DEBUG:
          print (self.z - old_z)
        self.verify_horiz_dendrite()

      remaining_x_rand = x_rand - intersection_x_rand
      remaining_y_rand = y_rand - intersection_y_rand
      remaining_length = np.sqrt(remaining_x_rand**2 + remaining_y_rand**2)

      theta_perp = np.arctan2(self.x, self.y)
      if m is None:
        theta_s = np.arctan2(1.0, 0.0)
      else:
        theta_s = np.arctan2(np.abs(m), np.sign(m))

      theta_new_dir = theta_perp + theta_s - np.pi / 2

      new_x_rand = remaining_length * math.cos(theta_new_dir)
      new_y_rand = remaining_length * math.sin(theta_new_dir)

      self._move_on_horiz_dendrite(new_x_rand, new_y_rand)
      self.verify_horiz_dendrite()

  def _move_on_horiz_dendrite(self, x_rand, y_rand):
    # x is the "up" direction on the horizontal dendrite.
    self.x = self.x + y_rand

    # (L / (2 pi R)) * (2 pi)
    theta_delta = x_rand / self.HORIZ_DENDRITE_RADIUS
    if np.abs(theta_delta) > np.pi:
      raise ValueError('Extremely large rotation. Consider adjusting k.')
    self.theta += theta_delta
    self.y = self.HORIZ_DENDRITE_RADIUS * math.cos(self.theta)
    self.z = self.HORIZ_DENDRITE_RADIUS * math.sin(self.theta)

  def move_on_horiz_dendrite(self, x_rand=None, y_rand=None):
    x_rand = x_rand or self.k * np.random.randn()
    y_rand = y_rand or self.k * np.random.randn()

    # x is the "up" direction on the horizontal dendrite, so we use y_rand
    # to increment.
    # NOTE: It might be worth changing x_rand, y_rand to right_rand, up_rand.
    new_x = self.x + y_rand
    # If we started behind the dendrite and stayed behind
    # (x < -VERT_DENDRITE_RADIUS) or started past (x > VERT_DENDRITE_RADIUS)
    # and stayed past, we are safe. Also if we are on the underside of the
    # horizontal cylinder (z <= 0), we are safe.
    if (self.z <= 0 or
        (self.x <= -self.VERT_DENDRITE_RADIUS and
        new_x <= -self.VERT_DENDRITE_RADIUS) or
        (self.x >= self.VERT_DENDRITE_RADIUS and
         new_x >= self.VERT_DENDRITE_RADIUS)):
      self._move_on_horiz_dendrite(x_rand, y_rand)
      # Since new_x stayed away from the vertical dendrite, we don't need
      # to change the point_type.
      self.verify_horiz_dendrite()
      return

    # In this case, our x-values either started in the same range as the
    # vertical dendrite or have crossed those x-values.
    self.move_from_horiz_dendrite_to_vert_dendrite(x_rand, y_rand)

  def _get_vert_dendrite_intersection(self, x0, theta0, m,
                                      curved_fn_squared):
    """Finds nearest intersect of the line and bounded curve.

    The line is in the x-theta plane and the curve is the boundary of the
    vertical dendrite.
    """
    def line_fn(theta):
      return x0 + m * (theta - theta0)

    def line_fn_squared(theta):
      return line_fn(theta)**2

    def intersect_fn(theta):
      return curved_fn_squared(theta) - line_fn_squared(theta)

    intersection_theta = None
    try:
      intersection_theta = float(mpmath.findroot(intersect_fn, theta0))
    except ValueError:
      pass

    # Get rid of false intersections due to squaring of arguments, the
    # only valid intersections must come in the domain of the
    # curve that defines the dendrite boundary.
    if (intersection_theta < self.BEGIN_MIXED_HORIZ_THETA or
        intersection_theta > self.END_MIXED_HORIZ_THETA):
      intersection_theta = None

    x_intersect = None
    if intersection_theta is not None:
      x_intersect = line_fn(intersection_theta)

    return intersection_theta, x_intersect

  def move_from_horiz_dendrite_to_vert_dendrite(self, x_rand, y_rand):
    """Attempts to determine move at boundary of dendrites in horizontal.

    May not actually leave the horizontal dendrite, but has been identified
    as a potential to do so.
    """
    x0 = self.x
    theta0 = self.theta
    new_x = x0 + y_rand

    def curved_fn_squared(theta):
      return (self.VERT_DENDRITE_RADIUS_SQUARED -
              self.HORIZ_DENDRITE_RADIUS_SQUARED * mpmath.cos(theta)**2)

    change_regions = False
    intersection_x_rand = x_rand
    intersection_y_rand = y_rand
    if x_rand == 0:
      # We don't move in the theta direction hence the only possible
      # intersection is if theta is already in line with the gap.
      if (self.BEGIN_MIXED_HORIZ_THETA <= theta0 <=
          self.END_MIXED_HORIZ_THETA):
        # In this case, we know either x started behind the vertical dendrite
        # (x < 0, the radius depends on theta) or past it (x > 0).
        if x0 < 0:
          # NOTE: We may need to use np.max(0, curved_fn_squared) here to
          #       avoid accidental square root of negatives.
          x_intersect = - float(np.sqrt(curved_fn_squared(theta0)))
          if new_x > x_intersect:  # On the right.
            change_regions = True
            # Need to reach point of intersection:
            # x0 + intersection_y_rand = x_intersect
            intersection_y_rand = x_intersect - x0
            # intersection_x_rand is already 0
        elif x0 > 0:
          x_intersect = float(np.sqrt(curved_fn_squared(theta0)))
          if new_x < x_intersect:  # On the left.
            change_regions = True
            intersection_y_rand = x_intersect - x0
            # intersection_x_rand is already 0
        else:
          raise ValueError('It is not possible for a point to reach here.')
    else:
      m = y_rand / x_rand
      intersection_theta, x_intersect = self._get_vert_dendrite_intersection(
          x0, theta0, m, curved_fn_squared)

      if x_intersect is not None:
        endpoints = sorted((x0, new_x))
        if endpoints[0] <= x_intersect <= endpoints[1]:
          change_regions = True
          intersection_y_rand = x_intersect - x0
          intersection_x_rand = intersection_theta - theta0
          # These should be the same fraction of the whole. Since
          # (intersection_theta, x_intersect) lies on a line through
          # (theta0, x0) this should be gauranteed.

    # We move up to the point of intersection (if it occurs) first since we
    # are still on the horizontal dendrite.
    self._move_on_horiz_dendrite(intersection_x_rand, intersection_y_rand)

    if not change_regions:
      # NOTE: This may cause problems due to ties.
      self.verify_horiz_dendrite()
    else:
      # Verify before moving within the new region so all the correct
      # values are set.
      try:
        self.verify_vert_dendrite()
      except ValueError:
        if DEBUG:
          print 'Not on vertical dendrite in exact sense.'
        # Failed test to be on vertical dendrite.
        if np.abs(self.y) <= self.VERT_DENDRITE_RADIUS:
          if DEBUG:
            old_x = self.x
            print 'Approximating by updating x by',
          # Force the x-value to make this occur if y is small enough.
          x_sq = self.VERT_DENDRITE_RADIUS_SQUARED - self.y**2
          self.x = np.sign(self.x) * np.sqrt(x_sq)
          if DEBUG:
            print (self.x - old_x)
        else:
          if DEBUG:
            old_x, old_y, old_z = self.x, self.y, self.z
            print 'Full reboot, setting y to the vertical radius.'
          self.x = 0.0
          self.y = np.sign(self.y) * self.VERT_DENDRITE_RADIUS
          z_sq = self.HORIZ_DENDRITE_RADIUS_SQUARED - self.y**2
          self.z = np.sign(self.z) * np.sqrt(z_sq)
          if DEBUG:
            print 'Total changes:', (old_x - self.x, old_y - self.y,
                                     old_z - self.z)
        # After putting the point on the boundary, verify it is
        # on the vertical dendrite.
        self.verify_vert_dendrite()

      remaining_x_rand = x_rand - intersection_x_rand
      remaining_y_rand = y_rand - intersection_y_rand
      remaining_length = np.sqrt(remaining_x_rand**2 + remaining_y_rand**2)

      theta_perp = np.arctan2(self.x, self.y)
      theta_absolute = np.arctan2(y_rand, x_rand)

      theta_new_dir = np.pi / 2 + theta_absolute - theta_perp
      new_x_rand = remaining_length * math.cos(theta_new_dir)
      new_y_rand = remaining_length * math.sin(theta_new_dir)

      self._move_on_vert_dendrite(new_x_rand, new_y_rand)
      self.verify_vert_dendrite()


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


def add_sphere_surface(ax):
  radius = Point.SPHERE_RADIUS
  center_x = 0.0
  center_y = 0.0
  center_z = Point.SPHERE_Z_CENTER
  z_cutoff = Point.SPHERE_CHANGE_POINT

  theta = np.linspace(0, 2 * np.pi, 100)
  # Need Z_CUTOFF = RADIUS cos(phi) + CENTER_Z
  phi_max = np.arccos((z_cutoff - center_z) / radius)
  phi = np.linspace(0, phi_max, 100)

  x = radius * np.outer(np.cos(theta), np.sin(phi)) + center_x
  y = radius * np.outer(np.sin(theta), np.sin(phi)) + center_y
  z = radius * np.outer(np.ones(np.size(theta)), np.cos(phi)) + center_z

  return ax.plot_surface(x, y, z, rstride=4, cstride=4, color='w')


def add_vert_dendrite_surface(ax):
  radius = Point.VERT_DENDRITE_RADIUS
  z_top = Point.SPHERE_CHANGE_POINT
  z_bottom = Point.DENDRITE_CHANGE_BOTTOM

  theta = np.linspace(0, 2 * np.pi, 100)
  z_mesh = np.linspace(z_bottom, z_top, 100)
  x = radius * np.outer(np.cos(theta), np.ones(np.size(z_mesh)))
  y = radius * np.outer(np.sin(theta), np.ones(np.size(z_mesh)))
  z = np.outer(np.ones(np.size(theta)), z_mesh)

  return ax.plot_surface(x, y, z, rstride=4, cstride=4, color='w')


def add_horiz_dendrite_surface(ax):
  radius = Point.HORIZ_DENDRITE_RADIUS
  x_max = Point.VERT_DENDRITE_RADIUS * 1.75
  x_min = - x_max

  theta = np.linspace(np.pi * 0.4, np.pi * 0.6, 100)
  x_mesh = np.linspace(x_min, x_max, 100)
  x = np.outer(np.ones(np.size(theta)), x_mesh)
  y = radius * np.outer(np.cos(theta), np.ones(np.size(x_mesh)))
  z = radius * np.outer(np.sin(theta), np.ones(np.size(x_mesh)))

  return ax.plot_surface(x, y, z, rstride=4, cstride=4, color='w')


def plot_simulation(num_points, num_frames=200, print_frequency=None,
                    interval=30, k=0.01, filename=None):
  points = [Point(0, 0, Point.SPHERE_Z_CENTER + Point.SPHERE_RADIUS, k)
            for i in xrange(num_points)]

  # Attaching 3D axis to the figure
  fig = plt.figure()
  ax = p3.Axes3D(fig)

  # Create lines with a single point as a scatter.
  all_points = [ax.plot([pt.x], [pt.y], [pt.z], c='b', marker='o')[0]
                for pt in points]
  add_sphere_surface(ax)
  add_vert_dendrite_surface(ax)
  add_horiz_dendrite_surface(ax)

  def update_plot(step_num):
    if print_frequency is not None and step_num % print_frequency == 0:
      print 'Step Number:', step_num

    for pt, point_container in zip(points, all_points):
      pt.move()
      color = PointType.COLOR_MAP[pt.point_type]
      point_container.set_color(color)
      point_container.set_data([pt.x], [pt.y])
      point_container.set_3d_properties([pt.z])

    return all_points

  # Setting the axes properties
  ax.set_xlim3d([-Point.SPHERE_RADIUS, Point.SPHERE_RADIUS])
  ax.set_xlabel('X')

  ax.set_ylim3d([-Point.SPHERE_RADIUS, Point.SPHERE_RADIUS])
  ax.set_ylabel('Y')

  ax.set_zlim3d([Point.DENDRITE_CHANGE_BOTTOM,
                 Point.SPHERE_Z_CENTER + Point.SPHERE_RADIUS])
  ax.set_zlabel('Z')

  anim = animation.FuncAnimation(fig, update_plot,
                                 repeat=False, frames=num_frames,
                                 interval=interval, blit=True)
  if filename is not None:
    anim.save(filename, writer='imagemagick_file')
  else:
    plt.show()

  return points


def create_gif():
  plot_simulation(100, num_frames=400, print_frequency=20,
                  filename='100pts_400steps_colored_points.gif')


points = plot_simulation(50, num_frames=200, print_frequency=20)
