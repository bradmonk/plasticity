import numpy as np
import os
import random
import scipy.io


def find_intersection(center0, direction0, center1, direction1):
    """Finds the intersection of two lines in R^3.

    Solves
        c0 + d0 t = c1 + d1 s
    and returns (t, s).

    Args:
        center0: 1D NumPy are with 3 elements.
        direction0: 1D NumPy are with 3 elements.
        center1: 1D NumPy are with 3 elements.
        direction1: 1D NumPy are with 3 elements.

    Raises:
        np.linalg.LinAlgError: If the matrix is not full rank.

    Returns:
        Two real values s and t.
    """
    # c0 + d0 t = c1 + d1 s
    # (-d0) t + (d1) s = c0 - c1
    # [-d0, d1] [t,s]^T = delta
    A = np.array([-direction0, direction1]).T
    delta = center0 - center1
    t, s = np.linalg.solve(A.T.dot(A), A.T.dot(delta))
    return t, s


class Face(object):

    def __init__(self, face_index, mesh_wrapper):
        self.face_index = face_index
        self.mesh_wrapper = mesh_wrapper

        self.points = {}

    @property
    def a(self):
        vertex_index = self.mesh_wrapper.triangles[self.face_index, 0]
        return self.mesh_wrapper.all_vertices[vertex_index, :]

    @property
    def b(self):
        vertex_index = self.mesh_wrapper.triangles[self.face_index, 1]
        return self.mesh_wrapper.all_vertices[vertex_index, :]

    @property
    def c(self):
        vertex_index = self.mesh_wrapper.triangles[self.face_index, 2]
        return self.mesh_wrapper.all_vertices[vertex_index, :]

    @property
    def face_opposite_a(self):
        return self.mesh_wrapper.neighbor_faces[self.face_index, 0]

    @property
    def face_opposite_b(self):
        return self.mesh_wrapper.neighbor_faces[self.face_index, 1]

    @property
    def face_opposite_c(self):
        return self.mesh_wrapper.neighbor_faces[self.face_index, 2]

    @property
    def w1(self):
        return self.mesh_wrapper.face_local_bases[self.face_index, :3]

    @property
    def w2(self):
        return self.mesh_wrapper.face_local_bases[self.face_index, 3:]

    def __str__(self):
        return 'Face(%d)' % self.face_index

    def __repr__(self):
        return str(self)

    def add_point(self, point):
        self.points[point.point_index] = point

    def remove_point(self, point):
        self.points.pop(point.point_index)

    def compute_angle(self, direction):
        """Computes angle of `direction` in current orthogonal coordinates.

        Assumes `direction` lies in plane spanned by w1 and w2. If this
        is not true, we'll still get a number, but it won't make sense.

        NOTE: We could check the validity of this by computing
                  ||direction|| cos(theta) w1 + ||direction|| sin(theta) w2
              and comparing it to `direction`.
        """
        scaled_cosine = self.w1.dot(direction)  # ||direction|| cos(theta)
        scaled_sine = self.w2.dot(direction)  # ||direction|| sin(theta)
        return np.arctan2(scaled_sine, scaled_cosine)

    @staticmethod
    def check_intersection(t, s):
        """Check that intersection is on the right line.

        t and s come from solving
            p + t * d = a + (b - a) * s
        where p is the point we are moving from, d is the normal direction
        we are moving and a and b are vertices of a side of a triangle.

        Check that (t >= 0), i.e. we are actually moving in the direction
        given by `d`. Also checks that (0 <= s <= 1), i.e. that the intersection
        is between `a` and `b`.
        """
        # Allow wiggle room nearby 0.
        if not (np.allclose(t, 0) or 0 <= t):
            return False

        # Allow wiggle room nearby 0.
        if not (np.allclose(s, 0) or 0 <= s):
            return False

        # Allow wiggle room nearby 1.
        if not (np.allclose(s, 1) or s <= 1):
            return False

        return True

    def compute_new_direction_theta(self, fixed_direction,
                                    direction_to_change, new_face_index):
        """Computes the angle of direction of motion in new face.

        Assumes `fixed_direction` is the vector along an edge of the triangle
        and `direction_to_change` should maintain the same angle between that
        vector on the current and next face. Also assumes `new_face_index`
        describes the other face with `fixed_direction` as an edge.
        """
        angle_change = (self.compute_angle(direction_to_change) -
                        self.compute_angle(fixed_direction))

        new_face = self.mesh_wrapper.faces[new_face_index]
        return angle_change + new_face.compute_angle(fixed_direction)

    def move_toward_side(self, list_of_moves,
                         particle_center, particle_direction,
                         side_first_vertex, side_second_vertex,
                         move_length, next_face_index):
        center_side_line = side_first_vertex
        direction_side_line = side_second_vertex - side_first_vertex

        t, s = find_intersection(particle_center, particle_direction,
                                 center_side_line, direction_side_line)

        if not self.check_intersection(t, s):
            # Do nothing; i.e. don't augment `list_of_moves`.
            return

        # If we haven't returned, this side is a valid choice.
        actual_move_length = move_length
        next_face = self.face_index
        # We don't need these attributes if the move stays on the same face.
        remaining_length = theta_new = None
        # If the move takes us past the intersection, we need to change to
        # a different face and can't use the full length of the move.
        if move_length > t:
            actual_move_length = t
            remaining_length = move_length - t
            next_face = next_face_index
            theta_new = self.compute_new_direction_theta(direction_side_line,
                                                         particle_direction,
                                                         next_face)

        next_point = particle_center + actual_move_length * particle_direction
        list_of_moves.append((next_face, next_point,
                              remaining_length, theta_new))

    def choose_move(self, point, move_choices):
        num_choices = len(move_choices)
        if num_choices == 1:
            return move_choices[0]

        # We expect to travel in the direction of either 1 or 2 sides (if
        # we are moving in the direction of a vertex).
        if num_choices != 2:
            raise ValueError('Unexpected number of possible moves on face.')

        first_choice, second_choice = move_choices
        no_move_first = np.allclose(first_choice[1], point)
        no_move_second = np.allclose(second_choice[1], point)

        if no_move_first and not no_move_second:
            return second_choice
        elif not no_move_first and no_move_second:
            return first_choice
        else:
            # NOTE: In the case of multiple moves, we expect that `next_point`
            #       and `L_new` should be identical for all tuples in
            #       `move_choices`, but don't check this here.
            # BEGIN: Temporary statements to show issues with `random.choice`.
            print 'Used random.choice'
            # END: Temporary statements to show issues with `random.choice`.
            # The goal in the case of a tie is to correctly determine which face
            # to go into *through* the vertex. We can do this by taking each
            # triangle and determining the angle contributed at the vertex.
            return random.choice(move_choices)

    def move(self, point, L, theta):
        particle_direction = np.cos(theta) * self.w1 + np.sin(theta) * self.w2
        # In case the point lies on a vertex, we may have more than
        # one possible move.
        move_choices = []

        # Side opposite A is B-->C.
        self.move_toward_side(move_choices, point, particle_direction,
                              self.b, self.c, L, self.face_opposite_a)
        # Side opposite B is C-->A.
        self.move_toward_side(move_choices, point, particle_direction,
                              self.c, self.a, L, self.face_opposite_b)
        # Side opposite C is A-->B.
        self.move_toward_side(move_choices, point, particle_direction,
                              self.a, self.b, L, self.face_opposite_c)

        return self.choose_move(point, move_choices)


class Mesh(object):

    def __init__(self, k, initial_point, initial_face_index,
                 all_vertices, triangles, face_local_bases, neighbor_faces):
        self.k = k
        self.initial_point = initial_point
        self.initial_face_index = initial_face_index

        self.all_vertices = all_vertices
        self.triangles = triangles
        self.face_local_bases = face_local_bases
        self.neighbor_faces = neighbor_faces

        self.faces = {}
        for face_index in xrange(triangles.shape[0]):
            self.faces[face_index] = Face(face_index, self)

        # Set initial_face object based on newly created faces.
        self.initial_face = self.faces[self.initial_face_index]
        # Counter to keep track of points on the mesh.
        self.current_point_index = -1

    def next_index(self):
        self.current_point_index += 1
        return self.current_point_index

    def serialize_mesh(self, filename):
        """Serializes current Mesh object to NPZ file.

        This allows faster loading on subsequent use and
        independence from dolfin. See an example usage in
        sample_code.save_serialized_mesh().
        """
        print 'Saving mesh to', filename
        if self.current_point_index != -1:
            print 'Points on mesh will not be serialized.'

        np.savez(filename, k=self.k, initial_point=self.initial_point,
                 initial_face_index=self.initial_face_index,
                 all_vertices=self.all_vertices, triangles=self.triangles,
                 face_local_bases=self.face_local_bases,
                 neighbor_faces=self.neighbor_faces)

    @classmethod
    def from_mesh(cls, mesh, initial_point, k):
        """Create a Mesh object (this class) from a dolfin mesh."""
        # Import here to prevent cyclic import since `dolfin_mesh_utils`
        # imports this `Mesh` class.
        import dolfin_mesh_utils
        return dolfin_mesh_utils.from_mesh(cls, mesh, initial_point, k)

    @classmethod
    def from_file(cls, filename):
        """Create a mesh from a previously serialized Mesh object.

        Serialization occurs via the serialize_mesh() method above.
        """
        constructor_args = _load_serialized_mesh(filename)
        return cls(*constructor_args)

    def __str__(self):
        return 'Mesh(num_faces=%d)' % len(self.faces)

    def __repr__(self):
        return str(self)


def get_random_components(k):
    """Compute length and angle of random vector.

    The vector is random with N(0, k^2) (normal mean 0, std. dev. k)
    components. Each of the x, y components are indenpently
    distributed N(0, k^2).
    """
    x_rand = k * np.random.randn()
    y_rand = k * np.random.randn()
    L = np.linalg.norm([x_rand, y_rand])
    theta = np.arctan2(y_rand, x_rand)
    return L, theta


class Point(object):

    def __init__(self, mesh_wrapper, point=None, face_index=None):
        self.point_index = mesh_wrapper.next_index()
        # NOTE: We don't need to store this since `self.face.mesh_wrapper`
        #       will also hold this value.
        self.mesh_wrapper = mesh_wrapper

        # Start with null `Face` object. This is required because `change_face`
        # will check the value (here in the constructor and on subsequent
        # calls as the Point moves).
        self.face = None
        if point is None or face_index is None:
            self.point = np.array(mesh_wrapper.initial_point)
            # Set the (previously null) `Face` object.
            self.change_face(mesh_wrapper.initial_face)
        else:
            self.point = point
            # Set the (previously null) `Face` object.
            self.change_face(mesh_wrapper.faces[face_index])

        self.k = mesh_wrapper.k

        self.move_counter = 0
        self.values = [(self.face.face_index, self.point)]

    def __str__(self):
        return 'Point(%d, face=%d)' % (self.point_index,
                                       self.face.face_index)

    def __repr__(self):
        return str(self)

    @property
    def x(self):
        return self.point[0]

    @property
    def y(self):
        return self.point[1]

    @property
    def z(self):
        return self.point[2]

    def change_face(self, face):
        """Updates the face on current object and adds point to face.

        Args:
            face: A Face object.
        """
        if self.face is not None:
            self.face.remove_point(self)

        self.face = face
        self.face.add_point(self)

    def _move(self, L, theta):
        # L_new will be 0 if the point finished moving on the
        # face it was already on. If not, it needs to change faces
        # and keep moving (and L_new > 0).
        next_face, next_point, L_new, theta_new = self.face.move(
            self.point, L, theta)
        self.point = next_point
        if next_face != self.face.face_index:
            self.change_face(self.mesh_wrapper.faces[next_face])
            # Continue to move until we stay on the same face by
            # calling this method recursively.
            self._move(L_new, theta_new)

    def move(self):
        L, theta = get_random_components(self.k)
        self._move(L, theta)

        self.move_counter += 1
        self.values.append((self.face.face_index, self.point))


def run_simulation(num_points, mesh_wrapper, num_steps=200,
                   print_frequency=None):
    points = [Point(mesh_wrapper) for _ in xrange(num_points)]

    for step_num in xrange(1, num_steps + 1):
        if print_frequency is not None and step_num % print_frequency == 0:
            print 'Step Number:', step_num

        for pt in points:
            pt.move()

    return points


def _load_serialized_mesh(filename):
    """Loads a previously serialized Mesh object."""
    print 'Loading mesh data from NPZ file', filename
    npzfile = np.load(filename)

    k = npzfile['k'].item()
    initial_point = npzfile['initial_point']
    initial_face_index = npzfile['initial_face_index'].item()

    all_vertices = npzfile['all_vertices']
    triangles = npzfile['triangles']
    face_local_bases = npzfile['face_local_bases']
    neighbor_faces = npzfile['neighbor_faces']

    return [k, initial_point, initial_face_index,
            all_vertices, triangles, face_local_bases, neighbor_faces]


def convert_mesh_to_matlab(mesh_filename):
    """Loads a previously serialized Mesh object and stores in MATLAB file."""
    (k, initial_point, initial_face_index,
     all_vertices, triangles,
     face_local_bases, neighbor_faces) = _load_serialized_mesh(mesh_filename)
    data = {
        'k': k,
        'initial_point': initial_point,
        'initial_face_index': initial_face_index,
        'all_vertices': all_vertices,
        'triangles': triangles,
        'face_local_bases': face_local_bases,
        'neighbor_faces': neighbor_faces,
    }

    root, ext = os.path.splitext(mesh_filename)
    matlab_filename = root + '.mat'
    scipy.io.savemat(matlab_filename, data)
    print 'Saved', matlab_filename
