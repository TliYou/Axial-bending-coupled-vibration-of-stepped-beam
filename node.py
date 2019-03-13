"""Node.py

The node of a beam.
Designed for Axial-bending coupled beam.
"""


import numpy as np
from numpy.linalg import inv, det, norm
from numpy.linalg import solve

INF = 1e12


def inv2(A):
    return inv(A)


class Node:
    """Define a Node of a beam.

    A beam can be divided into several segments, each segment has the same
    material and profile. Each node is the conjection of the segments.

    Parameters
    ----------
    beam: Beam
        To which beam the segment belongs.

    pos: float
        coordinate of the node.

    left: Segment, optional
        Segment to the left of the node.

    right: Segment, optional
        Segment to the right of the node.

    name: string, optional
        Name of the Node.
    """

    def __init__(self, beam, pos, left=None, right=None, name=''):
        self._beam = beam
        self._pos = pos
        self._left = left
        self._right = right
        self.name = name
        self._T_left = None
        self._T_right = None
        self.set_K()

    def __repr__(self):
        return 'Node %s' % self.name

    @property
    def beam(self):
        return self._beam

    @property
    def left(self):
        return self._left

    @property
    def right(self):
        return self._right

    @property
    def omega(self):
        return self._beam.omega

    @property
    def frequency(self):
        return self._beam.frequency

    def is_initialized(self):
        return self.beam.is_initialized()

    def initialize(self):
        self.calculate_T_left()
        self.calculate_T_right()

    def fixed(self, k=INF):
        """Fix the node.

        Apply the boundary condition 'fixed' to this node. It means u, v and
        theta are zero.
        """
        def K_fixed(w):
            K = np.diag([k, k, k]).astype('complex')
            return K
        self.set_K(K_fixed)
        return self

    def free(self):
        """Free the node.

        Apply the boundary condition 'free' to this node. It means the
        fx, fy and M on this node are zero.
        """
        def K_free(w):
            K = np.zeros([3, 3], dtype='complex')
            return K
        self.set_K(K_free)
        return self

    def pinned(self, k=INF):
        """Pin the node.

        Apply the boundary condition 'pinned' to this node. It means the
        u, v and M on this node are zero.
        """
        def K_pinned(w):
            K = np.diag([k, k, 0]).astype('complex')
            return K
        self.set_K(K_pinned)
        return self

    def set_H(self, func_H=None):
        """Set the flexibility function matrix of the boundary condition of the
        node.

        .. math:: (U V \Theta)^T = H \cdot (N_B Q_B M_B)^T

        Paremeters
        ----------
        func_H: function or None
            If func_H is a function, it should has one input w (rad/s) and
            output a 3x3 array which is the flexibility matrix of given
            frequency w. If func_H is None, it equals set the node free.
            (default = None)

        Returns
        -------
        self: Node
            return self.
        """
        if func_H is None:
            self.free()
        else:
            def func_K(w):
                return inv2(func_H(w))
            self.set_K(func_K)
        return self

    def set_K(self, func_K=None):
        """Set the stiffness function matrix of the boundary condition of the
        node.

        .. math:: (N_B Q_B M_B)^T = K \cdot (U V \Theta)^T

        Paremeters
        ----------
        func_K: function or None
            If func_K is a function, it should has one input w (rad/s) and
            output a 3x3 array which is the stiffness matrix of given frequency
            w. If func_K is None, it equals set the node free. (default = None)

        Returns
        -------
        self: Node
            return self.
        """
        # only necessary here because all other methods used for setting
        # boundary conditions call set_K
        self.beam.uninitialize()

        if func_K is None:
            self.free()
        else:
            self._func_K = func_K
        return self

    def parallel(self, func_K):
        """Add a sub-system as boundary condition of the node.

        The total stifness is the parallel of the original stiffness and the
        new stiffness.

        .. math:: K = \frac{K_1K_2}{K_1+K_2}

        Paremeters
        ----------
        func_K: function
            It should be a function with one input w (rad/s) and outputs a 2x2
            array which is the stiffness matrix of given frequency w.

        Returns
        -------
        self: Node
            return self.
        """
        func_K_1 = self._func_K
        func_K_2 = func_K

        def new_func(w):
            return func_K_1(w) + func_K_2(w)
        self.set_K(new_func)
        return self

    def add_spring(self, kx=0, ky=0, kr=0):
        """Add a spring to the node.

        Add a spring paralled to the original boundary condition.

        Paremeters
        ----------
        kx: float
            Axial stiffness. (default=0)

        ky: float
            Translational stiffness. (default=0)

        kr: float
            Rotational stiffness. (default=0)

        Returns
        -------
        self: Node
            return self.
        """
        def func_K(w):
            K = -np.diag([kx, ky, kr]).astype('complex')
            return K
        self.parallel(func_K)
        return self

    def add_damping(self, cx=0, cy=0, cr=0):
        """Add a damper to the node.

        Add a dashpot paralled to the original boundary condition.

        Paremeters
        ----------
        cx: float
            Axial damping. (default=0)

        cy: float
            Translational damping. (default=0)

        cr: float
            Rotational damping. (default=0)

        Returns
        -------
        self: Node
            return self.
        """
        def func_K(w):
            coeff = 1j * w
            K = -np.diag([cx, cy, cr]) * coeff
            return K.astype('complex')
        self.parallel(func_K)
        return self

    def add_mass(self, mt=0, mr=0):
        """Add a mass to the node.

        Paremeters
        ----------
        mt: float
            Mass. (default=0)

        mr: float
            Moment of inertia. (default=0)

        Returns
        -------
        self: Node
            return self.
        """
        def func_K(w):
            coeff = - w ** 2
            K = -np.diag([mt, mt, mr]) * coeff
            return K.astype('complex')
        self.parallel(func_K)
        return self

    def get_K(self):
        """Get the stiffness function matrix of the boundary condition of the
        node."""
        return self._func_K

    @property
    def T_left(self):
        """Get the matrix T on the left side of the node after
        initialization."""
        if not self.is_initialized():
            self.beam.initialize()
        return self._T_left

    @property
    def T_right(self):
        """Get the matrix T on the right side of the node after
        initialization."""
        if not self.is_initialized():
            self.beam.initialize()
        return self._T_right

    def calculate_T_left(self):
        """Get the matrix T on the left side of the node.

        Returns
        -------
        T: np.array
            Result. Note that the omega of self.left should be set correctly
            manually.
        """
        if self.left is None:
            self._T_left = None
            return None
        pos = self._pos
        R = self.left.get_response_matrix(pos)
        Ru = R[:3]
        Rf = R[3:]
#        Ru = self.left.get_displacement_response_matrix(pos)
#        Rf = self.left.get_force_response_matrix(pos)
        T = np.zeros([6, 6], dtype='complex')
        T[:3] = Ru
        T[3:] = Rf
        self._T_left = T
        return T

    def calculate_T_right(self):
        """Get the matrix T on the right side of the node.

        Returns
        -------
        T: np.array
            Result. Note that the omega of self.right should be set correctly
            manually.
        """
        if self.right is None:
            self._T_right = None
            return None
        pos = self._pos
        R = self.right.get_response_matrix(pos)
        Ru = R[:3]
        Rf = R[3:]
#        Ru = self.right.get_displacement_response_matrix(pos)
#        Rf = self.right.get_force_response_matrix(pos)
        K = self.get_K()(self.omega)
        T = np.zeros([6, 6], dtype='complex')
        T[:3] = Ru
        T[3:] = Rf + K.dot(Ru)
        self._T_right = T
        return T


class IntermediateNode(Node):
    """Define an intermediate Node of a beam.

    A beam can be divided into several segments, each segment has the same
    material and profile. Each node is the conjection of the segments.

    Parameters
    ----------
    pos: float
        coordinate of the node.

    left: Segment, optional
        Segment to the left of the node.

    right: Segment, optional
        Segment to the right of the node.

    name: string, optional
        Name of the Node.
    """

    def __init__(self, beam, pos, left, right, name=''):
        super().__init__(beam, pos, left, right, name)
        self._transfer_matrix = None

    def __repr__(self):
        return "IntermediateNode %s" % self.name

    def initialize(self):
        """Called by self.beam."""
        super().initialize()
        self.calculate_transfer_matrix()

    @property
    def transfer_matrix(self):
        """Get transfer matrix after initialization."""
        if not self.is_initialized():
            self.beam.initialize()
        return self._transfer_matrix

    def calculate_transfer_matrix(self):
        """Calculate transfer matrix.

        Calculate the transfer matrix from segment i-1 to segment i.

        Returns
        -------
        T: np.array
            Transfer matrix.
        """
        TL = self.T_left
        TR = self.T_right
        ####################
        T = self._solve(TR, TL)   # 此处存在计算精度问题！！！
        ####################
        self._transfer_matrix = T
        return T

    def _solve(sefl,A,B):
        D = []
        for row in A:
            D.append(np.sum(abs(row)))
        D = np.diag(D)
        An = solve(D,A)
        X = solve(An,solve(D,B))
        return X


class BoundaryNode(Node):
    """Define a boundary Node of a beam.

    Define left or right boundary node of a beam.

    Parameters
    ----------
    pos: float
        coordinate of the node.

    left: Segment, optional
        Segment to the left of the node.

    right: Segment, optional
        Segment to the right of the node.

    name: string, optional
        Name of the Node.

    Notes
    -----
    Either left or right should be defined, or ValueError will be raised.
    """

    def __init__(self, beam, pos, left=None, right=None, name=''):
        if left is None and right is None:
            raise ValueError('except left or right')
        if not (left is None or right is None):
            raise ValueError('except left or right')
        super().__init__(beam, pos, left, right, name)
        self._constraint_matrix = None

    def __repr__(self):
        return "BoundaryNode %s" % self.name

    def initialize(self):
        """Called by self.beam."""
        super().initialize()
        self.calculate_constraint_matrix()

    @property
    def constraint_matrix(self):
        """Get constraint matrix after initialization."""
        if not self.is_initialized():
            self.beam.initialize()
        return self._constraint_matrix

    # needs another version different from intermediate node because all
    # external force are assumed to applied on the right side of a node. But
    # when it comes to the boundary node on the right side, the force must
    # be applied on the left side of a node.
    def calculate_T_left(self):
        """Get the matrix T on the left side of the node. (for boundary node
        on the right only)

        Returns
        -------
        T: np.array
            Result. Note that the omega of self.right should be set correctly
            manually.
        """
        if self.left is None:
            self._T_left = None
            return None
        pos = self._pos
        R = self.left.get_response_matrix(pos)
        Ru = R[:3]
        Rf = R[3:]
#        Ru = self.left.get_displacement_response_matrix(pos)
#        Rf = self.left.get_force_response_matrix(pos)
        K = self.get_K()(self.omega)
        T = np.zeros([6, 6], dtype='complex')
        T[:3] = Ru
        T[3:] = Rf - K.dot(Ru)
        self._T_left = T
        return T

    def calculate_constraint_matrix(self):
        """Calculate constraint matrix.

        Calculate the constraint matrix on the left or right boundary.

        Returns
        -------
        T: np.array
            Constraint matrix for boundary conditions.
        """
        # 左侧边界条件
        if self.left is None:
            T = self.T_right
        # 右侧边界条件
        elif self.right is None:
            T = self.T_left
        else:
            raise Exception
        CM = T[3:, :]
        self._constraint_matrix = CM
        return CM
