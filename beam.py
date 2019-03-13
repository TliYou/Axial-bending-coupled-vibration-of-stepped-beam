"""beam.py

Build up the transverse vibration model for a multi-layered beam.
"""


from math import pi

import numpy as np
from numpy.linalg import det, inv, eig, norm, solve
from scipy.optimize import bisect, minimize, root
from scipy.signal import TransferFunction

from beam_property import ReadOnly
from segment import Segment, Section
from node import BoundaryNode, IntermediateNode
from generalized_solve import generalized_solve

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

TOL_NATURAL_FREQ = 1e-4


def find_root_1d(func, *args, **kwargs):
    return bisect(func, *args, **kwargs)


def find_root(func, *args, **kwargs):
    res = root(func, *args, method='lm', **kwargs)
    return res


def find_min(func, *args, **kwargs):
    res = minimize(func, *args, method='Powell',
                   options={'ftol': 1e-8},  **kwargs)
    return res


class Beam:
    """The transverse vibration model for a beam.

    Parameters
    ----------
    section: Section, optional
        Default section for segment.

    P: float, optional
        Default axial load. Positive means tension.

    name: string, optional
        Name of the Node.
    """

    def __init__(self, section=None, P=0):
        self._profile_default = section
        self._P_default = P
        self._coordinates = []
        self._segments = []
        self._nodes_intermediate = []
        self._node_left_boundary = None
        self._node_right_boundary = None
        self._omega = 0.0
        self._initialized = False

    def is_initialized(self):
        return self._initialized

    # called by segments or nodes
    def uninitialize(self):
        self._initialized = False

    def initialize(self):
        if not self.is_initialized():
            self._initialized = True
            try:
                # 必须先初始化segment，因为node的初始化依赖于segment的值。
                for s in self.segments:
                    s.initialize()
                for n in self.nodes:
                    n.initialize()

                self._S_nm1_j = None                     # S_{n-1,j}
                self._S_i_nm1 = None            # S_{i,n-1}

                self._H = None         # 频响矩阵
            except:
                self._initialized = False
                raise

    @property
    def segments(self):
        return self._segments

    @property
    def nodes(self):
        return ([self._node_left_boundary] +
                self._nodes_intermediate +
                [self._node_right_boundary])

    @property
    def left_boundary(self):
        return self._node_left_boundary

    @property
    def right_boundary(self):
        return self._node_right_boundary

    @property
    def H(self):
        return self._get_H()

    def append_segment(self, pos, section=None, P=None):
        self.uninitialize()
        section = self._section_default if section is None else section
        P = self._P_default if P is None else P
        self._coordinates.append(pos)
        self._segments.append(Segment(self, section, P,
                                      str(len(self._segments))))
        if len(self._segments) == 1:
            self._node_left_boundary = BoundaryNode(
                self, 0, right=self._segments[0], name='0')
        else:
            pos_left = self._coordinates[-2]
            node = IntermediateNode(self, pos_left,
                                    left=self._segments[-2],
                                    right=self._segments[-1],
                                    name=str(len(self._segments) - 1))
            self._nodes_intermediate.append(node)
        self._node_right_boundary = BoundaryNode(self, pos,
                                                 left=self._segments[-1],
                                                 name=str(len(self._segments)))
        self._displacement = np.zeros(
            (len(self._segments) + 1) * 3, dtype='complex')
        self._force = np.zeros((len(self._segments) + 1) * 3, dtype='complex')
        self._displacement[:] = np.nan
        self._force[:] = 0
        return self

    def reset_P(self, Ps):
        """Set axial load on each segment."""
        # 不需要调用uninitialize，因为它会在set_P中被调用
        for i,s in enumerate(self._segments):
            s.set_P(Ps[i])

    @property
    def omega(self):
        return self._omega

    @property
    def frequency(self):
        return self._omega / 2 / pi

    def set_omega(self, omega):
        """Set angular frequency of the beam. (rad/s)"""
        self._omega = omega
        self.uninitialize()

    def set_frequency(self, freq):
        """Set angular frequency of the beam. (Hz)"""
        # wrapper for set_omega
        omega = freq * 2 * pi
        self.set_omega(omega)

    def _matrix_boundary_condition(self):
        B = np.zeros([6, 6], dtype='complex')
        B[:3, :] = self._node_left_boundary.constraint_matrix
        B[3:, :] = self._node_right_boundary.constraint_matrix
        for i, n in enumerate(self._nodes_intermediate[::-1]):
            T = n.transfer_matrix
            B[3:, :] = B[3:, :].dot(T)
        return B

    def _det_imag(self, omega):
        self.set_omega(omega)
        B = self._matrix_boundary_condition()
        y = det(B)
        return y.imag

    def _det_real(self, omega):
        self.set_omega(omega)
        B = self._matrix_boundary_condition()
        y = det(B)
        return y.real

    def _det_real_imag(self, omega):
        self.set_omega(omega)
        B = self._matrix_boundary_condition()
        y = det(B)
        return y.real + y.imag

    def _det_abs(self, omegari):
        r, i = omegari
        self.set_omega(r + 1j * i)
        B = self._matrix_boundary_condition()
        y = det(B)
        return (abs(y))

    def _det(self, omegari):
        r, i = omegari
        self.set_omega(r + 1j * i)
        B = self._matrix_boundary_condition()
        y = det(B)
        return [y.real, y.imag]

    def natural_frequency_real(self, start=0, end=500, nInterval=500):
        """计算系统固有频率（仅考虑实部）。若系统有阻尼，则此函数仅可得到近似解。"""
        start *= 2 * pi
        end *= 2 * pi
        dw = (end - start) / nInterval
        if start < dw / 1e3:
            start = dw / 1e3
        points = np.linspace(start, end, nInterval + 1)
        y = np.zeros(nInterval + 1)
        y[0] = self._det_real(points[0])
        naturalFreq = []
        for i in range(nInterval):
            y[i + 1] = self._det_real(points[i + 1])
            if (y[i] <= 0 and y[i + 1] > 0) or (y[i] > 0 and y[i + 1] <= 0):
                #                print(points[i],points[i+1],y[i],y[i+1])
                f = find_root_1d(self._det_real, points[i], points[i + 1])
                naturalFreq.append(f)
#        plot(points/2/pi, y)
        return np.array(naturalFreq) / 2 / pi

    def natural_frequency(self, start=0, end=500, nInterval=500,
                          method='minimize', tol=TOL_NATURAL_FREQ):
        """计算系统的固有频率（考虑复模态）。
        此方法先通过natural_frequency_real求取估计固有频率，再以此为初始值，求取复模态频率。
        method为'minimize'或'root'，分别对应通过求解优化问题求根和直接求根。"""
        # root方法算虚部算不准
        # 阻尼很大时，minimize求不出根。
        if method not in ['minimize', 'root']:
            raise ValueError('method should be "minimize" or "root"')
        freqs_estimate = self.natural_frequency_real(start, end, nInterval)
        natural_frequencies = []
        for fe in freqs_estimate:
            x0 = [fe * 2 * pi, 0]
            if method == 'minimize':
                res = find_min(self._det_abs, x0)
            else:
                res = find_root(self._det, x0)
            if res.success:
                x = res.x[0] + res.x[1] * 1j
                if (len(natural_frequencies) != 0 and
                        abs(x - natural_frequencies[-1]) < tol):
                    pass
                else:
                    natural_frequencies.append(x)
            else:
                print('miss a root', res)
        return np.array(natural_frequencies) / 2 / pi

    def _get_all_C_from_C0(self, C0, nodeID=None, fx=0, fy=0, m=0):
        C = C0
        C_table = [C0]
        for i in range(len(self._nodes_intermediate)):
            C = (self._nodes_intermediate[i].
                 transfer_matrix.dot(C))
            if i + 1 == nodeID:
                C -= (inv(self._nodes_intermediate[i].T_right).
                      dot(np.array([0, 0, 0, fx, fy, m])))
            C_table.append(C)
        return C_table

    def _response(self, C, npoints=50):
        """通过给定各个segment的C来求取梁的响应。"""
        points = []
        response = []
        x_start = 0
        for i, seg in enumerate(self._segments):
            x_end = self._coordinates[i]
            xs = np.linspace(x_start, x_end, npoints + 1)
            x_start = x_end
            for x in xs:
                r = seg.get_response_matrix(x)
                points.append(x)
                response.append(r.dot(C[i]))
        x = np.array(points)
        y = np.array(response).T
        freq = self.frequency
        return Result(freq, x, y)

    def _node_response(self, C):
        """通过给定各个segment的C来求取梁在节点上的响应"""
        y = np.zeros([6, len(self._segments) + 1], dtype='complex')
        for i in range(len(self._segments)):
            pos = ([0] + self._coordinates)[i]
            y[:, i] = self._segments[i].get_response_matrix(pos).dot(C[i])
        pos = self._coordinates[-1]
        y[:, -1] = self._segments[-1].get_response_matrix(pos).dot(C[-1])
        x = [0] + self._coordinates
        return Result(self.frequency, x, y)

    def mode_shape(self, freq=None, npoints=50, tol=1e-6):
        """计算给定频率下的振型。"""
        if freq is not None:
            self.set_frequency(freq)
        B = self._matrix_boundary_condition()
        eigval, eigvect = eig(B)
        idx = np.argmin(abs(eigval))
        if tol < abs(eigval[idx]) / abs(eigval).max():
            raise ValueError("freq = %s(Hz) is not the "
                             "natural frequency of the beam" % (freq / 2 / pi))
        C0 = eigvect[:, idx]
        C = self._get_all_C_from_C0(C0)
        res = self._response(C, npoints)
        k = norm((res.u**2 + res.v**2)**.5, np.inf)
        res = Result(
            freq, res.x, [res.u / k, res.v / k, res.theta / k, res.N / k, res.Q / k, res.M / k])
        return res

    def modal(self, start=0, end=500, nInterval=500, npoints=50):
        """计算梁的模态。"""
        natural_frequency = self.natural_frequency(start, end, nInterval)
        reses = []
        for f in natural_frequency:
            res = self.mode_shape(f, npoints)
            reses.append(res)
        return reses

    def harmonic_response(self, nodeID, fx=0,fy=0,m=0, npoints=50):
        C0 = self._get_C0(nodeID, fx,fy,m)
        C_table = self._get_all_C_from_C0(C0, nodeID, fx,fy,m)
#        print(self._node_response(C_table))
        res = self._response(C_table, npoints)
        return res

    def transfer_function(self, node_input, type_input, node_output,
                          type_output, freqs):
        """Calculate the transfer function of given input and output node.
        Parameters
        ----------
        node_input: int
            Select which node to apply force.

        type_input: str
            Force or moment apply to the node, where 'Fx' stands for axial
            force, 'Fy' stands for transverse force and 'M' stands for moment.

        node_output: int
            Select which node to get displacement.

        type_displacement: str
            Displacement or rotational angle to calculate, where 'U' means
            axial displacement, 'V' stands for transverse displacement and
            'Theta' means rotational angle.

        freqs: array-like object
            The frequencies needed for transfer function.

        Returns
        -------
        tf: np.array
            Transfer function.
        """
        fx, fy, m = 0, 0,0
        if type_input.upper() == 'FX':
            fx = 1
        elif type_input.upper() == 'FY':
            fy = 1
        elif type_input.upper() == 'M':
            m = 1
        else:
            raise ValueError("type_force should be 'Fx', 'Fy', or 'M'.")

        h = []
        for freq in freqs:
            self.set_frequency(freq)
            C0 = self._get_C0(node_input,fx,fy,m)
            C_table = self._get_all_C_from_C0(C0, node_input, fx, fy, m)
            res = self._node_response(C_table)
            hi = getattr(res, type_output)[node_output]
            h.append(hi)
        return np.array(h)

    def _form_S_nm1_j(self):
        """Called by _S"""
        S = [None] * len(self._segments)
        S[-1] = np.eye(6, dtype='complex')
        for i in range(len(self._segments) - 1):
            j = len(self._segments) - 2 - i
            T = self.nodes[j + 1].transfer_matrix
            S[j] = S[j + 1].dot(T)
        invS = []
        for i in S:
            invS.append(inv(i))
        self._S_nm1_j = S                     # S_{n-1,j}
        self._S_i_nm1 = invS                  # S_{i,n-1}
        return S, invS

    def _S(self, i, j):
        if not self.is_initialized():
            self.initialize()
        if self._S_nm1_j is None:
            self._form_S_nm1_j()
        return self._S_i_nm1[i].dot(self._S_nm1_j[j])

    def _get_C0(self, nodeID, fx=0, fy=0, m=0):
        if 0 < nodeID < len(self._segments):
            node = self.nodes[nodeID]
            A1 = self._node_left_boundary.constraint_matrix
            A2 = self._node_right_boundary.constraint_matrix.dot(
                self._S(len(self._segments) - 1, 0))
            A = np.concatenate([A1, A2], axis=0)
            b2 = self._node_right_boundary.constraint_matrix.dot(
                    self._S(len(self._segments) - 1, nodeID))
            b2 = b2.dot(inv(node.T_right))
            b2 = b2.dot(np.array([0, 0,0, fx,fy,m]))
            b = np.array([0, 0,0, b2[0], b2[1],b2[2]], dtype='complex')
        if nodeID == 0:
            A1 = self._node_left_boundary.constraint_matrix
            A2 = self._node_right_boundary.constraint_matrix.dot(
                self._S(len(self._segments) - 1, 0))
            A = np.concatenate([A1, A2], axis=0)
            b = -np.array([fx,fy,m,0, 0, 0], dtype='complex')
        if nodeID == len(self._segments):
            A1 = self._node_left_boundary.constraint_matrix
            A2 = self._node_right_boundary.constraint_matrix.dot(
                self._S(len(self._segments) - 1, 0))
            A = np.concatenate([A1, A2], axis=0)
            b = np.array([0, 0,0,fx,fy,m], dtype='complex')
        return inv(A).dot(b)

    def _get_H(self):
        if not self.is_initialized():
            self.initialize()
        if self._H is not None:
            return self._H
        n = len(self._segments) + 1
        H = np.zeros([3 * n, 3 * n], dtype='complex')
        for i in range(n):
            # 单位力激励
            C = self._get_all_C_from_C0(self._get_C0(i, 1, 0,0), i, 1, 0,0)
            y = self._node_response(C).displacement
            H[:, i] = y.reshape(-1)
            C = self._get_all_C_from_C0(self._get_C0(i, 0, 1,0), i, 0, 1,0)
            y = self._node_response(C).displacement
            H[:, n+i] = y.reshape(-1)
            C = self._get_all_C_from_C0(self._get_C0(i, 0,0,1), i, 0,0,1)
            y = self._node_response(C).displacement
            H[:, 2*n+i] = y.reshape(-1)
        self._H = H
        return H

    def set_force(self, nodeID, fx=0, fy=0, m=0):
        self._force[nodeID] = fx
        self._force[nodeID + len(self._segments) + 1] = fy
        self._force[nodeID + 2 * len(self._segments) + 2] = m
        return self

    def set_displacement(self, nodeID, u=np.nan, v=np.nan, theta=np.nan):
        self._displacement[nodeID] = u
        self._displacement[nodeID + len(self._segments) + 1] = v
        self._displacement[nodeID + 2 * len(self._segments) + 2] = theta
        return self

    def _solve(self):
        n = len(self._segments) + 1
        disp, force = generalized_solve(self._displacement,
                                        self.H, self._force,
                                        leastsq=True)
        u = disp[:n]
        v = disp[n:n*2]
        theta = disp[2*n:]
        fx = force[:n]
        fy = force[n:2*n]
        m = force[2*n:]
        return u, v, theta, fx, fy, m

    def solve(self):
        x = [0] + self._coordinates
        return Result(self.frequency, x, self._solve())

    def response(self, npoints=50):
        u, v, theta, fx, fy, m = self._solve()
        C = []
        pos = [0] + self._coordinates
        for i in range(len(self._segments)):
            BL = self.segments[i].get_response_matrix(pos[i])[:3, :]
            BR = self.segments[i].get_response_matrix(pos[i + 1])[:3, :]
            B = np.concatenate([BL, BR], axis=0)
            C.append(
                solve(B, np.array([u[i], v[i],theta[i], u[i + 1], v[i + 1], theta[i+1]])))
        return self._response(C, npoints)


class Result(ReadOnly):
    def __init__(self, freq, x, y):
        super().__init__()
        self.freq = freq
        self.x = x
        self.u = y[0]
        self.v = y[1]
        self.theta = y[2]
        self.n = y[3]
        self.q = y[4]
        self.m = y[5]
        self._initialize()

    def __getattr__(self, name):
        return super().__getattribute__(name.lower())

    @property
    def displacement(self):
        return np.array([self.u, self.v, self.theta])

    @property
    def force(self):
        return np.array([self.n, self.q, self.m])

    def __repr__(self):
        return 'Result object @{freq:.2f}Hz'.format(freq=self.freq.real)


def plotbeam(res, target='u', scalar=None, type=None, title='', figsize=None):
    """绘制梁的振动计算结果。
    type可以为real,imag,abs或者angle。target可以是U,U1,U2,UR或F,F1,F2,FR。（参考
    ABAQUS viewer）"""
    if type is None:
        type = 'real'
        for i in res.u, res.v:
            if np.imag(i).any():
                import sys
                sys.stderr.write('plot real part discarding imaginary part.\n')
                break
    if type not in ['abs', 'real', 'imag', 'angle']:
        raise ValueError('"type" should be "real", "imag", "abs" or "angle"')
    u = eval('np.{func}(res.u)'.format(func=type))
    v = eval('np.{func}(res.v)'.format(func=type))
    theta = eval('np.{func}(res.theta)'.format(func=type))
    n = eval('np.{func}(res.n)'.format(func=type))
    q = eval('np.{func}(res.q)'.format(func=type))
    m = eval('np.{func}(res.m)'.format(func=type))

    x = res.x
    target = target.upper()
    if target == 'U':
        z = np.sqrt(u**2 + v**2)
    elif target == 'U1':
        z = u
    elif target == 'U2':
        z = v
    elif target == 'UR':
        z = theta
    elif target == 'F':
        z = sqrt(res.N**2 + res.Q**2)
    elif target == 'F1':
        z = N
    elif target == 'F2':
        z = Q
    elif target == 'FR':
        z = M
    else:
        raise ValueError('"target" should be "U", "U1", "U2", "UR" or "F", '
                         '"F1", "F2", "FR"')
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, aspect='equal')
    if scalar is None:
        scalar = 1 / np.max(abs(u)**2+abs(v)**2)**0.5 * 0.3
#        scalar = 1
    u = u * scalar
    v = v * scalar
    points = np.array([x + u, v]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap=plt.get_cmap('jet'),
                        norm=plt.Normalize(min(z), max(z)))
    lc.set_array(z)
    lc.set_linewidth(3)
    ax.add_collection(lc)
    ax.set_xlim(min(x+u), max(x+u))
    ax.set_ylim(min(v), max(v))
    ax.set_title(title)
    return fig


def plotres(res, target='v', type=None, title=''):
    """绘制梁的振动计算结果。
    type可以为real,imag,abs或者angle。target可以是u,v,theta,N,Q或者M。"""
    x = res.x
    y = getattr(res, target)
    if type is None:
        type = 'real'
        if np.imag(y).any():
            import sys
            sys.stderr.write('plot real part discarding imaginary part.\n')
    if type not in ['abs', 'real', 'imag', 'angle']:
        raise ValueError('"type" should be "real", "imag", "abs" or "angle"')
    y = eval('np.{func}(y)'.format(func=type))
    if type == 'angle':
        y *= 180 / np.pi
    fig = plt.figure()
    ax = fig.add_subplot(111)
#    points = np.array([x, y]).T.reshape(-1, 1, 2)
#    segments = np.concatenate([points[:-1], points[1:]], axis=1)
#    lc = LineCollection(segments, cmap=plt.get_cmap('jet'),
#                        norm=plt.Normalize(0, max(abs(y))))
#    lc.set_array(abs(y))
#    lc.set_linewidth(3)
#    ax.add_collection(lc)
    ax.plot(x,y)
    ax.set_xlim(0, max(x))
    ax.set_ylim(min(y), max(y))
    ax.set_title(title)
    return fig


if __name__ == '__main__':
    from pylab import *
    from cmath import sqrt
    from beam_property import *
    from segment import Segment
    from matplotlib.collections import LineCollection
    from matplotlib.colors import ListedColormap, BoundaryNorm
    E = 210e9 #*(1+0.02j)
    nu = 0.3
    b = 0.1
    sect1 = Section([-0.05,0.05], 210e9, 0.3, 0.1, 7800, 5/6, 'section-1')
    sect2 = Section([-0.05,0,0.05], [210e9,69e9], [0.3,0.33], [0.1,0.1], [7800,2700], [0.85,0.85], 'section-2')
    beam1 = Beam(section=sect2)
    beam1.append_segment(1.0, sect2, P=0)
    beam1.nodes[0].free()
    beam1.nodes[1].free()
    beam1.set_frequency(100)
    beam1.initialize()
    beam2 = Beam(section=sect2)
    beam2.append_segment(1.0, sect2, P=0)
    beam2.append_segment(2.0, sect2, P=0)
    beam2.append_segment(3.0, sect2, P=0)
    beam2.nodes[0].free()
    beam2.nodes[1].free().add_spring(1e6, 1e6, 1e6)
    beam2.nodes[2].free()
    beam2.nodes[3].free()
    beam2.set_frequency(100)
    beam2.set_force(0,1,0,0)
    beam2.set_force(1,0,1,0)
    beam2.set_force(2,0,0,1)
    beam2.set_force(3,1,1,1)
    beam2.initialize()
