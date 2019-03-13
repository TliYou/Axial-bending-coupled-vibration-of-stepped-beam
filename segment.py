"""segment.py

The Segment of a multi-layered beam.
"""

from math import pi
from cmath import sqrt
import numpy as np

from beam_property import Section

EPS = 1e-6


class Segment:
    """Define a segment of a beam.

    A beam can be divided into several segments, each segment has the same
    material and profile.

    Parameters
    ----------
    beam: Beam
        To which beam the segment belongs.

    section: Secion
        section of the segment

    P: float, optional
        Axial load on the segment. if P > 0, it means tension. (default=0)

    name: string, optional
        Name of the Segment.
    """

    def __init__(self, beam, section, P=0, name=''):
        self._beam = beam
        self.section = section
        self._P = P
        self.name = name

    @property
    def P(self):
        return self._P

    @property
    def beam(self):
        return self._beam

    @property
    def omega(self):
        return self._beam.omega

    @property
    def frequency(self):
        return self._beam.frequency

    def is_initialized(self):
        return self.beam.is_initialized()

    def set_P(self, P):
        """Set axial load on the segment"""
        self.beam.uninitialize()
        self._P = P

    def __repr__(self):
        return 'Segment %s' % self.name

    def __getattr__(self, name):
        if hasattr(self.section, name):
            return getattr(self.section, name)
        else:
            raise AttributeError(
                "'Segment' object has no attribute '%s'" % name)

    def _get_polycoeffs(self):
        """计算特征多项式。（计算的中间变量）"""
        s = self
        P = s._P
        w = s.omega
        m = s.m
        R = s.R
        J = s.J
        kuu = s.kuu
        kur = s.kur
        krr = s.krr
        kG = s.kG
        a11 = np.poly1d([-kuu, w**2 * m])
        a13 = np.poly1d([kur, -w**2 * R])
        a22 = np.poly1d([-kG - P, w**2 * m])
        a23a32 = np.poly1d([kG**2, 0])     # a23 * a32
        a31 = a13
        a33 = np.poly1d([-krr, w**2 * J - kG])

        p = a11 * a22 * a33 - a13 * a22 * a31 - a23a32 * a11
        return p

    def _get_lambda(self):
        """Calculate eigenvalues."""
        p = self._get_polycoeffs()
        roots = np.roots(p).astype('complex')
        self.k = np.zeros([6], dtype='complex')
        self.k[:3] = np.sqrt(roots)
        self.k[3:] = - self.k[:3]
        return self.k

    def _get_K(self):
        """计算位移间的比例系数（即非0解的系数比）。"""
        s = self
        P = s._P
        w = s.omega
        m = s.m
        R = s.R
        J = s.J
        kuu = s.kuu
        kur = s.kur
        krr = s.krr
        kG = s.kG
        k = self.k
        a11 = -k**2*kuu + w**2 * m
        a13 = k**2*kur-w**2 * R
        a22 = -k**2 * kG + w**2 * m - P * k**2
        a23 = 1j * k * kG

        # define the non-trial solution
        Cu, Cv, Cr = -a13 * a22, -a23 * a11, a11 * a22   # Kr stands for K_\theta
        # needs special treatment if a11 = a31 = 0
        for i in range(len(a11)):
            if abs(a13[i]) < EPS and abs(a11[i]) == np.min(abs(a11)):
                Cu[i] = 1
                Cv[i] = 0
                Cr[i] = 0

        # 归一化
        for i in range(6):
            k = max([abs(Cu[i]),abs(Cv[i]),abs(Cr[i])])
            Cu[i] /= k
            Cv[i] /= k
            Cr[i] /= k

        self.Cu = Cu
        self.Cv = Cv
        self.Cr = Cr
        return Cu, Cv, Cr

    def initialize(self):
        """Update all intermediate variables used in calculation.
        Called by self.beam.
        """
        self._get_lambda()
        self._get_K()

    def get_response_matrix(self,x):
        if not self.is_initialized():
            self.beam.initialize()
        W = np.exp(-1j * self.k * x)
        R = np.zeros([6, 6], dtype='complex')
        R[0] = self.Cu * W
        R[1] = self.Cv * W
        R[2] = self.Cr * W
        k = self.k
        P = self.P
        kuu = self.kuu
        kur = self.kur
        krr = self.krr
        kG = self.kG
        R[3] = (-1j * k * kuu * self.Cu +1j*k*kur*self.Cr)* W
        R[4] = (-kG * (1j * k * self.Cv + self.Cr) -
                 P * self.Cv * 1j * k) * W
        R[5] = (-1j * k * krr * self.Cr +1j*k*kur*self.Cu)* W
        return R


if __name__ == '__main__':
    from beam import Beam
    E = 210e9
    nu = 0.3
    b = 0.1
    sect1 = Section([-0.05,0.05], 210e9, 0.3, 0.1, 7800, 5/6, 'section-1')
    sect2 = Section([-0.05,0,0.05], [210e9,210e9], [0.3,0.3], [0.1,0.1], [7800,7800], [5/6,5/6], 'section-2')
    sect3 = Section([-0.05,0,0.05], [210e9,210e9], [0.3,0.3], [0.1,0.1], [7800,15000], [5/6,5/6], 'section-3')

    beam = Beam(section=sect1)
    beam.append_segment(1.0, sect2, P=0)
    beam.nodes[0].free()
    beam.nodes[1].free()
    def get_segment(freq):
        beam.set_frequency(freq)
        beam.initialize()
        segment = beam.segments[0]
        return segment
    segment = get_segment(100)
    BL = segment.get_response_matrix(0)[3:]
    BR = segment.get_response_matrix(1)[3:]
    B = np.zeros([6,6], dtype='complex')
    B[:3] = BL
    B[3:] = BR
