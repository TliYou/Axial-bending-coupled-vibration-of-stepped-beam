"""beam_property.py

The section for a multi-layered beam.
"""

from math import pi
import numpy as np

__all__ = ['Section']


class ReadOnly:
    """Calss for inherit.
    The variables will be read-only when initialized."""

    def __init__(self):
        super().__setattr__('_initialized', False)

    def _initialize(self):
        self._initialized = True

    def __setattr__(self, name, value):
        if self._initialized:
            raise AttributeError('Cannot set value after initilized.')
        else:
            super().__setattr__(name, value)


class Section(ReadOnly):
    """Define a Section of a beam.

    Paremeters
    ----------
    y: iterable
        Coordinate of each layer, length should be number of layers + 1

    E: float or iterable
        Young's modulus for each layer.

    nu: float or iterable between 0 and 0.5
        Possion's ratio for each layer.

    b: float or iterable
        Width of the corss-section for each layer.

    rho: float or iterable
        Density of the material. If the beam is a multi-layered beam, rho
        should be iterable containing the density of each layer. (from bottom
        to top)

    kappa: float or iterable
        Shear coefficient for each layer.

    name: string, optional
        Name of the Section.
    """

    def __init__(self, y, E, nu, b, rho, kappa, name=''):
        super().__init__()
        y = y if np.iterable(y) else [y]
        self.y = np.array(y)
        n = len(self.y) - 1
        E = E if np.iterable(E) else [E]
        self.E = np.array(E)
        nu = nu if np.iterable(nu) else [nu]
        self.nu = np.array(nu)
        b = b if np.iterable(b) else [b]
        self.b = np.array(b)
        rho = rho if np.iterable(rho) else [rho]
        self.rho = np.array(rho)
        kappa = kappa if np.iterable(kappa) else [kappa]
        self.kappa = np.array(kappa)
        if n != len(self.E):
            raise ValueError('length of E should be number of layers')
        if n != len(self.nu):
            raise ValueError('length of nu should be number of layers')
        if n != len(self.b):
            raise ValueError('length of b should be number of layers')
        if n != len(self.rho):
            raise ValueError('length of rho should be number of layers')
        if n != len(self.kappa):
            raise ValueError('length of kappa should be number of layers')
        self._calculate_attrs()
        self.name = name
        self._initialize()

    def _calculate_attrs(self):
        y = self.y
        E = self.E
        nu = self.nu
        b = self.b
        rho = self.rho
        kappa = self.kappa
        y = self.y
        h = np.diff(y)

        A = b * h
        S = 0.5 * b * (y[1:]**2 - y[:-1]**2)
        I = 1 / 3 * b * (y[1:]**3 - y[:-1]**3)
        G = E / 2 / (1 + nu)

        self.A = A
        self.S = S
        self.I = I
        self.G = G

        self.J = sum(rho * I)
        self.R = sum(rho * S)
        self.m = sum(rho * A)
        self.kuu = sum(E * A)
        self.kur = sum(E * S)
        self.kG = sum(kappa * G * A)
        self.krr = sum(E * I)

    def __repr__(self):
        return 'Section "%s"' % self.name


if __name__ == '__main__':
    E = 210e9
    nu = 0.3
    b = 0.1
    sect1 = Section([-0.05,0.05], 210e9, 0.3, 0.1, 7800, 5/6, 'section-1')
    sect2 = Section([-0.05,0,0.05], [210e9,210e9], [0.3,0.3], [0.1,0.1], [7800,7800], [5/6,5/6], 'section-2')
    sect3 = Section([-0.05,0,0.05], [210e9,210e9], [0.3,0.3], [0.1,0.1], [7800,15000], [5/6,5/6], 'section-3')
#    sect4 = Section([-0.05,0,0.05], [210e9,210e9], [0.3,0.3], [0.1,0.1], [7800,15000], [5/6,5/6], 'section-3')
