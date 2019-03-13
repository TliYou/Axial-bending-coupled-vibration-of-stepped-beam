#!/usr/bin/env python3

from math import sin, cos, pi

def convert(E11,E22,G12,G23,G31, nu12, varphi):
    S11 = 1 / E11
    S22 = 1 / E22
    Q66 = G12
    S12 = -nu12 / E11
    Q11 = S22 / (S11*S22-S12**2)
    Q22 = S11 / (S11*S22-S12**2)
    Q12 = -S12 / (S11*S22-S12**2)
    Q44 = G23
    Q55 = G31
    varphi = varphi/180*pi
    s = sin(varphi)
    c = cos(varphi)
    Q11b = Q11 * c**4+Q22*s**4+2*(Q12+2*Q66)*s**2*c**2
#    Q22b = Q11 *s**4+Q22*c**4+2*(Q12+2*Q66)*s**2*c**2
#    Q12b = Q12*(s**4+c**4)+(Q11+Q22-4*Q66)*s**2*c**2
    Q55b = Q55 * c**2 + Q44 * s**2
    return Q11b, Q55b


E=convert(144.84e9,9.65e9,4.14e9,3.45e9,4.14e9,0.3,0)

from beam import *
h = 0.0254 / 4
E1 = 144.8e9
L = 0.381
rho = 1389.23
G12 = 4.14e9
G13 = 4.14e9
E0,G0 = convert(E1,9.65e9,G12,3.45e9,G13,0.3,0)
E90,G90 = convert(E1,9.65e9,G12,3.45e9,G13,0.3,90)

sect = Section([-2*h,-h,0,h,2*h], [E0,E90,E90,E0],
               [E0/(2*G0)-1,E90/(2*G90)-1,E90/(2*G90)-1,E0/(2*G0)-1],
               [0.0254]*4, [rho]*4, [5/6]*4, 'section-0')
beam0 = Beam()
beam0.append_segment(L, sect, P=0)
beam0.nodes[0].fixed()
beam0.nodes[1].free()
f = beam0.natural_frequency(0,20000)
nf=f*L**2*(rho/(E1*(h*4)**2))**0.5 * 2 * pi