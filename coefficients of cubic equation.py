#!/usr/bin/env python3


from sympy import *

var('k w J R m kuu kur kG krr P')

a11 = -k**2*kuu + w**2 * m
a13 = k**2*kur-w**2 * R
a22 = -k**2 * kG + w**2 * m - P * k**2
a23 = 1j * k * kG
a31 = a13
a32 = -a23
a33 = -k**2*krr-kG+w**2
a23a32 = a23 * a32

p = a11 * a22 * a33 - a13 * a22 * a31 - a23a32 * a11
p = p.expand()
print('6次项：')
print(p.coeff(k,6).factor())
print()
print('4次项：')
print(p.coeff(k,4).factor(w))
print()
print('2次项：')
print(p.coeff(k,2).factor(w))
print()
print('0次项：')
print(p.coeff(k,0).factor())
print()