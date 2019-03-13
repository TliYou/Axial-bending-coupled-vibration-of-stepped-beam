#!/usr/bin/env python3

import numpy as np
from numpy.linalg import solve


def generalized_solve(y, H, f, leastsq=False):
    b_y, A, b_f = _form_A_b_by_sorting(y, H, f)
    x = _solve_generalized(b_y, A, b_f, leastsq)
    f = f.copy()
    y = y.copy()
    idx = 0
    for i in range(len(f)):
        if np.isnan(f[i]):
            f[i] = x[idx]
            idx += 1
    for i in range(len(y)):
        if np.isnan(y[i]):
            y[i] = x[idx]
            idx += 1
    return y, f


def _solve_generalized(b_y, A, b_f, leastsq):
    m = np.zeros(2, dtype='int')
    m[1] = len(b_y)
    m[0] = A.shape[0] - m[1]
    m0, m1 = m
    n = np.zeros(2, dtype='int')
    n[1] = len(b_f)
    n[0] = A.shape[1] - n[1]
    n0, n1 = n
    _A = np.empty([2, 2], dtype='object')
    _A[0, 0] = A[:m0, :n0]
    _A[0, 1] = A[:m0, n0:]
    _A[1, 0] = A[m0:, :n0]
    _A[1, 1] = A[m0:, n0:]
    # 生成左端矩阵B
    B = np.zeros([m0 + m1, m0 + n0]).astype(A.dtype)
    # 第0行
    B[:m0, :n0] = _A[0, 0]
    B[:m0, n0:] = -np.eye(m0)
    # 第1行
    B[m0:, :n0] = _A[1, 0]
    B[m0:, n0:] = 0
    # 生成右端矩阵C
    C = np.zeros([m0 + m1, m1 + n1]).astype(A.dtype)
    # 第0行
    C[:m0, :n1] = -_A[0, 1]
    C[:m0, n1:] = 0
    # 第1行
    C[m0:, :n1] = -_A[1, 1]
    C[m0:, n1:] = np.eye(m1)
    b = np.concatenate([b_f, b_y])
    # 求解
    if B.shape[0] < B.shape[1]:
        import sys
        sys.stderr.write(
            'Number of inputs is smaller than unknown, result may be erroneous.\n')
    if leastsq:
        return solve(B.T.dot(B),
                     B.T.dot(C).dot(b))
    else:
        return solve(B, C.dot(b))


def _form_A_b_by_sorting(y, H, f):
    A = np.zeros_like(H)
    idx_nan = []
    idx_isn = []
    b_y = []
    for i in range(len(y)):
        if np.isnan(y[i]):
            idx_nan.append(i)
        else:
            idx_isn.append(i)
            b_y.append(y[i])
    A[:len(idx_nan), :] = H[idx_nan, :]
    A[len(idx_nan):, :] = H[idx_isn, :]
    b_y = np.array(b_y)
    A2 = np.zeros_like(A)
    idx_nan = []
    idx_isn = []
    b_f = []
    for i in range(len(f)):
        if np.isnan(f[i]):
            idx_nan.append(i)
        else:
            idx_isn.append(i)
            b_f.append(f[i])
    A2[:, :len(idx_nan)] = A[:, idx_nan]
    A2[:, len(idx_nan):] = A[:, idx_isn]
    b_f = np.array(b_f)
    return b_y, A2, b_f


if __name__ == '__main__':
    from numpy.linalg import det
    H = np.array([1, 1, 1, 1, 1, 1,
                  0, 1, 1, 3, 4, 2,
                  0, 0, 1, 1, 2, 1,
                  0, 0, 0, 1, 0, 3,
                  0, 0, 0, 0, 1, 2,
                  2, 0, 0, 0, 0, 1]).reshape(6, 6)
    H = H + H.T
    f0 = np.array([0, 1, 2, 3, 4, 5])
    y0 = H.dot(f0)   # array([25, 39, 21, 26, 26, 31])
    f = np.array([np.nan, 1,   2,     np.nan, np.nan,   np.nan])
    y = np.array([25,   39,   21, 26,  26,  np.nan])
    y, f = generalized_solve(y, H, f, leastsq=True)
    print(y)
    print(f)
