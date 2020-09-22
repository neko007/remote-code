import numpy as np
import matplotlib.pyplot as plt
import time

def magic(n):
    '''
    生成n阶魔方矩阵
    '''
    row,col = 0, n//2
    magic = []
    for i in range(n):
        magic.append([0]*n)
    magic[row][col] = 1
    for i in range(2, n*n+1):
        r, l = (row-1+n)%n, (col+1)%n
        if(magic[r][l]==0):
            row, col=r, l
        else:
            row = (row+1)%n
        magic[row][col] = i
    marray = np.array(magic)
    return marray

def triangulor(A, b):
    # 先转化上三角
    N = A.shape[0]
    phi = np.full_like(b, np.nan)
    for i in range(1, N):
        e = A[i][i-1] / A[i-1][i-1]
        A[i][i-1] = 0
        A[i] = A[i] - e * A[i-1]
        b[i] = b[i] - e * b[i-1]

    # 在转化为对角矩阵
    phi[-1] = b[-1] / A[-1][-1]
    for i in range(N-2, -1, -1):
        phi[i] = (b[i] - A[i][i+1]*phi[i+1]) / A[i][i]

    return phi

def initialize(N):
    phi_0, phi_20 = 2, 3
    A = np.zeros((N, N))
    b = np.zeros(N)

    # 赋初值
    A[0, 0], A[-1, -1] = 1, 1
    b[0], b[-1] = phi_0, phi_20
    for i in range(1, N-1):
        A[i, i-1:i+2] = [1, -2, 1]
    return A, b

def inv(A, b):
    A, b = A.copy(), b.copy()
    t1 = time.time()
    phi = np.linalg.inv(A) @ b
    t2 = time.time()
    return phi, t2-t1

def tri(A, b):
    A, b = A.copy(), b.copy()
    t1 = time.time()
    phi = triangulor(A, b)
    t2 = time.time()
    return phi, t2-t1



if __name__ == '__main__':
    # phi_{i+1} - 2phi_{i} + phi_{i-1} = 0
    N = 21
    A, b = initialize(N)
    phi1, t1 = inv(A, b)
    phi2, t2 = tri(A, b)

    print(t1, t2)

    nx_max = 1000
    t = np.full((3, nx_max), np.nan)
    for N in range(2, nx_max):
        print(N)
        A, b = initialize(N+1)
        _, t[0, N] = inv(A, b)
        _, t[1, N] = tri(A, b)

        A = magic(N)
        b = np.arange(1, N+1)
        _, t[2, N] = inv(A, b)




