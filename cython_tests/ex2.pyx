import numpy as np

def f(int n):
    a = np.zeros((n, n, 3), dtype = "float32")
    cdef int i
    cdef int j
    for i in range(n):
        for j in range(n):
            a[i, j] = [1,1,1]

