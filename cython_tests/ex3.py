import numpy as np

def f( n):
    #cdef int i
    #cdef int j
    a = [[[0]*3 for i in range(n)] for j in range(n)]
    for i in range(n):
        for j in range(n):
            a[i][j] = [1,1,1]

    print(a[100][100], a[999][999])
