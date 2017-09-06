import numpy as np

def f(n):
    a = np.zeros((n, n, 3), dtype = "float32")
    for i in range(n):
        for j in range(n):
            a[i, j] = [1,1,1]

