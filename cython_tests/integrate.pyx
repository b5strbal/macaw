cdef extern from "math.h":
    double sin(double)

cdef double f(double x) except *:
    return sin(x**2)

def integrate_f(double a, double b, int N):
    cdef double s = 0
    cdef double dx
    cdef int i
    cdef list l = []
    dx = (b-a)/N
    for i in range(N):
        s += f(a+i*dx)
        l.append(s)
    return s * dx
