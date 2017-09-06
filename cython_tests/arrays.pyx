from libcpp.vector cimport vector

def f(int n):
    l = []
    cdef int i
    for i in range(n):
        l.append(i)

def g(int n):
    cdef vector[int] l
    cdef int i
    for i in range(n):
        l.push_back(i)

def h(n):
    l = []
    for i in range(n):
        l.append(i)
