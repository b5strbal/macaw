def addone(x):
    x += 1


def addone_all(vec):
    for i in xrange(len(vec)):
        x[i] += 1


Numpy with C ints:

sage: a = numpy.array(range(10000),dtype=int)
sage: a.dtype
dtype('int64')
sage: type(a[0])
<type 'numpy.int64'>
sage: %timeit addone(a)
100000 loops, best of 3: 10 µs per loop

Numpy with python ints:

sage: a = numpy.array(range(10000),dtype=object)
sage: a.dtype
dtype('O')
sage: type(a[0])
<type 'int'>
sage: %timeit addone(a)
10000 loops, best of 3: 119 µs per loop

List with python int, in place using for loop:

sage: l = range(10000)
sage: type(l[0])
<type 'int'>
sage: %timeit addone_all(l)
1000 loops, best of 3: 1.52 ms per loop

List with Sage Integer, in place using for loop:

sage: l = [Integer(x) for x in range(10000)]
sage: type(l[0])
<type 'sage.rings.integer.Integer'>
sage: %timeit addone_all(l)
1000 loops, best of 3: 1.57 ms per loop

Numpy with Sage Integer:

sage: a = numpy.array([Integer(x) for x in range(10000)],dtype=Integer)
sage: a.dtype
dtype('O')
sage: type(a[0])
<type 'sage.rings.integer.Integer'>
sage: %timeit addone(a)
100 loops, best of 3: 2.87 ms per loop

List with Sage Integer, not in plane using map():

sage: l = [Integer(x) for x in range(10000)]
sage: type(l[0])
<type 'sage.rings.integer.Integer'>
sage: %timeit map(lambda x: x+1, l)
100 loops, best of 3: 3.45 ms per loop

Sage vectors:

sage: v = vector(range(10000))
sage: type(v[0])
<type 'sage.rings.integer.Integer'>
sage: %timeit addone_all(v)
100 loops, best of 3: 3.97 ms per loop

List with python int, not in place using map():

sage: l = range(10000)
sage: type(l[0])
<type 'int'>
sage: %timeit map(lambda x: x+1, l)
100 loops, best of 3: 6.62 ms per loop
