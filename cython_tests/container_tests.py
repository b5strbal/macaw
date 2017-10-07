def addone(x):
    x += 1


def addone_all(vec):
    for i in xrange(len(vec)):
        x[i] += 1


Numpy with C ints:

>>> a = numpy.array(range(10000),dtype=int)
>>> a.dtype
dtype('int64')
>>> type(a[0])
<type 'numpy.int64'>
>>> %timeit addone(a)
100000 loops, best of 3: 10 µs per loop

Numpy with python ints:

>>> a = numpy.array(range(10000),dtype=object)
>>> a.dtype
dtype('O')
>>> type(a[0])
<type 'int'>
>>> %timeit addone(a)
10000 loops, best of 3: 119 µs per loop

List with python int, in place using for loop:

>>> l = range(10000)
>>> type(l[0])
<type 'int'>
>>> %timeit addone_all(l)
1000 loops, best of 3: 1.52 ms per loop

List with Sage Integer, in place using for loop:

>>> l = [Integer(x) for x in range(10000)]
>>> type(l[0])
<type 'sage.rings.integer.Integer'>
>>> %timeit addone_all(l)
1000 loops, best of 3: 1.57 ms per loop

Numpy with Sage Integer:

>>> a = numpy.array([Integer(x) for x in range(10000)],dtype=Integer)
>>> a.dtype
dtype('O')
>>> type(a[0])
<type 'sage.rings.integer.Integer'>
>>> %timeit addone(a)
100 loops, best of 3: 2.87 ms per loop

List with Sage Integer, not in plane using map():

>>> l = [Integer(x) for x in range(10000)]
>>> type(l[0])
<type 'sage.rings.integer.Integer'>
>>> %timeit map(lambda x: x+1, l)
100 loops, best of 3: 3.45 ms per loop

Sage vectors:

>>> v = vector(range(10000))
>>> type(v[0])
<type 'sage.rings.integer.Integer'>
>>> %timeit addone_all(v)
100 loops, best of 3: 3.97 ms per loop

List with python int, not in place using map():

>>> l = range(10000)
>>> type(l[0])
<type 'int'>
>>> %timeit map(lambda x: x+1, l)
100 loops, best of 3: 6.62 ms per loop
