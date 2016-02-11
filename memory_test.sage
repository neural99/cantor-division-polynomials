load('CanHyperCurve.sage')

import wrapper

def random_eq(p):
    x = SR.var('x')
    f = x^5
    for i in xrange(0, 5):
	f += ZZ.random_element(p) * x^i
    return f

x = SR.var('x')
p = 29
W = wrapper.CanHyperCurveWrapper(100)
print get_memory_usage()

for i in xrange(0, 10^4):
    d = W.compute_div_points(3, random_eq(p), p)

del W
print get_memory_usage()

import gc
gc.collect()

print get_memory_usage()
