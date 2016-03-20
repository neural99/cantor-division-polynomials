import itertools

load('CanHyperCurve.sage')

def enumerate_polys_all(q, g):
    finite_field = GF(q, conway=True, prefix='u') 
    units = set(GF(q, conway=True, prefix='u')) - set([0])
    it = units
    args = [ it ]
    args.extend( [ finite_field for i in xrange(0, g) ])
    it = itertools.product(*args)
    for tup in it:
        pol = 0
        for i in xrange(0, g+1):
            pol += tup[i] * x^(g-i)
        yield finite_field['x'](pol)

q = 9
N = 3
        
s = 0
i = 0

for f in enumerate_polys_all(q,5):
    if f.discriminant() != 0:
        g = SR(str(f))
        i += 1
        H = CanHyperCurve(2, g)
        div_points = len(list(H.division_points_naive(N, q))) - 1
        s += div_points / (q*(q-1)^2)
        print "f  = " + str(g) + "d = " + str(div_points) + " i = " + str(i)
        del H
print "Total. q= " + str(q) + " N=" + str(N) + " : " + str(s)
