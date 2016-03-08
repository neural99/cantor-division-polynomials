import itertools

load('CanHyperCurve.sage')

def enumerate_polys(p):
    for a5 in xrange(1, p):
        for a3 in xrange(0, p):
            for a2 in xrange(0, p):
                for a1 in xrange(0, p):
                    for a0 in xrange(0, p):
                        pol = a5 * x^5 + a3 * x^3 + a2 * x^2 + a1 * x + a0
                        yield GF(p)['x'](pol)

def enumerate_polys_all(p, g):
    it = xrange(1, p)
    args = [ it ]
    args.extend( [ xrange(0, p) for i in xrange(0, g) ])
    it = itertools.product(*args)
    for tup in it:
        pol = 0
        for i in xrange(0, g+1):
            pol += tup[i] * x^(g-i)
        yield GF(p)['x'](pol)

p = 3
N = 3
        
s = 0
i = 0
#print p^5-p^3 # Total number of non-singular f
for f in enumerate_polys_all(p,5):
    if f.discriminant() != 0:
        g = SR(str(f))
        i += 1
        H = CanHyperCurve(2, g)
        div_points = len(list(H.division_points(N, p))) - 1
        s += div_points
        print "f  = " + str(g) + "d = " + str(div_points) + " i = " + str(i)
        del H
print "Total. p=3 N=3: " + str(s)
