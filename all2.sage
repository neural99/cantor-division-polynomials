import itertools

# DOES NOT WORK FOR p=5 since we need to divide by 5!!

load('CanHyperCurve.sage')

def enumerate_polys(p):
    for a5 in xrange(1, p):
        for a3 in xrange(0, p):
            for a2 in xrange(0, p):
                for a1 in xrange(0, p):
                    for a0 in xrange(0, p):
                        pol = GF(p)['x'](a5 * x^5 + a3 * x^3 + a2 * x^2 + a1 * x + a0)
                        if pol.discriminant() != 0:
                            yield pol


N = 3
p = 3

x = GF(p)['x'].gen()

s = 0
i = 0
print p^5-p^3 # Total number of non-singular f
for f in enumerate_polys(p):
    g = SR(str(f))
    i += 1
    H = CanHyperCurve(2, g)
    #div_points = len(list(H.division_points_naive(N, p))) - 1
    div_points = 1
    s += div_points
    print "f  = " + str(g) + "d = " + str(div_points) + " i = " + str(i)
    del H
print "Total. p=3 N=3: " + str(s)
