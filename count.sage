load('CanHyperCurve.sage')

p = 7
x = GF(p)['x'].gen()
N = 3

def enumerate_polys(p):
    for a5 in xrange(1, p):
        for a3 in xrange(0, p):
            for a2 in xrange(0, p):
                for a1 in xrange(0, p):
                    for a0 in xrange(0, p):
                        pol = a5 * x^5 + a3 * x^3 + a2 * x^2 + a1 * x + a0
                        yield GF(p)['x'](pol)

s = 0
for f in enumerate_polys(p):
    if f.discriminant() == 0:
        continue
    g = SR(str(f))
    print "f = " + str(g)
    H = CanHyperCurve(2, g)
    div_points = len(list(H.division_points(N, p))) - 1
    print "# div points = " + str(div_points)
    s += div_points/(p-1)^2
print "Total. N=3, p=7. " + s
