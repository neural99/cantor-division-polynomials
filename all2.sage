import itertools

# DOES NOT WORK FOR p=5 since we need to divide by 5!!

load('CanHyperCurve.sage')

def enumerate_polys(q):
    finite_field = GF(q, conway=True, prefix='u') 
    units = set(GF(q, conway=True, prefix='u')) - set([0])
    for a5 in units:
        for a3 in finite_field:
            for a2 in finite_field:
                for a1 in finite_field:
                    for a0 in finite_field:
                        pol = finite_field['x'](a5 * x^5 + a3 * x^3 + a2 * x^2 + a1 * x + a0)
                        if pol.discriminant() != 0:
                            yield pol


N = 3
q = 9

finite_field = GF(q, conway=True, prefix='u') 
x = finite_field['x'].gen()

s = 0
i = 0
for f in enumerate_polys(q):
    g = SR(str(f))
    i += 1
    H = CanHyperCurve(2, g)
    div_points = len(list(H.division_points_naive(N, q))) - 1
    #div_points = 1
    s += div_points/(q-1)^2
    print "f  = " + str(g) + "d = " + str(div_points) + " i = " + str(i)
    del H
print "Total. N=" + str(N) + " q = " + str(q) + ": " + str(s)
