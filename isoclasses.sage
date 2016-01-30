import itertools

load('CanHyperCurve.sage')

p=3
x = GF(p)['x'].gen()

def get_polys():
    prod = itertools.product(list(GF(p)), list(GF(p)), list(GF(p)), list(GF(p)))
    for (a4, a6, a8, a10) in prod:
        pol = x^5 + a4 * x^3 + a6 * x^2  + a8 * x + a10
        if pol.discriminant() != 0:
            yield(pol, a4, a6, a8, a10)

def orbits():
    lst = []
    flat = []
    for (f, a4, a6, a8, a10) in get_polys():
        if f not in flat:
            orbit = set()
            orbit.add(f)
            for alpha in GF(p):
                if alpha != 0:
                    b4 = alpha^4 * a4
                    b6 = alpha^6 * a6
                    b8 = alpha^8 * a8
                    b10 = alpha^10 * a10

                    orbit.add(x^5 + b4 * x^3 + b6 * x^2  + b8 * x + b10)
            lst.append(list(orbit))
            flat.extend(list(orbit))
    return lst

N=3

s = 0
print len(orbits())
i = 0
for o in orbits():
    f = SR(str(o[0]))
    i += 1
    H = CanHyperCurve(2, f)
    div_points = len(list(H.division_points(N, p))) - 1
    s += div_points*len(o)
    print "f = " + str(f) + " div points = " + str(div_points)

    del H
print "Totalt N=3, p=3 : " + str(s)
