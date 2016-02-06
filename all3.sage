import itertools

# DOES NOT WORK WITH in char=5 since we need to divide with 5
# Total number of polynomials (including singular):
# f = 3(q-1)q^2 + 2(q-1)q^3
# f(q=7) = 4998
# f(q=11) = 30,250
# f(q=13) = 58,812
# f(q=19) = 266,418

load('CanHyperCurve.sage')

def enumerate_polys_type1(p):
    for a5 in xrange(1, p):
    	for a1 in xrange(0, p):
    	    for a0 in xrange(0, p):
      	          pol = GF(p)['x'](a5 * x^5 + x^3 + x^2 + a1 * x + a0)
                  if pol.discriminant() != 0:
                      yield pol
# a3=0
def enumerate_polys_type2_i(p):
    for a5 in xrange(1, p):
    	for a1 in xrange(0, p):
            for a2 in xrange(1, p):
    	        for a0 in xrange(0, p):
      	            pol = GF(p)['x'](a5 * x^5 + a2 * x^2 + a1 * x + a0)
                    if pol.discriminant() != 0:
                        yield pol

# a2=0
def enumerate_polys_type2_ii(p):
    for a5 in xrange(1, p):
    	for a1 in xrange(0, p):
            for a3 in xrange(1, p):
    	        for a0 in xrange(0, p):
      	            pol = GF(p)['x'](a5 * x^5 + a3 * x^3 + a1 * x + a0)
                    if pol.discriminant() != 0:
                        yield pol
# a2=a3=0
def enumerate_polys_type2_iii(p):
    for a5 in xrange(1, p):
    	for a1 in xrange(0, p):
    	    for a0 in xrange(0, p):
      	        pol = GF(p)['x'](a5 * x^5 + a1 * x + a0)
                if pol.discriminant() != 0:
                    yield pol

def enumerate_polys_type2(p):
    return itertools.chain(enumerate_polys_type2_i(p), enumerate_polys_type2_ii(p), enumerate_polys_type2_iii(p)) 

def count1():
    i = 0
    s = 0

    # Quadratic residue case
    for f in enumerate_polys_type1(p):
        g = SR(str(f))
        H = CanHyperCurve(2, g)
        #div_points = len(list(H.division_points_naive(N,p))) - 1
        div_points = 1
        s += div_points/2
        print "f = " + str(g) + " d = " + str(div_points) + " i = " + str(i)
        i += 1
        del H

    non_residues = [x for x in xrange(p) if kronecker(x,p)==-1]
    r = non_residues[0]
    # Non quadratic residue case
    for f in enumerate_polys_type1(p):
        g = SR(r * f)
        H = CanHyperCurve(2, g)
        #div_points = len(list(H.division_points_naive(N,p))) - 1
        div_points = 1
        s += div_points/2
        print "f = " + str(g) + " d = " + str(div_points) + " i = " + str(i)
        i += 1
        del H

    for f in enumerate_polys_type2(p):
        g = SR(str(f))
        H = CanHyperCurve(2, g)
        #div_points = len(list(H.division_points_naive(N,p))) - 1
        div_points = 1
        s += div_points * (1 / (p-1)^2)
        print "f = " + str(g) + " d = " + str(div_points) + " i = " + str(i)
        i += 1
        del H
        
    return s

#for p in [3, 7, 11, 13, 19]:
for p in [3]:
    N = 3
    x = GF(p)['x'].gen()
    c = count1()
    print "!! Total N = " + str(N) + " p = " + str(p) + ": " + str(c)
