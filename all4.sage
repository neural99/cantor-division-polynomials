import itertools

# DOES NOT WORK WITH in char=5 since we need to divide with 5

load('CanHyperCurve.sage')

# a3, a2 != 0
def enumerate_polys_type1(p):
    for a5 in xrange(1, p):
    	for a1 in xrange(0, p):
    	    for a0 in xrange(0, p):
      	          pol = GF(p)['x'](a5 * x^5 + x^3 + x^2 + a1 * x + a0)
                  if pol.discriminant() != 0:
                      yield (pol, 1/2, True)
# a3=0, a2 != 0
def enumerate_polys_type2_i(p):
    # a1 != 0
    for a5 in xrange(1, p):
        for a0 in xrange(0, p):
      	    pol = GF(p)['x'](a5 * x^5 + x^2 + x + a0)
            if pol.discriminant() != 0:
                yield (pol, 1/2, True)
    # a1 = 0
    for a5 in xrange(1, p):
        for a2 in xrange(1, p):
            pol = GF(p)['x'](a5 * x^5 + a2 * x^2 + 1)
            if pol.discriminant() != 0:
                yield (pol, 1/(2*(p-1)), True)

# a2=0, a3 != 0
def enumerate_polys_type2_ii(p):
    # a1, a0 != 0
    for a5 in xrange(1, p):
    	for a3 in xrange(1, p):
      	    pol = GF(p)['x'](a5 * x^5 + a3 * x^3 + x + 1)
            if pol.discriminant() != 0:
                yield (pol, 1/2, True)
    # a1 = 0, a0 != 0
    for a5 in xrange(1, p):
        for a3 in xrange(1, p):
            pol = GF(p)['x'](a5 * x^5 + a3 * x^3 + 1)
            if pol.discriminant() != 0:
                yield (pol, 1/(2*(p-1)), True)
    # a1 != 0, a0 = 0
    for a5 in xrange(1, p):
        for a3 in xrange(1, p):
            pol = GF(p)['x'](a5 * x^5 + a3 * x^3 + x)
            if pol.discriminant() != 0:
                yield (pol, 1/(p-1), False)
    # a1=a0=0 is always singular
    
# a2=a3=0
def enumerate_polys_type2_iii(p):
    # a1, a0 != 0
    for a5 in xrange(1, p):
      	pol = GF(p)['x'](a5 * x^5 + x + 1)
        if pol.discriminant() != 0:
            yield (pol, 1/2, True)
    # a1 = 0, a0 != 0
    for a5 in xrange(1, p):
        pol = GF(p)['x'](a5 * x^5 + 1)
        if pol.discriminant() != 0:
            yield (pol, 1/(2*(p-1)), True)
    # a0 = 0, a1 != 0
    for a5 in xrange(1, p):
        pol = GF(p)['x'](a5 * x^5 + x)
        if pol.discriminant() != 0:
            yield(pol, 1/(p-1), False)
    # a0=a1=0 is alsways singular
                         
    
def enumerate_polys(p):
    return itertools.chain(enumerate_polys_type1(p),
                           enumerate_polys_type2_i(p),
                           enumerate_polys_type2_ii(p),
                           enumerate_polys_type2_iii(p)) 

def count1():
    i = 0
    s = 0

    non_residues = [x for x in xrange(p) if kronecker(x,p)==-1]
    r = non_residues[0]
    
    for f, weight, res in enumerate_polys(p):
        g = SR(str(f))
        H = CanHyperCurve(2, g)
        div_points = len(list(H.division_points_naive(N,p))) - 1
        s += div_points * weight
        print "f = " + str(g) + " d = " + str(div_points) + " i = " + str(i)
        i += 1
        del H

        # If res is true, then we need to deal with quadratic non-residues
        if res:
            g = SR(r * f)
            H = CanHyperCurve(2, g)
            div_points = len(list(H.division_points_naive(N,p))) - 1
            s += div_points * weight
            print "f = " + str(g) + " d = " + str(div_points) + " i = " + str(i)
            i += 1
            del H
        
    return s

for p in [3, 7, 11, 13, 17, 19]:
#for p in [17]:
    N = 3
    x = GF(p)['x'].gen()
    c = count1()
    print "!! Total N = " + str(N) + " p = " + str(p) + ": " + str(c)
