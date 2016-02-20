import itertools

# DOES NOT WORK WITH in char=5 since we need to divide with 5

load('CanHyperCurve.sage')

import wrapper

# a3, a2 != 0
def enumerate_polys_type1(p):
    print "A"
    for a5 in xrange(1, p):
    	for a1 in xrange(0, p):
    	    for a0 in xrange(0, p):
      	          pol = GF(p)['x'](a5 * x^5 + x^3 + x^2 + a1 * x + a0)
                  if pol.discriminant() != 0:
                      yield (pol, 1/2, True)
# a3=0, a2 != 0
def enumerate_polys_type2_i(p):
    print "B.1"
    # a1 != 0
    for a5 in xrange(1, p):
        for a0 in xrange(0, p):
      	    pol = GF(p)['x'](a5 * x^5 + x^2 + x + a0)
            if pol.discriminant() != 0:
                yield (pol, 1/2, True)
    print "B.2"
    # a1 = 0
    for a5 in xrange(1, p):
        for a2 in xrange(1, p):
            pol = GF(p)['x'](a5 * x^5 + a2 * x^2 + 1)
            if pol.discriminant() != 0:
                yield (pol, 1/(2*(p-1)), True)

# a2=0, a3 != 0
def enumerate_polys_type2_ii(p):
    print "C.1"
    # a1, a0 != 0
    for a5 in xrange(1, p):
    	for a3 in xrange(1, p):
      	    pol = GF(p)['x'](a5 * x^5 + a3 * x^3 + x + 1)
            if pol.discriminant() != 0:
                yield (pol, 1/2, True)
    print "C.2"
    # a1 = 0, a0 != 0
    for a5 in xrange(1, p):
        for a3 in xrange(1, p):
            pol = GF(p)['x'](a5 * x^5 + a3 * x^3 + 1)
            if pol.discriminant() != 0:
                yield (pol, 1/(2*(p-1)), True)
    print "C.3"
    # a1 != 0, a0 = 0
    for a5 in xrange(1, p):
        for a3 in xrange(1, p):
            pol = GF(p)['x'](a5 * x^5 + a3 * x^3 + x)
            if pol.discriminant() != 0:
                yield (pol, 1/(p-1), False)
    # a1=a0=0 is always singular
    
# a2=a3=0
def enumerate_polys_type2_iii(p):
    print "D.1"
    # a1, a0 != 0
    for a5 in xrange(1, p):
      	pol = GF(p)['x'](a5 * x^5 + x + 1)
        if pol.discriminant() != 0:
            yield (pol, 1/2, True)
    print "D.2"
    # a1 = 0, a0 != 0
    for a5 in xrange(1, p):
        pol = GF(p)['x'](a5 * x^5 + 1)
        if pol.discriminant() != 0:
            yield (pol, 1/(2*(p-1)), True)
    print "D.3"
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

def compute_s(W, N, f, p, weight, i):
    g = SR(str(f))
    div_points = SR(W.compute_div_points(N, g, p))
    s = div_points * weight
    print "f = " + str(g) + " d = " + str(div_points) + "w = " + str(weight) +  " i = " + str(i)
    return s

def count1():
    i = 0
    s = 0

    non_residues = [x for x in xrange(p) if kronecker(x,p)==-1]
    r = non_residues[0]

    W = wrapper.CanHyperCurveWrapper(10)
    
    for f, weight, res in enumerate_polys(p):
        s += compute_s(W, N, f, p, weight, i)
        i += 1

        # If res is true, then we need to deal with quadratic non-residues
        if res:
            g = SR(r * f)
            s += compute_s(W, N, g, p, weight, i)
            i += 1
        
    return s

#for p in [3,7,11, 13, 17, 19, 23, 29, 31,37,41,43,47,53,59,61,67]:
for p in [31,37,41,43,47,53,59,61,67]:
    N = 3
    x = GF(p)['x'].gen()
    c = count1()
    print "!! Total N = " + str(N) + " p = " + str(p) + ": " + str(c)

    import gc
    gc.collect()
