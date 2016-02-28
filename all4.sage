import itertools
import os
import time

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

def email_error(s):
    msg = "Subject: all4.sage error\n\n"
    msg += s

    os.system('echo "' + msg + '" | ssmtp daniel.lannstrom@gmail.com')

def compute_s(W, N, f, p, weight, i, true_ind, curr_sum):
    while True:
        emailed_error = False
        try:
            g = SR(str(f))
            div_points = SR(W.compute_div_points(N, g, p))
            s = div_points * weight
        except Exception as exc:
            name = type(exc).__name__
            s = "Exception: " + name + "\n" + str(exc)
            print s

            if not emailed_error:
                email_error(s)

            # Try again
            continue
        break

    # Print status line
    print "f = " + str(g) + " d = " + str(div_points) + "w = " + str(weight) +  " i = " + str(i) + " t_ind = " + str(true_ind) + " c_sum = " + str(curr_sum)

    return s

def count1(q, algo='cantor',starting_ind=0, instance_limit=10):
    i = 0
    s = 0

    non_residues = [x for x in GF(q) if kronecker(x,q)==-1]
    r = non_residues[0]

    W = wrapper.CanHyperCurveWrapper(instance_limit, algo)
    
    for true_ind, (f, weight, res) in enumerate(itertools.islice(enumerate_polys(q), starting_ind, None)):
        s += compute_s(W, N, f, q, weight, i, true_ind, s)
        i += 1

        # If res is true, then we need to deal with quadratic non-residues
        if res:
            g = SR(r * f)
            s += compute_s(W, N, g, q, weight, i, true_ind, s)
            i += 1
        
    return s

# Calculate the division points for the representative with index ind
def calculate_specific(ind, algo='cantor',instance_limit=10):
    i = 0
    s = 0

    non_residues = [x for x in xrange(p) if kronecker(x,p)==-1]
    r = non_residues[0]

    W = wrapper.CanHyperCurveWrapper(instance_limit, algo)
    
    for f, weight, res in itertools.islice(enumerate_polys(p), ind, ind + 1):
        print compute_s(W, N, f, p, weight, i)

        if res:
            g = SR(r * f)
            print compute_s(W, N, g, p, weight, i)

# Find the index in enumerate_polys of a given polynomial
def find_specific(poly):
    i = 0
    s = 0

    non_residues = [x for x in xrange(p) if kronecker(x,p)==-1]
    r = non_residues[0]

    for ind, (f, weight, res) in enumerate(enumerate_polys(p)):
        print str(ind) + "  " + str(f)
        if SR(f) == SR(poly):
            return ind
            
        if res:
            g = SR(r*f)
            if g == SR(poly):
                return ind

def email_success(N, q, c, wtime):
    msg = "Subject: Calcuation Done! N = " + str(N) + " q = " + str(q) + " c = " + str(c) + "\n\n"
    msg += "Wall time needed: " + str(wtime) 

    os.system('echo "' + msg + '" | ssmtp daniel.lannstrom@gmail.com')

def calculate_sum(alg, instance_limit, N, qlist):
    last_clock = time.clock()
    #for p in [3,7,11, 13, 17, 19, 23, 29, 31,37,41,43,47,53,59,61,67]:
    for q in qlist:
        x = GF(q)['x'].gen()
        c = count1(q, algo=alg, instance_limit=instance_limit)
        wtime = time.clock() - last_clock
        print "!! Total N = " + str(N) + " q = " + str(q) + ": " + str(c) + ", wtime = " + str(wtime)
        email_success(N, q, c, wtime)

        # Not sure if this helps
        import gc
        gc.collect()

        # Reset timer
        last_clock = time.clock()

def usage():
    print "Usage: sage all4.sage <algo> <instance-limit> <N> <list of q values to calculate>"
    sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        usage()
    try:
        algo = sys.argv[1]
        if algo != 'cantor' and algo != 'naive' :
            print "algo must be either cantor or naive!"
            usage()
        instance_limit = int(sys.argv[2])
        N = int(sys.argv[3])
        qlist = [ int(s) for s in sys.argv[4].split(",") ]
        print "Input: N = " + str(N) + " qlist = " + str(qlist)
        calculate_sum(algo, instance_limit, N, qlist)
    except ValueError:
        usage()
