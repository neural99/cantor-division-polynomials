import itertools
import os
import time

EMAIL = 'change@me.com'

# DOES NOT WORK in char=5 since we need to divide with 5

load('CanHyperCurve.sage')

import wrapper

# a3, a2 != 0
def enumerate_polys_type1(q):
    print "A"
    finite_field = GF(q, conway=True, prefix='u') 
    units = set(GF(q, conway=True, prefix='u')) - set([0])
    for a5 in units:
    	for a1 in finite_field:
    	    for a0 in finite_field:
      	          pol = finite_field['x'](a5 * x^5 + x^3 + x^2 + a1 * x + a0)
                  if pol.discriminant() != 0:
                      yield (pol, 1/2, True)
# a3=0, a2 != 0
def enumerate_polys_type2_i(q):
    print "B.1"
    finite_field = GF(q, conway=True, prefix='u') 
    units = set(GF(q, conway=True, prefix='u')) - set([0])

    # a1 != 0
    for a5 in units:
        for a0 in finite_field:
      	    pol = finite_field['x'](a5 * x^5 + x^2 + x + a0)
            if pol.discriminant() != 0:
                yield (pol, 1/2, True)
    print "B.2"
    # a1 = 0
    for a5 in units:
        for a2 in units:
            pol = finite_field['x'](a5 * x^5 + a2 * x^2 + 1)
            if pol.discriminant() != 0:
                yield (pol, 1/(2*(q-1)), True)

# a2=0, a3 != 0
def enumerate_polys_type2_ii(q):
    print "C.1"
    finite_field = GF(q, conway=True, prefix='u') 
    units = set(GF(q, conway=True, prefix='u')) - set([0])

    # a1, a0 != 0
    for a5 in units:
    	for a3 in units:
      	    pol = finite_field['x'](a5 * x^5 + a3 * x^3 + x + 1)
            if pol.discriminant() != 0:
                yield (pol, 1/2, True)
    print "C.2"
    # a1 = 0, a0 != 0
    for a5 in units:
        for a3 in units:
            pol = finite_field['x'](a5 * x^5 + a3 * x^3 + 1)
            if pol.discriminant() != 0:
                yield (pol, 1/(2*(q-1)), True)
    print "C.3"
    # a1 != 0, a0 = 0
    for a5 in units:
        for a3 in units:
            pol = finite_field['x'](a5 * x^5 + a3 * x^3 + x)
            if pol.discriminant() != 0:
                yield (pol, 1/(q-1), False)
    # a1=a0=0 is always singular
    
# a2=a3=0
def enumerate_polys_type2_iii(q):
    print "D.1"
    finite_field = GF(q, conway=True, prefix='u') 
    units = set(GF(q, conway=True, prefix='u')) - set([0])

    # a1, a0 != 0
    for a5 in units:
      	pol = finite_field['x'](a5 * x^5 + x + 1)
        if pol.discriminant() != 0:
            yield (pol, 1/2, True)
    print "D.2"
    # a1 = 0, a0 != 0
    for a5 in units:
        pol = finite_field['x'](a5 * x^5 + 1)
        if pol.discriminant() != 0:
            yield (pol, 1/(2*(q-1)), True)
    print "D.3"
    # a0 = 0, a1 != 0
    for a5 in units:
        pol = finite_field['x'](a5 * x^5 + x)
        if pol.discriminant() != 0:
            yield(pol, 1/(q-1), False)
    # a0=a1=0 is alsways singular
                         
def enumerate_polys(q):
    return itertools.chain(enumerate_polys_type1(q),
                           enumerate_polys_type2_i(q),
                           enumerate_polys_type2_ii(q),
                           enumerate_polys_type2_iii(q)) 

def email_error(s):
    if email_flag == 'true':
        msg = "Subject: all4.sage error\n\n"
        msg += s
        os.system('echo "' + msg + '" | ssmtp ' + EMAIL)

def compute_s(W, N, f, q, weight, i, true_ind, curr_sum, starting_index):
    while True:
        emailed_error = False
        try:
            g = SR(str(f))
            div_points = SR(W.compute_div_points(N, g, q))
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
    print "f = " + str(g) + " d = " + str(div_points) + "w = " + str(weight) +  " i = " + str(i) + " t_ind = " + str(true_ind+starting_index) + " c_sum = " + str(curr_sum)

    return s

def count1(N, q, algo='cantor',starting_ind=0, instance_limit=10, c_sum=0):
    i = 0
    s = c_sum 

    finite_field = GF(q, conway=True, prefix='u')
    units = set(finite_field) - set([0])
    generator = finite_field.primitive_element()
    non_residues = [ generator^(2*k+1) for k in xrange(1, (q-1)/2 + 1) ]
    r = non_residues[0]

    W = wrapper.CanHyperCurveWrapper(instance_limit, algo)
    
    for true_ind, (f, weight, res) in enumerate(itertools.islice(enumerate_polys(q), starting_ind, None)):
        s += compute_s(W, N, f, q, weight, i, true_ind, s, starting_ind)
        i += 1

        # If res is true, then we need to deal with quadratic non-residues
        if res:
            g = SR(r * f)
            s += compute_s(W, N, g, q, weight, i, true_ind, s, starting_ind)
            i += 1
        
    return s

# Calculate the division points for the representative with index ind
def calculate_specific(N, q, ind, algo='cantor',instance_limit=10):
    i = 0
    s = 0

    finite_field = GF(q, conway=True, prefix='u')
    units = set(finite_field) - set([0])
    generator = finite_field.primitive_element()
    non_residues = [ generator^(2*k+1) for k in xrange(1, (q-1)/2 + 1) ]
    r = non_residues[0]

    W = wrapper.CanHyperCurveWrapper(instance_limit, algo)
    
    for f, weight, res in itertools.islice(enumerate_polys(p), ind, ind + 1):
        print compute_s(W, N, f, q, weight, i)

        if res:
            g = SR(r * f)
            print compute_s(W, N, g, q, weight, i)

# Find the index in enumerate_polys of a given polynomial
def find_specific(poly):
    i = 0
    s = 0

    finite_field = GF(q, conway=True, prefix='u')
    units = set(finite_field) - set([0])
    generator = finite_field.primitive_element()
    non_residues = [ generator^(2*k+1) for k in xrange(1, (q-1)/2 + 1) ]
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
    if email_flag == 'true':
        msg = "Subject: Calcuation Done! N = " + str(N) + " q = " + str(q) + " c = " + str(c) + "\n\n"
        msg += "Wall time needed: " + str(wtime) 
        os.system('echo "' + msg + '" | ssmtp ' + EMAIL)

def calculate_sum(alg, instance_limit, N, qlist, index=0, c_sum=0):
    last_clock = time.clock()
    for q in qlist:
        finite_field = GF(q, conway=True, prefix='u')
        x = finite_field['x'].gen()
        c = count1(N, q, algo=alg, instance_limit=instance_limit, starting_ind=index, c_sum=c_sum)
        wtime = time.clock() - last_clock
        print "!! Total N = " + str(N) + " q = " + str(q) + ": " + str(c) + ", wtime = " + str(wtime)
        email_success(N, q, c, wtime)

        # Not sure if this helps
        import gc
        gc.collect()

        # Reset timer
        last_clock = time.clock()

def usage():
    print "Usage: sage all4.sage <email-result?> <algo> <instance-limit> <N> <list of q values to calculate> <starting index> <starting sum>"
    sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 8:
        usage()
    try:
        email_flag = sys.argv[1]
        algo = sys.argv[2]
        if algo != 'cantor' and algo != 'naive' :
            print "algo must be either cantor or naive!"
            usage()
        instance_limit = int(sys.argv[3])
        N = int(sys.argv[4])
        qlist = [ int(s) for s in sys.argv[5].split(",") ]
        t_ind = int(sys.argv[6])
        s_sum = int(sys.argv[7])
        print "Input: N = " + str(N) + " qlist = " + str(qlist) + " starting t_ind = " + str(t_ind) + " starting sum = " + str(s_sum)
    except ValueError as e:
        usage()
    calculate_sum(algo, instance_limit, N, qlist, index=t_ind, c_sum=s_sum)
