import itertools
import random

# Test Cantor's division point algorithm on some random curves of genus 2
# over different prime fields F_p with p < 100
def test_randoms():
    def test_find_V(H, p):
        for U in enumerate_monic_polys(p, 2):
            if sorted(list(H._naive_enum_mumford_V(U,p))) != sorted(H._genus_2_enum_mumford_V(U, p)):
                print "test_find_V failed: p = " + str(p) + " U = " + str(U)
        print "test_V OK"

    def random_eq(p):
        x = SR.var('x')
        f = x^5
        for i in xrange(0, 5):
            f += ZZ.random_element(p) * x^i
        return f

    levels = [ 3 ]

    P = Primes()
    p = P.next(3)
    while p < 30:
        found = False
        f = None
        while not found:
            # Test that it is non singular 
            try:
                f = random_eq(p)
                HyperellipticCurve(GF(p)['x'](f))
                found = True
            except ValueError:
                continue
        level = random.choice(levels)
        print "p = " + str(p) + " n = " + str(level) + " f = " + str(f)
        H = CanHyperCurve(2, f)
#        test_find_V(H, p)
        lst1 = sorted(map(lambda k: str(k), list(set(H.division_points(level, p)))))
        lst2 = sorted(map(lambda k: str(k), list(set(H.division_points_naive(level, p)))))
        print "Cantor: " + str(lst1)
        print "Naive: " + str(lst2)
        if lst1 != lst2:
            print "failed: p=" + str(p) + " f = " + str(f)
        else:
            print "OK"
        p = P.next(p)

def cantor_benchmark():
    # Make them visible below
    global H
    global p
    import sage.misc.sage_timeit
    x = SR.var('x')
    H = CanHyperCurve(2, x^5+1)
    lst = []
    P = Primes()
    p = P.next(3)
    while p < 100:
        try:
            t = sage.misc.sage_timeit.sage_timeit('list(H.division_points(5,p))', globals_dict=globals(), number=1, repeat=1, seconds=True)
            print p, t
            lst.append([p, t])
        except ValueError:
            print "Skipping singular curve"
        p = P.next(p)
    return lst

def naive_benchmark():
    # Make them visible below
    global H
    global p
    import sage.misc.sage_timeit
    x = SR.var('x')
    H = CanHyperCurve(2, x^5+1)
    lst = []
    P = Primes()
    p = P.next(3)
    while p < 100:
        try:
            t = sage.misc.sage_timeit.sage_timeit('list(H.division_points_naive(5,p))', globals_dict=globals(), number=1, repeat=1, seconds=True)
            print p, t
            lst.append([p, t])
        except ValueError:
            print "Skipping singular curve"
        p = P.next(p)
    return lst

# Write symmetric polynoimal in elementary symmetric polynomials.
# Adapted from https://groups.google.com/forum/#!msg/sage-support/mOhavaXWIjQ/bSbkK6_2BC0J
def symmetric_poly_on_elementary_basis(poly):
    R.<x0,x1,s1,s2> = PolynomialRing(QQ,order = TermOrder('degrevlex',2)+TermOrder('degrevlex',2))
    sym = []
    xvars = [x0,x1]
    for i in range(1,3):
        temp = 0
        for p in itertools.combinations([0,1],i):
            temp = temp + prod([xvars[p[j]] for j in range(i)])
        sym.append(sage_eval('s' + str(i), locals=locals()) - temp)
    symid = R.ideal(sym + [R(poly)]).elimination_ideal([x0,x1])
    return symid.gens()[0]

# Extract the coeff of kth power in the coeffients list outputed by coefficients()
def _get_series_coeff(coeffs, k):
    for c in coeffs:
        if c[1] == k: return c[0]
    return 0

# Reduce the expression from a symbolic polynomial in x over Z to a poly over F_p
def _reduce_poly(exp, p):
    x = SR.var('x')
    R = GF(p)['t']
    t = R.gen()
    poly = R(exp.substitute(x=t))
    return poly

# We have r*(x,y)=(delta(X), epsilon(X)) where the right hand side is semi-Mumford meaning that
# delta(X) is not monic. This function gives a standard Mumford representation.
def _get_monic_cantor_pol(delta, epsilon):
    leading_coeff =  delta.list()[-1]
    delta =  delta / leading_coeff
    # If we have switched the sign of delta, switch the sign of epsilon too
    if leading_coeff < 0:
        epsilon = -epsilon
    return delta, epsilon

# Enumerate monic polys in F_p[x] of degree <= g
def enumerate_monic_polys(p, g):
    lst = []
    for d in xrange(0, g + 1): # 0 <= d <= g
        d_polys = _enumerate_polys_rec(p, d, [0], 0)
        lst.extend(d_polys)
    return lst

def enumerate_polys(p, g):
    yield 0
    for poly in enumerate_monic_polys(p, g):
        for n in xrange(1, p):
            yield n * poly

def _enumerate_polys_rec(p, d, lst, i):
    if i > d: return lst
    elif i == d:
        lst2 = map(lambda p: p + x^i, lst)
        return _enumerate_polys_rec(p, d, lst2, i + 1)
    else:
        lst2 = []
        for poly in lst:
            lst2.extend([ poly + x^i * coeff for coeff in xrange(0, p) ])
        return _enumerate_polys_rec(p, d, lst2, i + 1)

class CanHyperCurve(sage.structure.sage_object.SageObject):
    # g = genus of curve
    # G = equation i.e. the curve is given by y^2=G(x) where G is a polynomial of degree 2g+1
    def __init__(self, g, G):
        self.g = g
        self.G = G

        x = SR.var('x')
        z = SR.var('z')
        self.F = G.substitute(x=(x-z))
        # Note the (-1)^(g+1) is there to make the sign of the leading coefficient of the div polys positive
        self.S = (-1)^(g+1) * self.F.sqrt()

        if g == 2:
            self._enum_mumford_V = self._genus_2_enum_mumford_V
        else:
            self._enum_mumford_V = self._naive_enum_mumford_V

    # The following function computes the cantor div polys (psi) from the definition given in Cantors article. I.e. using
    # Pade approximants.

    @cached_method
    def _get_S_cached(self, l):
        z = SR.var('z')
        return self.S.series(z, l).coefficients(z)

    # Get the Hankel Matrix corresponding with the series of S
    # I.e. H_m,n(S)
    def _compute_hankel_matrix(self, m, n):
        M = matrix(SR, n, n)
        coeffs = self._get_S_cached(m+n)
        for i in xrange(0, n):
            for j in xrange(0, n):
                x1 = m-n+1+i+j
                M[j,i] = _get_series_coeff(coeffs, x1)
        return M

    # Compute f_r, i.e. the unnormalized division polynomial
    def compute_F(self, r):
        if r < self.g:
            return 0

        mr = floor((r+1+self.g)/2)
        nr = floor((r+1-self.g-1)/2)
        M = self._compute_hankel_matrix(mr, nr)
        return M.determinant()

    # Compute psi_r, i.e. the normalized division polynomial directly from the determinant
    def compute_psi_from_determinant(self, r):
        vr = (r^2 - r -self.g^2 + self.g)/2
        psi = self.compute_F(r) * 2^vr * self.G^(vr/2)
        return psi

    def _compute_B(self, r):
        if r <= self.g: return 0

        m = floor((r+self.g)/2)
        n = floor((r-self.g-1)/2)

        M = matrix(SR, n + 1, n + 1)
        coeffs = self._get_S_cached(m + n + 1)
        for i in xrange(0, n + 1):
            for j in xrange(0, n):
                x1 = m-n+1+i+j
                M[j,i] = _get_series_coeff(coeffs, x1)
        # Fill out last row
        z = SR.var('z')
        for i in xrange(0, n + 1):
            M[n, i] = z^(n-i)
        return M.determinant()

    def _compute_S_poly(self, coeffs, k):
        z = SR.var('z')
        poly = 0
        for i in xrange(0, k + 1):
            poly = poly + z^i * _get_series_coeff(coeffs, i)
        return poly

    def _compute_A(self, r):
        z = SR.var('z')
        if r <= self.g: return -z^r

        m = floor((r+self.g)/2)
        n = floor((r-self.g-1)/2)

        M = matrix(SR, n + 1, n + 1)
        coeffs = self._get_S_cached(m + n + 1)
        for i in xrange(0, n + 1):
            for j in xrange(0, n):
                x1 = m-n+1+i+j
                M[j,i] = _get_series_coeff(coeffs, x1)
        # Fill out last row
        for i in xrange(0, n + 1):
            M[n, i] = z^(n-i) * self._compute_S_poly(coeffs, m - n + i)
        return M.determinant()

    def _compute_D(self, r):
        z = SR.var('z')
        D = -(self._compute_A(r)^2 - self._compute_B(r)^2 * self.F) / z^r
        return D

    def _compute_delta(self, r):
        z = SR.var('z')
        D = self._compute_D(r)
        # Normalize D
        vr = (r^2 - r -self.g^2 + self.g)/2
        delta =  D.substitute(z=4*self.G*z) * 2^(2*vr) * self.G^(vr)
        return delta.factor()

    def _compute_epsilon(self, r):
        y = SR.var('y')
        z = SR.var('z')
        num = y * z * (self.compute_psi(r-1)^2 * self._compute_delta(r+1) - self.compute_psi(r+1)^2 * self._compute_delta(r-1))
        den = self.compute_psi(r-1) * self.compute_psi(r)^2 * self.compute_psi(r+1)
        R = ZZ['x,y']
        F = PolynomialRing(R, 'z')
        delta = F(self._compute_delta(r))
        quot = F(num) / R(den)
        q, rem = quot.quo_rem(delta)
        return rem

    # Compute psi using the recursive determinant formula
    def _psi_recursion(self, r, s):
        z = SR.var('z')
        M = matrix(SR, self.g + 1, self.g + 1, 0)
        # Fill out the first 2 columns
        for i in xrange(0, 2):
            for j in xrange(0, self.g + 1):
                M[j, i] = (self.compute_psi(s - self.g + j + i) * self.compute_psi(r - j + i))
        # Compute the remaining columns for each row
        for j in xrange(0, self.g+1):
            gp = self.compute_gamma(s - self.g + 1 + j) * self.compute_gamma(r + 1 - j)
            coeffs = gp.factor().coefficients(z)
            for i in xrange(0, self.g - 1):
                M[j, i+2] = _get_series_coeff(coeffs, i)

        D = M.determinant()

        if D == 0: raise ValueError("Determinant is zero in recursive formula")
        d = self.compute_psi(s - r)
        for k in xrange(2, self.g + 1):
            d *= self.compute_psi(r - self.g + k) * self.compute_psi(s - self.g + k)
        if d == 0: raise ValueError("Right hand side zero in recursion")
        return D / d

    def compute_psi(self, r):
        # Trivial cases
        if r < 0:
            raise ValueError("Psi_r only defined for non-negative r")
        elif r < self.g: return 0
        elif r == self.g: return 1

        # Base case
        if r < 5 * self.g:
            return self.compute_psi_from_determinant(r).factor()
        else:
            # Recursive case
            if (r - self.g) % 2 == 0:
                a = (r - self.g) / 2
                b = (r + self.g) / 2
            else:
                a = (r - self.g - 1) / 2
                b = (r + self.g - 1) / 2 + 1
            return self._psi_recursion(a,b).factor()

    def compute_P(self, r):
        psi = self.compute_psi(r)
        if (r - self.g) % 2 == 0:
            return psi
        else:
            y = SR.var('y')
            return (psi / (2*y)^self.g).substitute(y=self.G^(1/2))

    # Find n-torsion points by calculating Cantor's div polys
    def torsion_points(self, n, p):
        if n % 2 == 0:
            return itertools.chain(self._find_2_torsion_points(n,p),
                                   self._non_order_2_torsion_points(n,p))
        else:
            return self._non_order_2_torsion_points(n,p)

    def _non_order_2_torsion_points(self, n, p):
        gc = _reduce_poly(self.compute_P(n - self.g + 1), p)
        i = 1
        while n - self.g + 1 + i < n + self.g:
            gc = gc.gcd(_reduce_poly(self.compute_P(n - self.g + 1 + i), p))
            i += 1

        xvalues = []
        for irr in gc.factor():
            if irr[0].degree() == 1:
                x = irr[0].roots(multiplicities=False)[0]
                xvalues.append(x)
        return self._get_points_from_X_values(n, xvalues, p)

    def _get_points_from_X_values(self, n, xvalues, p):
        R = GF(p)['x']
        H = HyperellipticCurve(R(self.G))
        J = H.jacobian()

        torsion = set()
        for x in xvalues:
            nx = R(x)
            points = H.lift_x(nx, all=True)
            for point in points:
                 torsion.add(point)
        return list(torsion)

    # Need to handle 2 torsion points seperately
    def _find_2_torsion_points(self, n, p):
        R = GF(p)['x']
        H = HyperellipticCurve(R(self.G))

        lst = []

        for irr in R(self.G).factor():
            if irr[0].degree() == 1:
                x = irr[0].roots(multiplicities=False)[0]
                lst.append(H(x,0))
        return lst

    def _compute_C(self, r):
        m = floor((r + self.g) / 2)
        n = floor((r - self.g-1) / 2)

        coeffs = self._get_S_cached(m + n + self.g + 2)

        M = matrix(SR, n + 1, n + 1)
        C = 0

        if r <= self.g:
            C = 1
        else:
            z = SR.var('z')
            # The first n rows does not depend on j
            for xind in xrange(0, n + 1):
                for yind in xrange(0, n):
                    x1 = m-n+1 + xind + yind
                    s = _get_series_coeff(coeffs, x1)
                    M[yind,xind] = s
            for j in xrange(0, self.g + 1):
                # Set up last row
                for xind in xrange(0, n + 1):
                    s = _get_series_coeff(coeffs, m + 1 + j + xind)
                    M[n,xind] = s
                C += z^(m+n+j+1) * M.determinant()
            C = C / z^r
        return C

    def compute_gamma(self, r):
        z = SR.var('z')
        C = self._compute_C(r)
        # Normalize C
        vr = ((r+1)^2 - (r+1) -self.g^2 + self.g)/2
        gamma =  C.substitute(z=4*self.G*z) * 2^vr * self.G^(vr/2)
        return gamma

    # Find n-torsion points on the hyperelliptic curve y^2=f over F_p
    # Using "brute-forcing", i.e. looping over every element
    def torsion_points_brute_force(self, n, p):
        R = GF(p)['x']
        poly = R(self.G)
        H = HyperellipticCurve(poly)
        J = H.jacobian()

        torsion = set()

        for x in xrange(0, p):
            points = H.lift_x(x, all=True)
            for point in points:
                Jp = J(point)
                if n*Jp == 0:
                    torsion.add(point)
        return list(torsion)

    # Enumerate possible polynomials V corresponding to U such that (U,V) is a well-defined divisor in Mumford representation
    # using a brute-force approach
    def _naive_enum_mumford_V(self, U, p):
        R = GF(p)['x']
        x = var('x')
        for V in enumerate_polys(p, U.degree(x) - 1):
            q, r = R(V^2 - self.G).quo_rem(R(U))
            if r == 0:
                yield R(V)

    # More sophisticated approach in genus 2 case
    def _genus_2_enum_mumford_V(self, U, p):
        assert self.g == 2

        #print U

        def square_root(K, v):
            s = sqrt(v)
            if s in K:
                return s
            else:
                raise ValueError("Square root not in field!")

        # Returns (y0, y1), (-y0, -y1), (-y0, y1), (y0, -y1)
        def permute_signs(y0, y1, F):
            for i in [0, 1]:
                for j in [0, 1]:
                    yield (F((-1)^i*y0), F((-1)^j*y1))

        # The double root case needs another iteration to get a poly such that V_1^2 - f = 0 (mod (x-x0)^2)
        def double_root(U, f, x0, p):
            x = var('x')
            R = GF(p)
            K = R['x']
            try:
                y0 = square_root(R, R(f(x=x0)))
            except ValueError:
                return []
            if y0 == 0:
                return []
            W_1 = R(y0)
            P = K((W_1^2 - f)/(x-x0))
            k = -P(x=x0)/(2*y0)
            W_2 = W_1 + k*(x-x0)

            return list(set([ K(W_2), K(-W_2)]))

        # Do the computations in the splitting field of U
        def single_irreducible_factor(U, f, p):
            lst = []
            x, t = var('x'), var('t')
            R = GF(p)['t']
            K = GF(p^2, 't', modulus=U)
            F = K['x']
            a = R(f(x=t))
            t1 = K.gen()
            t2 = -F((U)/(x-t1)).list()[0]
            p1, p2 = K(a(t=t1)), K(a(t=t2))
            try:
                y0, y1 = square_root(K, p1), square_root(K, p2)
            except:
                return []

            assert K(a(t=t1)) == K(y0^2)
            assert K(a(t=t2)) == K(y1^2)
            assert K(R(U(x=t)).substitute(t=t1)) == 0
            assert K(R(U(x=t)).substitute(t=t2)) == 0

            x = F.gen()
            for (V_1, V_2) in permute_signs(y0, y1, F):
                V = crt(V_1, V_2, F(x-t1), F(x-t2))
                if V in GF(p)['x']:
                    lst.append(V)
            return list(set(lst))

        # In the case U splits into two linear factors, no field extension is neccessary
        def two_linear_factors(U, f, p, factors):
            lst = []
            K = GF(p)

            if factors[0][1] == 2: # Double root
                x = -factors[0][0].list()[0]
                return double_root(U, f, x, p)
            else:
                x0 = -factors[0][0].list()[0]
                x1 = -factors[1][0].list()[0]
            F = K['x']
            x = F.gen()
            try:
                y0 = square_root(K, K(f(x=x0)))
                y1 = square_root(K, K(f(x=x1)))
            except ValueError:
                return []
            for (V_1, V_2) in permute_signs(y0, y1, F):
                V = crt(F(V_1), F(V_2), F(x-x0), F(x-x1))
                if V in GF(p)['x']:
                    lst.append(V)
            return list(set(lst))

        def single_linear_factor(U, f, p, factors):
            K = GF(p)
            x0 = -factors[0][0].list()[0]
            F = K['x']
            x = F.gen()
            try:
                y0 = square_root(K, K(f(x=x0)))
            except ValueError:
                return []
            return list(set([ F(y0), F(-y0)]))

        lst = []
        R = GF(p)['x']
        factors = list(R(U).factor())
        n_factors = sum(map(lambda x: x[1], factors))
        # U is either (i) a constant,  (ii) singe linear factor, (iii) a single irreducible factor,
        # or (iv) factors into two linear factors
        if n_factors == 0:
            return [ R(0) ]
        elif n_factors == 1 and factors[0][0].degree() == 2:
                return single_irreducible_factor(U, R(self.G), p)
        elif n_factors == 1 and factors[0][0].degree() == 1:
                return single_linear_factor(U, R(self.G), p, factors)
        elif n_factors == 2:
            return two_linear_factors(U, R(self.G), p, factors)
        else:
            raise ValueError("Invalid U supplied")

    def _U_candidates(self, M, p):
        x = var('x')
        x0, x1, y0, y1 = var('x0', 'x1', 'y0', 'y1')

        s_num = SR(M.determinant().numerator())
        d = max(s_num.degree(y0), s_num.degree(y1))
        for k in xrange(1, d/2+1):
            s_num = s_num.substitute(y0^(2*k) == (SR(self.G).substitute(x=x0))^k,
                                     y1^(2*k) == (SR(self.G).substitute(x=x1))^k)
        # s_num is anti-symmetric, divide with the Vandemonde determinant
        assert bool(s_num.gcd(x0-x1) != 1)
        s_num = s_num / (x0-x1)

        # s_num is a symmetric polynomial in x0, x1.
        a = symmetric_poly_on_elementary_basis(s_num.factor()).factor()
        #import sympy
        #import sympy.polys.polyfuncs
        #a2 = sympy.polys.polyfuncs.symmetrize(s_num.factor(), formal=True)[0].factor()
        #print bool(SR(a1 == a2))
        #print "a1"
        #print a1
        #print "a2"
        #print a2
        
        R = GF(p)['s1', 's2']
        for i in GF(p):
            for j in GF(p):
                if R(a).substitute(s1=i, s2=j) == 0:
                    yield x^2 - i * x + j

    # Get r-divison points for hyperelliptic curve of genus 2
    # I.e. divisors D of J(C)/F_p such that r*D~0.
    def division_points(self, r, p):
        # Only works for genus 2 curves
        assert self.g == 2
        # Algorithm does not work for even r since the formula (8.35) does not hold for order 2 points.
        assert r % 2 == 1

        R2 = GF(p)['x']
        H = HyperellipticCurve(R2(self.G))
        J = H.jacobian()

        # (1) = (\infty) is always a division point
        i0 = [ J(0) ]
        # Divisor points of the form D=(x,y)
        i1 = itertools.imap(lambda x: J(x), self.torsion_points(r, p))
        # Divisor points of the form D=2(x,y)
        i2 = itertools.imap(lambda x: 2*J(x), self.torsion_points(2*r, p))
        # Drop the second redundant (1)
        i2 = itertools.islice(i2, 1, None)
        # Divisor points of the form D=(x0,y0)+(x1,y1), x0!=x1
        i3 = self._two_term_division_points(r, p)

        return itertools.chain(i0, i1, i2, i3)

    # Find divisor points on the form D=(x0,y0)+(x1,y1), where y0!=0 and y1!=0
    def _two_term_division_points(self, r, p):
        assert self.g == 2

        x, y, z, X = var('x','y', 'z', 'X')
        x0, x1, y0, y1 = var('x0', 'x1', 'y0', 'y1')
        e1, e2 = var('e1', 'e2')

        R = QQ
        K = R['x0,y0,x1,y1']
        F = Frac(K)
        G = PolynomialRing(F, 'X')

        d1 = G(self._compute_delta(r).substitute(z=(x-X)/(4*y^2)).substitute(x=x0, y=y0).factor())
        d2 = G(self._compute_delta(r).substitute(z=(x-X)/(4*y^2)).substitute(x=x1, y=y1).factor())
        e1 = G(self._compute_epsilon(r).substitute(z=(x-X)/(4*y^2)).substitute(x=x0, y=y0).factor())
        e2 = G(self._compute_epsilon(r).substitute(z=(x-X)/(4*y^2)).substitute(x=x1, y=y1).factor())

        d1m, e1m = _get_monic_cantor_pol(d1, e1)
        d2m, e2m = _get_monic_cantor_pol(d2, e2)

        M3 = matrix(F, 2, 2)
        M3[0, 0] = d1m.coefficients(X)[0]
        M3[0, 1] = d2m.coefficients(X)[0]
        M3[1, 0] = d1m.coefficients(X)[1]
        M3[1, 1] = d2m.coefficients(X)[1]

        M4 = matrix(F, 2, 2)
        M4[0, 0] = e1m.coefficients(X)[0]
        M4[0, 1] = e2m.coefficients(X)[0]
        M4[1, 0] = e1m.coefficients(X)[1]
        M4[1, 1] = e2m.coefficients(X)[1]

        # (delta_0, epsilon_0) and (delta_1, epsilon_1) are propotional if r*D=r*(delta_0, epsilon_0) + r*(delta_1, epsilon_1)=0
        lst = []
        R2 = GF(p)['x']
        H = HyperellipticCurve(R2(self.G))
        J = H.jacobian()
        set1 = set(self._U_candidates(M3, p))
        #print set1
        set2 = set(self._U_candidates(M4, p))
        #print set2

        # Vanishing of both determinants is a neccessary condition, so the intersection gives us a candidate set
        U_candidates = list(set1.intersection(set2))
        #print U_candidates
        #print "!!len: " + str(len(U_candidates))
        for U in U_candidates:
            for V in self._enum_mumford_V(U, p):
                jp = J(R2(U), R2(V))
                a = r * jp
                #print "(U,V)=" + str((U,V)) + " r * J(U,V)=" + str(a)
                if a == 0:
                    yield jp

    # Find division points using a naive, brute-force algorithm
    def division_points_naive(self, r, p):
        R = GF(p)['x']
        x= R.gen()
        H = HyperellipticCurve(R(self.G))
        J = H.jacobian()

        lst = []
        for U in enumerate_monic_polys(p, self.g):
            for V in self._enum_mumford_V(U, p):
                d = J(R(U), R(V))
                if r*d == 0:
                    if not d in lst:
                        lst.append(d)
        return lst

    # unit tests

    # Check we get the data in the example in Section 10
    def _test_cantor_example(self, tester):
        x = SR.var('x')
        z = SR.var('z')
        tester.assert_(self.g == 2)
        tester.assert_(self.G == x^5 + 1)
        tester.assert_(self.compute_psi(1) == 0)
        tester.assert_(self.compute_psi(2) == 1)
        tester.assert_(self.compute_psi(3) == 4*(x^4 - x^3 + x^2 - x + 1)*(x + 1))
        tester.assert_(self.compute_psi(4) == 10*(x^5 - 4)^2*x^2)
        tester.assert_(self.compute_psi(5) == 20*(x^10 - 108*x^5 + 16)*(x^5 - 4)*(x^4 - x^3 + x^2 - x + 1)*(x + 1)*x)
        tester.assert_(self.compute_psi(6) == 5*(7*x^20 - 4872*x^15 - 7408*x^10 - 12672*x^5 - 768)*(x^5 - 4)^2*x^2)
        tester.assert_(self.compute_psi(7) == 8*(7*x^40 - 22344*x^35 + 224896*x^30 - 9451008*x^25 + 170240*x^20 - 57028608*x^15 - 3272704*x^10 - 4767744*x^5 - 32768)*(x^4 - x^3 + x^2 - x + 1)*(x + 1))
        tester.assert_(self.compute_psi(8) == 84*x^60 - 981792*x^55 - 8072256*x^50 - 5136395520*x^45 - 9367057920*x^40 - 90223337472*x^35 - 114657705984*x^30 - 117268611072*x^25 - 155547729920*x^20 - 229528043520*x^15 - 38803603456*x^10 + 1686110208*x^5 - 16777216)

    # Test the recursion formula. The first recursion case occurs at r=5*g.
    def _test_recursion(self, tester):
        tester.info("testing to recursively compute psi_" + str(5*self.g) + "\n")
        tester.assert_(self.compute_psi(5*self.g) == self.compute_psi_from_determinant(5*self.g).factor())
#        tester.info("testing to recursively compute psi_" + str(5*self.g+1) + "\n")
#        tester.assert_(self.compute_psi(5*self.g + 1) == self.compute_psi_from_determinant(5*self.g + 1).factor())

    # Find the torsion points using Cantor's div polys and compare with the naive brute-force algorithm using Cantor's algorithm
    def _test_compare_torsion_algos(self, tester):
        p = 3
        while p < 20:
            p = next_prime(p)
            for n in xrange(2, 10):
                try:
                    torsion1 = self.torsion_points_brute_force(n, p)
                    torsion2 = self.torsion_points(n, p)
#                    tester.info(str(p) + "," + str(n) + "\n" + str(torsion1) + "\n" + str(torsion2) + '\n')
                    tester.assert_(sorted(torsion1) == sorted(torsion2))
                except ValueError:
#                    tester.info('skipping singular curve')
                    continue

    # Test constant and leading coefficients of delta. According to (8.11), (8.12) the following assertions should hold
    def _test_delta_constant_and_leading_coeff(self, tester):
        n = 3
        z = SR.var('z')
        constant_coeff = _get_series_coeff(self._compute_delta(n).coefficients(z), 0)
        exp1 = -(self.compute_psi(n-1)*self.compute_psi(n+1))
#        print constant_coeff; print exp1;
        tester.assert_((constant_coeff == 0 and exp1 == 0) or (constant_coeff.factor() == exp1.factor()))

        leading_coeff = _get_series_coeff(self._compute_delta(n).coefficients(z), self.g).factor()
        exp2 = (-(4*self.G)^(self.g) * self.compute_psi(n)^2).factor()
#       print leading_coeff; print exp2;
        tester.assert_((leading_coeff == 0 and exp2 == 0) or leading_coeff.factor() == exp2.factor())

    # In the case g=1 we should get back our classical elliptic curve divsion polynomials
    def _test_genus1(self, tester):
        x = SR.var('x')
        y = SR.var('y')
        z = SR.var('z')

        H = CanHyperCurve(1, x^3+1)
        assert H.compute_psi(0) == 0
        assert H.compute_psi(1) == 1
        assert H.compute_psi(2) == 2 * sqrt(x^3+1)
        assert H.compute_psi(3) == 3*(x^3+4)*x
        assert H.compute_psi(4) == 4*(x^4 - 2*x^3 + 6*x^2 + 4*x + 4)*sqrt(x^3 + 1)*(x^2 + 2*x - 2)
        assert H._compute_delta(3).factor() == (H.compute_psi(2)*H.compute_psi(4) - 4 * H.G * z * H.compute_psi(3)^2).factor()
        assert SR(H._compute_epsilon(3)).factor() == ((H.compute_psi(5)*y - 4 * (x^4 - 2*x^3 + 6*x^2 + 4*x + 4)^2* y * (x^2 + 2*x - 2)^2)/(H.compute_psi(3)^3)).factor()

    def _test_epsilon_delta_formula(self, tester):
        for n in [3, 5, 7]:
            R = QQ['x']
            x = R.gen()
            H1= HyperellipticCurve(x^5+17)
            J = H1.jacobian()
            U,V = n*J(H1(2,7))
            #print U; print V

            X = SR.var('X')
            x = SR.var('x')
            y = SR.var('y')
            z = SR.var('z')

            H = CanHyperCurve(2, x^5 + 17)
            delta = R(H._compute_delta(n).substitute(z=(x-X)/(4*y^2)).substitute(x=2, y=7).factor().substitute(X=x))
            epsilon = R(H._compute_epsilon(n).substitute(z=(x-X)/(4*y^2)).substitute(x=2, y=7).factor().substitute(X=x))

            delta, epsilon = _get_monic_cantor_pol(delta, epsilon)
            #print delta; print epsilon
            tester.assert_(U == delta)
            tester.assert_(V == epsilon)

    def _test_division_points(self, tester):
        x = SR.var('x')
        H = CanHyperCurve(2, x^5+1)
        tester.assert_(str(list(H.division_points(5,7))) == "[(1), (x, y + 6), (x, y + 1), (x^2, y + 6), (x^2, y + 1)]")

    def _test_cmp_division_points_algos(self, tester):
        tester.assert_(str(sorted(list(self.division_points_naive(3, 7)))) == str(sorted(list(self.division_points(3,7)))))
        tester.assert_(str(sorted(list(self.division_points_naive(3, 11)))) == str(sorted(list(self.division_points(3,11)))))
        tester.assert_(str(sorted(list(self.division_points_naive(5, 11)))) == str(sorted(list(self.division_points(5,11)))))

    def _test_enum_mumford_V1(self, tester):
        H=CanHyperCurve(2, x^5 + 4*x+5)
        assert sorted(list(H._naive_enum_mumford_V(x^2+x, 7))) == sorted(H._genus_2_enum_mumford_V(x^2+x,7))
        assert sorted(list(H._naive_enum_mumford_V(x^2, 7))) == sorted(H._genus_2_enum_mumford_V(x^2,7))

        H=CanHyperCurve(2, x^5 + 7*x^3 + 4*x^2+1)
        assert sorted(list(H._naive_enum_mumford_V(x^2+1,19))) == sorted(H._genus_2_enum_mumford_V(x^2+1, 19))
        assert sorted(list(H._naive_enum_mumford_V(x^2,19))) == sorted(H._genus_2_enum_mumford_V(x^2, 19))
        assert sorted(list(H._naive_enum_mumford_V(x^2-4*x+4,19))) == sorted(H._genus_2_enum_mumford_V(x^2-4*x+4, 19))

        H=CanHyperCurve(2, x^5 + 11*x^3 + 5*x^2+7*x+1)
        assert sorted(list(H._naive_enum_mumford_V(x^2+1,19))) == sorted(H._genus_2_enum_mumford_V(x^2+1, 19))
        assert sorted(list(H._naive_enum_mumford_V(x^2+3,19))) == sorted(H._genus_2_enum_mumford_V(x^2+3, 19))
        assert sorted(list(H._naive_enum_mumford_V(x^2+2*x,19))) == sorted(H._genus_2_enum_mumford_V(x^2+2*x, 19))
        assert sorted(list(H._naive_enum_mumford_V(x^2,19))) == sorted(H._genus_2_enum_mumford_V(x^2, 19))
        assert sorted(list(H._naive_enum_mumford_V(x^2+4*x+4,19))) == sorted(H._genus_2_enum_mumford_V(x^2+4*x+4, 19))

    def _test_enum_mumford_V2(self, tester):
        H=CanHyperCurve(2, x^5 + 7*x^3 + 4*x^2+1)
        for U in enumerate_monic_polys(5, 2):
            assert sorted(list(H._naive_enum_mumford_V(U,5))) == sorted(H._genus_2_enum_mumford_V(U, 5))
        for U in enumerate_monic_polys(7, 2):
            assert sorted(list(H._naive_enum_mumford_V(U,7))) == sorted(H._genus_2_enum_mumford_V(U, 7))
        H=CanHyperCurve(2, x^5 + x + 5)
        for U in enumerate_monic_polys(11, 2):
            assert sorted(list(H._naive_enum_mumford_V(U,11))) == sorted(H._genus_2_enum_mumford_V(U, 11))
