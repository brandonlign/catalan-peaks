"""Exact check of the theorem for small n (Section 9 of the paper).

For each n from 1 to N (default 420) let f = C_{n+1} (newer row) and g = C_n
(older row), as integer polynomials in X = t - 1 by the closed formula. Run the
Euclidean algorithm on f, g with negated remainders, rescaling each remainder by
a positive rational so that it has coprime integer coefficients. By Lemma 9.1 the rows
interlace strictly, with the rightmost zero in the newer row, if every
polynomial in this chain has positive leading coefficient and each has degree
one less than the one before it (after first reducing g modulo f when both have
the same degree). Simple negative zeros follow because the coefficients in t are
nonnegative and C_n(0) > 0. Expected output: progress every 50 and a final line. Plain Python integers
work; with gmpy2 installed the full run is much faster.
"""
import sys
from rows import coeffs_X

try:  # gmpy2 (GMP integers) makes the run about an order of magnitude faster
    from gmpy2 import mpz, gcd
except ImportError:
    from math import gcd
    mpz = int


def primitive(p):
    g = 0
    for c in p:
        g = gcd(g, c)
        if g == 1:
            return p
    return [c // g for c in p] if g > 1 else p


def remainder(a, b):
    """A positive multiple of the remainder of a modulo b (b has positive leading coefficient)."""
    a = list(a)
    lb = b[-1]
    while len(a) >= len(b) and any(a):
        la, k = a[-1], len(a) - len(b)
        a = [c * lb for c in a]
        for i, c in enumerate(b):
            a[k + i] -= la * c
        a.pop()
        while len(a) > 1 and a[-1] == 0:
            a.pop()
    return primitive(a)


def interlaces(f, g):
    if len(f) == len(g):
        g = remainder(g, f)
    while len(f) > 1:
        if len(g) != len(f) - 1 or g[-1] <= 0:
            return False
        f, g = g, [-c for c in remainder(f, g)]
    return True


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 420
    row = lambda n: primitive([mpz(c) for c in coeffs_X(n)])
    old = row(1)
    for n in range(1, N + 1):
        new = row(n + 1)
        assert interlaces(new, old), n
        if n % 50 == 0:
            print(f"pairs 1..{n}: ok", flush=True)
        old = new
    print(f"C_n and C_(n+1) interlace strictly for all 1 <= n <= {N}")
