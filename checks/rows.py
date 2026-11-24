"""The short-peak polynomials C_n(t) of Catalan words, computed three ways.

A Catalan word of length n starts at 0 and has w_{i+1} <= w_i + 1, all w_i >= 0.
A short peak is an index i with w_{i-1} < w_i > w_{i+1}. C_n(t) counts words by
short peaks. coeffs_X uses the closed formula of Lemma 2.2 of the paper; brute and
series_rows give the same rows by direct enumeration and from the generating
function of Lemma 2.1.
"""
from fractions import Fraction
from math import comb, factorial


def coeffs_X(n):
    """Coefficients of C_n(1 + X) in powers of X (Lemma 2.2)."""
    if n == 0:
        return [1]
    return [factorial(2 * n - 3 * i) // (factorial(i) * factorial(n - 2 * i - 1) * factorial(n - 2 * i)
                                          * (n - i) * (n - i + 1)) for i in range((n - 1) // 2 + 1)]


def shift(p, s):
    """Coefficients of p(x + s) from those of p(x), exactly."""
    out = [0] * len(p)
    for i, c in enumerate(p):
        for k in range(i + 1):
            out[k] += c * comb(i, k) * s ** (i - k)
    return out


def row(n):
    """Coefficients of C_n(t) in powers of t."""
    return shift(coeffs_X(n), -1)


def brute(n):
    counts = {}

    def go(w):
        if len(w) == n:
            k = sum(1 for i in range(1, n - 1) if w[i - 1] < w[i] > w[i + 1])
            counts[k] = counts.get(k, 0) + 1
            return
        for v in range(w[-1] + 2):
            go(w + [v])

    go([0])
    return [counts.get(k, 0) for k in range(max(counts) + 1)]


def series_rows(N):
    """C_n(t) for n <= N from G = (1 - a x^2 - sqrt((1 + a x^2)^2 - 4x)) / (2x(1 - ax)), a = 1 - t.
    Power series in x with coefficients polynomials in t (lists), by exact arithmetic."""
    from polyarith import padd, pmul, pscale
    a = [1, -1]
    # D = (1 + a x^2)^2 - 4x as series in x.
    D = [[0] for _ in range(N + 3)]
    D[0] = [1]
    D[1] = [-4]
    D[2] = pscale(a, 2)
    D[4] = pmul(a, a)
    # sqrt(D) = S with S_0 = 1: 2 S_0 S_k = D_k - sum_{0<i<k} S_i S_{k-i}.
    S = [[1]]
    for k in range(1, N + 3):
        acc = D[k] if k < len(D) else [0]
        for i in range(1, k):
            acc = padd(acc, pscale(pmul(S[i], S[k - i]), -1))
        S.append([Fraction(c, 2) for c in acc])
    num = [pscale(s, -1) for s in S]
    num[0] = padd(num[0], [1])
    num[2] = padd(num[2], pscale(a, -1))
    # Divide by 2x, then by (1 - a x).
    Q = [pscale(num[k + 1], Fraction(1, 2)) for k in range(N + 1)]
    out = []
    for k in range(N + 1):
        cur = Q[k]
        if k:
            cur = padd(cur, pmul(a, out[-1]))
        out.append(cur)
    return out
