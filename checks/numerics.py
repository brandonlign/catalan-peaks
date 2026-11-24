"""Optional floating-point illustration of Sections 4-7 (needs mpmath).

1. The normalized values a0^(-j) p_j(a0 + v) of Lemma 4.1, against 11/16 - kappa v/h^j.
2. Proposition 5.1: p_n(a) against beta_m/(pi b^m) Re(q^(-m) Z_q S_m), for
   several n and angles; the ratio should be 1.
3. |S_m - 1|, the phase gap (W_(m+1) - W_m)/alpha, W_m, and the difference between
   arg Z_q and the formula of Lemma 7.1, for m = 421, 422 and 1000.

None of this is used in the proof; constants.py checks the inequalities exactly.
Runs in a few minutes.
"""
from fractions import Fraction as F
import mpmath as mp
from rows import coeffs_X

mp.mp.dps = 40


def p_exact(n, a):
    """p_n(a) = C_n(1 - a) = sum_i c_i (-a)^i, with c_i the coefficients of C_n(1 + X)."""
    return sum(c * (-a) ** i for i, c in enumerate(coeffs_X(n)))


def to_mp(x):
    return mp.mpf(x.numerator) / x.denominator


def params(al):
    c = mp.cos(al)
    R = 1 + 2 * c
    b = 2 * (1 + c) * c ** 2 / R ** 2
    a = 1 / (R * b ** 2)
    return c, R, b, a, a * b, mp.expj(al), -c + 1j * mp.sqrt(R ** 2 - c ** 2)


def Z_and_L(d, roots, u):
    Z = 1
    for h in roots:
        if h != d:
            Z *= mp.sqrt(1 - d / h)
    Z /= 1 - u * d

    def L(y):
        v = (1 - u * d) / (1 - u * d * (1 + y))
        for h in roots:
            if h != d:
                v *= mp.sqrt(1 - y * d / (h - d))
        return v
    return Z, L


def S_of(m, al):
    c, R, b, a, u, q, eta = params(al)
    roots = [q, mp.conj(q), eta, mp.conj(eta)]
    Zq, Lq = Z_and_L(q, roots, u)
    Ze, Le = Z_and_L(eta, roots, u)
    beta = mp.beta(mp.mpf(3) / 2, m - mp.mpf(1) / 2)
    E = lambda L: mp.quad(lambda y: mp.sqrt(y) * (1 + y) ** (-m - 1) * L(y),
                          [0, mp.mpf(1) / m, mp.mpf(10) / m, 1, mp.inf]) / beta
    S = E(Lq) + (q / eta) ** m * (Ze / Zq) * E(Le)
    return Zq, S, beta, b, q, a


def phi(al):
    c, R, b, a, u, q, eta = params(al)
    H = 1 + (2 * c * q + q ** 2) / R ** 2
    return al / 2 - mp.pi / 4 + mp.arg(H) / 2 - mp.arg(1 - u * q)


if __name__ == "__main__":
    a0, h = F(27, 16), F(3, 4)
    mp.mp.dps = 300
    kappa = 16 * mp.sqrt(3) / 27
    worst = 0
    for j in (420, 421, 500):
        for k in range(5):
            v = h ** j * F(k, 4)
            val = to_mp(p_exact(j, a0 + v)) / to_mp(a0) ** j
            worst = max(worst, abs(val - mp.mpf(11) / 16 + kappa * k / 4))
    print("1. Lemma 4.1: largest deviation for j in {420, 421, 500}:", mp.nstr(worst, 3))
    mp.mp.dps = 40

    print("2. Proposition 5.1: p_n(a) / formula")
    for n in (5, 40, 420):
        m = n + 1
        out = []
        for al in (mp.mpf("0.05"), mp.mpf("0.6"), mp.mpf("1.3")):
            Zq, S, beta, b, q, a = S_of(m, al)
            lhs = to_mp(p_exact(n, F(mp.nstr(a, 40))))
            rhs = beta / (mp.pi * b ** m) * mp.re(q ** (-m) * Zq * S)
            out.append(mp.nstr(lhs / rhs, 8))
        print(f"   n={n}: ratios {out}")

    print("3. Sections 6-7")
    for m in (421, 422, 1000):
        print(f"   m={m}")
        for al in (mp.mpf(3) / m, mp.mpf(5) / m, mp.mpf("0.1"), mp.mpf("0.5"), mp.mpf("1.5"),
                   mp.acos(mp.mpf(2) / m)):
            Zq, S1, *_ = S_of(m, al)
            _, S2, *_ = S_of(m + 1, al)
            f = phi(al)
            W1 = m * al - f - mp.arg(S1)
            W2 = (m + 1) * al - f - mp.arg(S2)
            print("     alpha=%.4f  |S_m-1|=%.4f  gap/alpha=%.4f  W_m=%.3f  argZ-phi=%.1e"
                  % (al, abs(S1 - 1), (W2 - W1) / al, W1, mp.arg(Zq) - f))
