"""Exact checks of every numerical inequality in Sections 4-7 of the paper.

All arithmetic is with Python Fractions. Square roots, pi, sin, cos and exp
enter only through rational lower and upper bounds (truncated Taylor series
with their alternating-series remainders, integer square roots), so each
line below is a proof of the stated inequality. Expected output: one line per
group of checks, then "all constants ok".
"""
from fractions import Fraction as F
from math import isqrt

SCALE = 10 ** 40
PI_LO = F(314159265358979323846, 10 ** 20)
PI_HI = F(314159265358979323847, 10 ** 20)

count = 0


def ok(cond, label):
    global count
    if not cond:
        raise AssertionError(label)
    count += 1


# ---------- rational bounds for irrational quantities ----------

def sqrt_lo(x):
    x = F(x)
    return F(isqrt(x.numerator * SCALE ** 2 // x.denominator), SCALE)


def sqrt_hi(x):
    return sqrt_lo(x) + F(1, SCALE)


def _series(x, start, terms):
    # sum_{k<terms} (-1)^k x^(start+2k) / (start+2k)!, and the next term's size
    s, t, k = F(0), F(1), 0
    for j in range(1, start + 1):
        t *= F(x) / j
    for k in range(terms):
        s += t if k % 2 == 0 else -t
        n = start + 2 * k
        t = t * F(x) * F(x) / ((n + 1) * (n + 2))
    return s, t


def sin_bounds(x, terms=12):
    # valid for 0 <= x <= 2: the terms decrease from the second one on
    s, nxt = _series(x, 1, terms)
    return s - nxt, s + nxt


def cos_bounds(x, terms=12):
    s, nxt = _series(x, 0, terms)
    return s - nxt, s + nxt


def exp_lo(x, terms=60):
    s, t = F(0), F(1)
    for k in range(terms):
        s += t
        t = t * F(x) / (k + 1)
    return s  # x >= 0: partial sums are lower bounds


SQ3_LO, SQ3_HI = sqrt_lo(3), sqrt_hi(3)
KAPPA_LO, KAPPA_HI = 16 * SQ3_LO / 27, 16 * SQ3_HI / 27


def moment(m, j):
    """E_m y^j = prod_{i=1}^j (2i+1)/(2m-1-2i), the moment formula of Section 5."""
    r = F(1)
    for i in range(1, j + 1):
        r *= F(2 * i + 1, 2 * m - 1 - 2 * i)
    return r


def poly_mul(p, q):
    out = [F(0)] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            out[i + j] += x * y
    return out


def poly(*factors):
    p = [F(1)]
    for f in factors:
        p = poly_mul(p, [F(c) for c in f])
    return p


def expect(m, p):
    return sum(c * moment(m, j) for j, c in enumerate(p))


# ---------- Section 4: the first two zeros ----------

def micro():
    a0, h, rho = F(27, 16), F(3, 4), F(4, 9)
    ok(rho * a0 == h, "rho a0 = h")
    x0 = 1 / a0
    ok(a0 * x0 ** 2 == F(16, 27), "a0 x0^2")
    ok(1 - F(9, 4) * x0 == F(-1, 3) and 1 + x0 / 2 + F(9, 16) * x0 ** 2 == F(121, 81), "D(a0,x0) factors")
    num = 1 - a0 * x0 ** 2 - F(-1, 3) * F(11, 9)
    ok(num == F(22, 27) and num / (2 * x0) == F(11, 16), "residue 11/16")
    # the factorization of D(a0, x), checked at sample points (both sides are quartics; 5 points suffice)
    for x in (F(0), F(1, 3), F(2, 7), F(-5, 11), F(7, 2)):
        ok((1 + a0 * x * x) ** 2 - 4 * x == (1 - F(9, 4) * x) ** 2 * (1 + x / 2 + F(9, 16) * x * x), "factorization of D(a0, x)")
    t = lambda r: F(r) / 2 + F(9, 16) * F(r) ** 2
    # radius 0.8
    r = F(8, 10)
    ok(r < F(4, 3), "0.8 < 4/3")
    ok(r * a0 - 1 >= F(35, 100), "|1-a0x| >= 0.35")
    sD = (1 + r / rho) * sqrt_hi(1 + t(r))
    ok(sD < F(38, 10), "|sqrt D| < 3.8 on 0.8")
    ok((1 + r * r * a0 + F(38, 10)) / (2 * r * F(35, 100)) < 13, "|G| < 13 on 0.8")
    ok(13 + F(11, 16) / F(35, 100) < 16, "|E1| < 16")
    # radius 0.52: the pole of 1/sqrt(D) at rho
    ok(a0 * rho ** 2 == F(1, 3) and 1 + rho / 2 + F(9, 16) * rho ** 2 == F(4, 3) and 1 - a0 * rho == F(1, 4), "values at rho")
    # principal part: -rho(1+a0 rho^2)/(2(1-a0 rho) * (2/sqrt3)) = -16 sqrt3/27
    ok(rho * (1 + a0 * rho ** 2) / (2 * (1 - a0 * rho)) == F(32, 27), "kappa = (32/27)(sqrt3/2)")
    r = F(52, 100)
    ok(r < 1 / a0, "0.52 < 1/a0")
    ok(1 - r * a0 > F(12, 100), "|1-a0x| > 0.12")
    lo = (r / rho - 1) * sqrt_lo(1 - t(r))
    hi = (1 + r / rho) * sqrt_hi(1 + t(r))
    ok(lo > F(1, 10) and hi < F(31, 10), "0.1 < |sqrt D| < 3.1 on 0.52")
    G = (1 + r * r * a0 + F(31, 10)) / (2 * r * F(12, 100))
    ok(G < 38, "|G| < 38 on 0.52")
    Ga = (r * r + r * r * (1 + a0 * r * r) / F(1, 10)) / (2 * r * F(12, 100)) + r * 38 / F(12, 100)
    ok(Ga < 250, "|G_a| < 250 on 0.52")
    ok(250 + KAPPA_HI / (r / rho - 1) < 300, "|E2| < 300")
    # radius 0.4 and |a - a0| <= 1e-3
    r, e = F(4, 10), F(1, 1000)
    D0 = (1 - r / rho) ** 2 * (1 - t(r))
    ok(D0 == F(71, 10000), "|D(a0,x)| >= 0.0071")
    dD = e * r * r * (2 + (2 * a0 + e) * r * r)
    ok(dD < F(42, 100000) < D0, "|D - D0| < 0.00042")
    amax = a0 + e
    ok((1 + r * r * amax) ** 2 + 4 * r < F(324, 100), "|D| < 3.24")
    ok(1 - r * amax >= F(32, 100), "|1-ax| >= 0.32")
    ok((1 + F(17, 100) * amax + F(18, 10)) / (2 * r * F(32, 100)) < 20, "|G| < 20 on the polydisk")
    # the three error terms at j = 420
    b1, b2, b3 = 1 / (F(8, 10) * a0), rho / F(52, 100), h * h / (F(4, 10) * a0)
    ok(b1 < 1 and b2 < 1 and b3 < 1, "bases < 1")
    j = 420
    ok(1000 * h ** j < F(1, 2), "1e3 h^j < 1/2")
    ok(16 * b1 ** j + 300 * b2 ** j + 4 * 10 ** 7 * b3 ** j < F(1, 1000), "micro error < 1e-3")
    # Corollary 4.2
    ok(F(11, 16) - F(58, 100) * KAPPA_HI > F(1, 1000), "11/16 - 0.58 kappa > 1e-3")
    ok(F(58, 75) <= 1 and F(11, 16) - F(58, 75) * KAPPA_LO < -F(1, 1000), "11/16 - (0.58/0.75) kappa < -1e-3")
    n = 420
    ok(F(58, 100) * h ** n * (n + 1) ** 2 < 3 * a0, "s (n+1)^2 < 3 a0")
    ok(h * F(n + 2, n + 1) ** 2 < 1, "s (n+1)^2 decreasing")
    ok(1 - F(1, 2) + F(1, 24) - F(1, 720) >= F(1, 3), "1 - cos x >= x^2/3 on [0,1]")
    print(f"Section 4 (first two zeros): ok")


# ---------- Section 6: bounds for S ----------

def middle():
    m = 421
    cl, ch = cos_bounds(F(3, 2))
    sl, _ = sin_bounds(F(1, 4))
    ok(1 + 2 * cl > F(114, 100), "R >= 1.14")
    ok(sl > F(247, 1000), "sin(1/4) > 0.247")
    ok(1 / (2 * F(247, 1000)) < 3 and 1 / F(14, 100) < 8 and 1 / F(247, 1000) < 5, "constants for q")
    ok(expect(m, poly([1, 3], [1, 8], [1, 8], [1, 5])) - 1 < F(93, 1000), "E(1+3y)(1+8y)^2(1+5y) - 1 < 0.093")
    # c/R = c/(1+2c) <= 1/3, so sin(theta) >= sqrt(8)/3
    st = sqrt_lo(8) / 3
    ok(1 / (2 * st) < 1 and F(114, 100) / F(14, 100) < 10 and 1 / st < 2, "constants for eta")
    ok(expect(m, poly([1, 10], [1, 10], [1, 10], [1, 10])) < 2, "E(1+10y)^4 < 2")
    ok(6 / sqrt_lo(F(247, 1000)) < 13, "amplitude < 13")
    ok(13 * 2 * (1 / F(114, 100)) ** m < F(1, 10 ** 20), "second term < 1e-20")
    print("Proposition 6.4 (middle range): ok")


def large():
    m = 421
    cl, ch = cos_bounds(F(3, 2))
    sl, _ = sin_bounds(F(3, 2))
    ok(ch < F(71, 1000) and 1 + 2 * ch < F(115, 100), "c < 0.071, R < 1.15")
    ok(F(2, m) < cl, "alpha_* > 3/2")
    ok(sl > F(997, 1000), "sin(3/2) > 0.997")
    ok(1 / F(115, 100) > F(86, 100), "sin(theta) > 0.86")
    ok(all(x <= 2 for x in (1 / (2 * F(997, 1000)), 1 / F(199, 100), 1 / F(997, 1000),
                            1 / (2 * F(86, 100)), F(115, 100) / F(199, 100), 1 / F(86, 100))), "constants <= 2")
    ok(F(115, 100) / 4 <= F(288, 1000), "R/(2c) <= 0.288 m")
    ok(expect(m, poly([1, 2], [1, 2], [1, 2], [1, F(288, 1000) * m])) - 1 < F(48, 100), "E(1+2y)^3(1+0.288my) - 1 < 0.48")
    ok(2 * F(115, 100) / sqrt_lo(F(997, 1000)) < F(231, 100), "amplitude < 2.31")
    ok(F(4 * m, m + 2) > F(398, 100), "j log R > 3.98")
    ok(exp_lo(F(398, 100)) > 1 / F(19, 1000), "exp(-3.98) < 0.019")
    ok(F(48, 100) + F(19, 1000) * F(231, 100) * F(148, 100) < F(57, 100), "|S-1| < 0.57")
    print("Proposition 6.5 (large angles): ok")


def small():
    m = 421
    cl, _ = cos_bounds(F(1, 4))
    ok(cl >= F(31, 32), "cos(1/4) >= 31/32")
    ok(1 - F(1, 16) / 6 > F(25, 26), "sin(alpha) > 25 alpha/26")
    c = F(31, 32)
    ok(1 + 2 * c == F(47, 16), "R >= 47/16")
    ok((1 + 2 * c) / (2 * (1 + c) * c * c) <= F(4, 5), "u <= 4/5")
    ok(F(26, 75) / 4 == F(13, 150), "1/(4 sin alpha) <= 13m/150")
    ok(1 / (F(47, 16) - 1) < 1, "1/(R-1) < 1")
    ok(1 - F(4, 5) * F(9, 8) == F(1, 10) and F(4, 5) / F(1, 10) == 8, "pole constant 8 for y <= 1/8")
    ok(F(26, 75) < 2, "1/sin(alpha) < 2m")
    ok(F(2, 31) < 1, "int y^(1/2)(1+y)^(-17) <= int (1+y)^(-33/2) = 2/31 < 1")
    # (1 + 1/j)^(j+1) <= 4 for j >= 1 gives beta_j >= (2/3) j^(-3/2) / 4 = 1/(6 j^(3/2))
    ok(F(422, 421) ** 422 <= 4, "(1+1/j)^(j+1) <= 4")
    def tau(m):
        # 6 (m+1)^(3/2) (1+2m)^2 (8/9)^(m-20); bound (m+1)^(3/2) by (m+1)^2
        return 6 * (m + 1) ** 2 * (1 + 2 * m) ** 2 * F(8, 9) ** (m - 20)
    ok(tau(m) < F(1, 10 ** 7), "tau_421 < 1e-7")
    P = poly([1, 8], [1, 1], [1, 1], [1, F(13 * m, 150)])
    P[0] -= 1
    EP = expect(m, P)
    ok(EP < F(175, 1000), "E P_m < 0.175")
    DP = expect(m, poly_mul([F(3, 2 * m - 1), F(1)], P))
    ok(F(m, 3) * DP < F(237, 1000), "(m/3) E[(3/(2m-1)+y) P_m] < 0.237")
    ok(F(16, 31) * F(47, 16) == F(47, 31) and F(47, 31) <= 2 and 3 / sqrt_lo(8) <= 2, "eta constants <= 2")
    ok(expect(m, poly([1, 2], [1, 2], [1, 2], [1, 2])) < 2, "E(1+2y)^4 < 2")
    ok(36 * F(26, 75) * m < 16 * m * m, "6 sqrt(26m/75) < 4m")
    sec = 8 * m * F(16, 47) ** m
    ok(F(175, 1000) + F(1, 10 ** 7) + sec < F(176, 1000), "|S-1| < 0.176")
    ok(F(m, 3) * (tau(m) + 2 * sec) < F(1, 1000), "(m/3)(tau + 16m(16/47)^m) < 1e-3")
    # the last bound decreases in m: each factor's ratio m -> m+1 is below 1
    ok(F(m + 1, m) * F(m + 2, m + 1) ** 2 * F(2 * m + 3, 2 * m + 1) ** 2 * F(8, 9) < 1, "tau term decreasing")
    ok(F(m + 1, m) ** 2 * F(16, 47) < 1, "secondary term decreasing")
    ok(F(237, 1000) + F(1, 1000) <= F(238, 1000), "|S_{m+1} - S_m| < 0.238 alpha")
    print("Proposition 6.6 (small angles): ok")


# ---------- Section 7: the phase ----------

def phase():
    # arcsin bounds via sin(x) > r  =>  arcsin(r) < x
    for r, x in ((F(94, 1000), F(942, 10000)), (F(176, 1000), F(177, 1000)), (F(57, 100), F(607, 1000))):
        lo, _ = sin_bounds(x)
        ok(lo > r, f"arcsin({r}) < {x}")
    ok(2 * F(942, 10000) < F(19, 100), "middle: 2 arcsin 0.094 < 0.19")
    ok(F(1, 4) - F(19, 100) >= F(6, 100) and F(3, 2) + F(19, 100) <= F(169, 100) < PI_LO, "middle gap in (0.06, 1.69)")
    ok(2 * F(607, 1000) < F(122, 100), "large: 2 arcsin 0.57 < 1.22")
    ok(F(3, 2) - F(122, 100) >= F(28, 100) and PI_HI / 2 + F(122, 100) <= F(280, 100) < PI_LO, "large gap in (0.28, 2.80)")
    ok(F(238, 1000) / F(824, 1000) < F(29, 100), "small: 0.238/0.824 < 0.29")
    # Lemma 7.3 (alpha = 3/m)
    m = 421
    ok(4 * F(16, 47) ** 2 < F(47, 100), "Arg H < 0.47 alpha")
    ok(F(4, 5) / F(1, 5) == 4, "-Arg(1-uq) <= 4 alpha")
    ok(F(1, 2) + F(47, 200) + 4 <= F(475, 100), "phi <= -pi/4 + 4.75 alpha")
    lowW = 3 + PI_LO / 4 - F(1425, 100) / m - F(177, 1000)
    highW = 3 + PI_HI / 4 + F(177, 1000)
    ok(lowW > F(357, 100) > PI_HI and highW < F(397, 100), "3.57 < W_m < 3.97")
    ok(F(129, 100) * 3 / m < F(1, 100) and F(397, 100) + F(1, 100) < 3 * PI_LO / 2, "W_{m+1} < 3 pi/2")
    # Lemma 7.4 (alpha = alpha_*)
    x = F(2, m)
    asin_hi = x / sqrt_lo(1 - x * x)  # arcsin x <= x / sqrt(1 - x^2)
    ok(asin_hi < F(5, 1000), "arcsin(2/m) < 0.005")
    ok(PI_HI / 8 + PI_HI - (PI_LO / 2 - F(5, 1000)) < F(197, 100), "phi < 1.97")
    ok((m + F(1, 2)) * asin_hi < F(201, 100), "(m+1/2) arcsin(2/m) < 2.01")
    ok(F(201, 100) + F(197, 100) + F(61, 100) == F(459, 100) and F(607, 1000) < F(61, 100), "total 4.59")
    ok(3 * PI_LO / 4 - F(459, 100) > -PI_LO, "n = 2k case")
    ok(5 * PI_LO / 4 - F(459, 100) > -PI_HI / 2, "n = 2k+1 case")
    print("Section 7 (phase): ok")


if __name__ == "__main__":
    micro()
    middle()
    large()
    small()
    phase()
    print(f"all constants ok ({count} inequalities)")
