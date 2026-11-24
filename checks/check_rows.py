"""Check the closed formula of Lemma 2.2 for the rows C_n(t).

Compares it with a direct count of Catalan words by short peaks for n <= 12 and
with the power series of the generating function (Lemma 2.1) for n <= 60, and
prints the first rows (Table 1). Expected output: the table and two "ok" lines.
"""
from rows import row, brute, series_rows

for n in range(1, 9):
    print(f"C_{n}(t) coefficients: {row(n)}")
for n in range(1, 13):
    assert row(n) == brute(n), n
print("closed formula = direct count of Catalan words for 1 <= n <= 12: ok")
S = series_rows(60)
for n in range(0, 61):
    r = [c for c in S[n]]
    while len(r) > 1 and r[-1] == 0:
        r.pop()
    expect = [1] if n == 0 else row(n)
    assert r == expect, n
print("closed formula = generating function for 0 <= n <= 60: ok")
