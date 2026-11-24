# Short peaks in Catalan words: real zeros and interlacing

Brandon Li

Paper and verification code proving a conjecture of Alexandersson on the
polynomials $C_n(t)$ that count Catalan words of length $n$ by peaks (short
peaks in the terminology of Baril, Harris, Harry, McClinton and Ramírez). The
same triangle counts Dyck paths of semilength $n$ by occurrences of $UDUU$
([OEIS A116424](https://oeis.org/A116424)). The conjecture is listed among the
[open problems on real-rootedness](https://www.symmetricfunctions.com/realRootedOpenProblems.htm)
of the symmetric functions catalog.

**Theorem.** Every $C_n(t)$ has only simple, negative real zeros. For every
$n\ge1$ the zeros of $C_n$ and $C_{n+1}$ strictly interlace, and the largest
zero among them belongs to $C_{n+1}$.

For $n\ge420$ the proof is analytic. With $a=1-t$, a contour integral around
the four branch cuts of the generating function writes each row as an explicit
oscillating cosine, and the phases of consecutive rows differ by less than
$\pi$; together with a local argument near $a=27/16$, this produces
$\deg C_n+\deg C_{n+1}$ alternating sign changes. For $n\le420$ the theorem is
checked by an exact Sturm-sequence computation.

- [Paper (PDF)](paper/main.pdf)
- [LaTeX source](paper/)
- [Verification code](checks/)

## Reproduce

```sh
bash checks/run_all.sh
```

This needs Python 3 (exact integer and rational arithmetic). It runs:

- `check_rows.py`: Lemma 2.2, comparing the closed formula for the rows with a
  direct count of Catalan words ($n\le12$) and with the generating function
  ($n\le60$).
- `constants.py`: every numerical inequality of Sections 4–7, in exact
  rational arithmetic with rational bounds for $\pi$, square roots, sine,
  cosine and the exponential.
- `finite_interlacing.py`: the Sturm-sequence test of Section 9 for every pair
  $(C_n,C_{n+1})$ with $1\le n\le420$.

`rows.py` (the rows $C_n$ and the closed formula of Lemma 2.2) and
`polyarith.py` (exact polynomial arithmetic) are shared helpers.

The last step uses large integer gcds. With
[gmpy2](https://pypi.org/project/gmpy2/) installed it takes about twenty minutes; set
`PYTHON` to that interpreter, for example
`PYTHON=.venv/bin/python bash checks/run_all.sh`. Plain Python works but is
much slower. The output is in [`checks/output.txt`](checks/output.txt).

`numerics.py` (needs `mpmath`) evaluates both sides of the integral formula of
Proposition 5.1 and tabulates $|S_m-1|$ and the phase gaps; the proof does not
use it. The PDF builds with Tectonic (`cd paper && tectonic main.tex`) or
`pdflatex`.

## License

Code: MIT. Paper: CC BY 4.0.
