#!/usr/bin/env bash
# Runs every exact check for the paper and writes checks/output.txt. Needs Python 3;
# set PYTHON to an interpreter with gmpy2 installed to speed up the last step.
set -euo pipefail
cd "$(dirname "$0")"
{
  echo "== Lemma 2.2: closed formula for the rows"
  "${PYTHON:-python3}" check_rows.py
  echo
  echo "== Sections 4-7: numerical inequalities, in exact arithmetic"
  "${PYTHON:-python3}" constants.py
  echo
  echo "== Section 9: interlacing of C_n and C_(n+1) for 1 <= n <= 420"
  "${PYTHON:-python3}" finite_interlacing.py 420
  echo
  echo "all checks passed"
} | tee output.txt
