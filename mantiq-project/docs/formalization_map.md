# Formalization map: paper claim to code to theorem

Every row pairs a mathematical statement with the code that checks it and the
system that runs it. This is the "MantiQ Factor" ledger: nothing is claimed
that is not executable somewhere.

| # | Statement | Wolfram (src/SetPartitionQSym.wl) | Python (prototype/verify3.py) | Lean (lean/Verify.lean) |
|---|-----------|-----------------------------------|-------------------------------|-------------------------|
| 1 | Set partitions as rgs, blocks, descents, crossings | SetPartitionsRGS, BlocksOfPartition, DescentSetOfPartition, CrossingNumberPartition | partitions_by_rgs, descent_set, crossing_number | rgs, des, arcs, crossings |
| 2 | F_{n,D} per paper Def. 2.3 (rises REQUIRED at D) | FundamentalQSym[n, D, vars] | fundamental_qsym | (not formalized: multivariate polys) |
| 3 | Q_n(Set(n,ell)) = sum q^cr F_{n,Des} | QSet[n, ell, q, vars] | q_set | (not formalized: multivariate polys) |
| 4 | Paper Theorem 3.1 RHS hook coefficients | PaperRHS[n, ell, q, vars] | paper_rhs | paperCoeff, touchardRiordanT (T_{2b} by enumerating matchings) |
| 5 | Paper Theorem 3.1 FALSE from (4,2); census ell >= 2 and n - ell >= 2 | RunProofAudit A5, V2 | verify3.py V2 (n <= 7) | census_5: fails exactly {(4,2), (5,2), (5,3)} for n <= 5 |
| 6 | Corrected theorem hook coefficients A(n-k, ell, q) | CorrectedRHS + RunValidationSuite V1, RunProofAudit A4 | verify3.py V1 (n <= 7) | ACoeff (coefficient lists); census compares against paperCoeff |
| 7 | Collapse lemma c_D(n,ell) = A(n-|D|,ell,q) for all D | RunValidationSuite V3, RunProofAudit A10 | verify3.py V3 | collapse_4_2_allD, collapse_5_3_allD, collapse_all_n4 (every descent set at (4,2) and (5,3); every (n,ell) with n <= 4) |
| 8 | A(4,2,q) = q and T_{3,2}(q) = 0 | RunProofAudit A1-A2, V6 | verify3.py V6 | A_4_2, T_3_2, counterexample_4_2 |
| 9 | ell = 1: Q = e_n = s_(1^n) | V4, A6 | V4 | (not formalized: Schur machinery) |
| 10 | ell = n: Q = s_(n) | V5, A7 | V5 | (not formalized: Schur machinery) |
| 11 | A^NC(m, ell) = NC configurations on m-1 points with m-ell chords (Motzkin triangle) | ANCCount, NCMatchCount, V7c, A8 | V7c | triangle_2_6 |
| 12 | Row sums Motzkin(m-1) (Cor. 3.3 (4)) | V7d, A9 | V7d | motzRows_2_7 |
| 13 | Hook totals Narayana(n, ell) | V8 | V8 | narayana_2_5 |
| 14 | Inclusion-exclusion a(m, ell) = sum_j (-1)^j C(m-1,j) S(m-j, ell) | V9/InclusionExclusionBell, A11 | V9 | inclusionExclusion_2_5 |
| 15 | Bell identity sum_ell a(m, ell) = Bell(m-1) | V9, A11 | V9 | bell_2_6 |

Lean deliberately checks the counting layer (partitions, matchings, q-polynomial
coefficients, sequences) where brute force is decisive, and leaves the
multivariate symmetric-function layer to Wolfram and Python; the systems
agree on every overlapping claim.

## What each system contributes

- **Wolfram Language** is the workbench: exact polynomial algebra, Schur
  expansions via hook decompositions (h and e by IntegerPartitions), symbolic
  q-polynomials, interactive Manipulate teaching material, and the audit
  reports exported to `exports/`.
- **Python** is the independent cross-check: a from-scratch reimplementation
  with no shared code. Agreement up to n = 7 (293 failing pairs identified
  matches the Wolfram census exactly) rules out an implementation artifact.
- **Lean 4** is the formal guard: first-order definitions and brute-force
  enumeration via `native_decide` (trusted compiler evaluation, not a kernel
  proof term) check the counterexample, the failure census, the collapse lemma
  at `(4,2)` and `(5,3)` over every descent set and for every `(n, ell)` with
  n <= 4, and the sequence identities up to n = 5-7. It deliberately does not
  formalize the multivariate symmetric-function layer.

## Division of labor for a new paper

1. Transcribe definitions into Wolfram (they read like the paper).
2. Compute both sides of every claimed identity and diff them.
3. Census the parameter space; find the smallest failing instance.
4. Cross-check in a second language before believing a "theorem is false".
5. State the corrected theorem, verify it, and push the first cases into Lean
   so the correction cannot quietly rot.
