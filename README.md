# MantiQ Factor: a computer-aided audit of Bagno-Garber (arXiv:2606.13149)

This repository contains a complete, reproducible verification workflow for the
paper **"Touchard-Riordan Polynomials and Schur-positivity of Set Partitions"**
by E. Bagno and D. Garber ([arXiv:2606.13149](https://arxiv.org/abs/2606.13149)),
run with the MantiQ Factor pipeline: **Wolfram Mathematica** as the symbolic
workbench, **Python** as an independent cross-check, and **Lean 4** as a formal
guard. The paper is not included here (link above); nothing else is needed to
reproduce everything.

## Headline results

1. **Theorem 3.1 of the paper is false in general.** Smallest counterexample:
   `(n, ell) = (4, 2)`, where the coefficient of the Schur function `s_(4)` in
   the paper's central sum is `q`, while the paper's formula gives
   `T_{3,2}(q) = 0` (there are no matchings on 3 points with 4 chords).
2. **Failure census** (verified in Wolfram to n = 6, Python to n = 7, Lean to
   n = 5): the identity holds exactly when `ell = 1`, `ell = n` or
   `n - ell = 1`, and fails exactly when `ell >= 2 and n - ell >= 2`.
3. **Corrected theorem** (new): the true hook coefficients are
   `A(m, ell, q)`, the crossing polynomials of the **adjacency-free**
   partitions of `[m]` into `ell` blocks:

   ```
   Q_n(Set(n, ell)) = sum_{k=0}^{n-ell} A(n-k, ell, q) * s_(n-k, 1^k)
   ```

   This is manifestly Schur-positive with nonnegative q-coefficients.
4. **Collapse lemma** (new): for every descent set `D`,
   `c_D(n, ell) = A(n - |D|, ell, q)`: the descent set only shifts.
5. **Corollary 3.3 (3) and (4) of the paper are true**, with an explicit
   Motzkin triangle identification: `A^NC(m, ell) = NCMatchCount(m-1, m-ell)`,
   row sums are Motzkin numbers, hook totals are Narayana numbers.
6. **Two new identities**: inclusion-exclusion
   `a(m, ell) = sum_j (-1)^j C(m-1, j) S(m-j, ell)` and
   `sum_ell a(m, ell) = Bell(m-1)`.

Full reasoning: [`mantiq-project/docs/plain_english_proof.md`](mantiq-project/docs/plain_english_proof.md).
Claim-by-claim ledger: [`mantiq-project/docs/formalization_map.md`](mantiq-project/docs/formalization_map.md).

## Repository layout

```
mantiq-project/
  src/SetPartitionQSym.wl      the Wolfram package (all definitions + audits)
  TouchardRiordanNotebook.nb   teaching notebook with interactive Manipulates
  build_project.wls            full build: suite + audit + figures + reports
  proof_audit.wls              audit only (fast)
  test_notebook.wls            notebook validation (parse + evaluate all cells)
  exports/                     generated reports, CSVs and figures
  docs/                        plain-English proof and formalization map
  tools/gen_notebook.py        generates the .nb (edit this, not the .nb)
prototype/
  verify3.py                   independent Python verification (n <= 7)
  verify3_results.json         its output for n = 7
lean/
  Verify.lean                  formal checks, pure Lean core (no Mathlib)
  lakefile.lean, lean-toolchain
```

## Prerequisites

- **Mathematica 14 or 15** with `wolframscript` on PATH (Wolfram Research).
- **Python 3** with nothing beyond the standard library.
- **Lean 4** via [elan](https://leanprover.github.io/lean4) (only for the
  formal lane; `lean-toolchain` pins the version).

## Quickstart

```text
:: 1. Wolfram lane (about 20 s): suite, audit, figures, reports
wolframscript -file mantiq-project\build_project.wls

:: 2. Teaching notebook: open in Mathematica and Evaluate Notebook
mantiq-project\TouchardRiordanNotebook.nb

:: 3. Independent Python cross-check (n <= 7)
python prototype\verify3.py 7

:: 4. Formal Lean lane
cd lean
lake build
```

Expected outcome: the Wolfram build prints a validation suite (all PASS) and a
10-check proof audit with one intentional WARN documenting the counterexample;
`verify3.py` prints V1-V9 all PASS; `lake build` completes with all theorems
proved by evaluation.

## The teaching notebook

`TouchardRiordanNotebook.nb` is built for classroom use: it loads the package,
walks through set partitions, descent sets and crossings with an interactive
arc-diagram explorer, states and tests Theorem 3.1, exhibits the counterexample
with a hook-by-hook comparison explorer, proves the collapse lemma, states the
corrected theorem, tabulates `A(m, ell, q)`, explores the Motzkin/Narayana
noncrossing story interactively, and closes with key findings. All interactive
controls are plain `Manipulate` objects; no extra packages are required.

## How the workflow runs (MantiQ Factor, 6 steps)

1. **Transcribe** the paper's definitions 1:1 into a Wolfram package.
2. **Compute** both sides of every claimed identity exactly and diff them.
3. **Census** the parameter space; find the smallest failing instance.
4. **Cross-check** in a second language (Python) before believing that a
   published theorem is false.
5. **Correct and strengthen**: state the theorem that is actually true, verify
   it, and derive the structural lemma (collapse) that explains the failure.
6. **Formalize** the decisive small cases in Lean 4 so the correction cannot
   quietly rot.

## Status of the failure pattern

The characterization "fails exactly when `ell >= 2` and `n - ell >= 2`" is
verified up to `n = 7` in Python and Wolfram and up to `n = 5` in Lean; beyond
that it is a conjecture. Proving it in general is a natural follow-up problem,
and the collapse lemma reduces it to comparing `A(m, ell, q)` with
Touchard-Riordan values.
