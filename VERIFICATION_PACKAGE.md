# Verification Package: Theorem 3.1 of Bagno-Garber (arXiv:2606.13149) is False

This package is written for an independent verifier (human or agent). It contains
the exact statement being challenged, the definitions needed to interpret it, a
smallest counterexample worked out by hand, and four independent machine checks
with their expected outputs. Nothing here requires trusting the authors of this
audit: every claim can be re-derived from the definitions alone, and every
command below reruns in seconds to minutes on a laptop.

Paper: Eli Bagno and David Garber, "Touchard-Riordan Polynomials and
Schur-positivity of Set Partitions", arXiv:2606.13149v1.
The challenge: **Theorem 3.1, Equation (1), is false as stated.**

---

## 1. The statement being challenged (quoted from the paper)

Section 3 of the paper defines, for a set partition `pi` of `[n]`:

    Des(pi) = { i in [n-1] : i and i+1 are in the same block }

and states:

> **Theorem 3.1.**
> sum over `pi in Set(n, ell)` of `q^cr(pi) F_{n,Des(pi)}`
>   = sum over `k = 0 .. n-1` of `T_{n-1-k, n-k-ell}(q) s_{(n-k, 1^k)}`    (1)

Supporting definitions from the paper:

- **Definition 2.3 (F).** For `D` a subset of `[n-1]`,
  `F_{n,D}(x) = sum over weakly increasing i_1 <= ... <= i_n in [n]` with
  `i_j < i_{j+1}` for every `j in D`, of `x_{i_1} ... x_{i_n}`.
  (Rises are REQUIRED at the positions of `D`; extra rises are allowed. This is
  Gessel's fundamental quasisymmetric function, Section 2 of the paper.)

- **Linear graph and arcs (Section 2, before Example 2.5).** A partition `pi`
  of `[n]` is drawn as `n` points `1..n` on a horizontal line. For every block
  `B = {j_1 < ... < j_t}`, every two **consecutive** elements `j_k, j_{k+1}` of
  the block are joined by a semicircular arc; when `j_{k+1} = j_k + 1` they are
  joined by a straight line instead of an arc.

- **Definition 2.6 (cr).** For `r < s`, two arcs `(i_r, j_r)` and `(i_s, j_s)`
  cross iff `i_r < i_s < j_r < j_s`. `cr(pi)` is the number of such crossings.

- **Section 2.2 (T).** `T_{2j}(q) = sum over perfect matchings M on 2j points
  of q^cr(M)` is the Touchard-Riordan polynomial (j chords), and
  `T_{n,j}(q) = C(n, 2j) T_{2j}(q)` extends it to matchings with unmatched
  points (choose which 2j of the n points carry the chords).

**Convention note (this closes the only plausible escape route).** A straight
segment joining `i` and `i+1` can never participate in a crossing under
Definition 2.6, because no integer lies strictly between `i` and `i+1`: a
crossing with `(i, i+1)` as the outer arc needs an inner endpoint strictly
between them. So counting only the arcs of length >= 2 (and giving straight
segments zero crossings) is *equivalent* to applying Definition 2.6 literally
to every pair. The paper's own Example 2.5 (`pi = {{1,3,7},{2,4,5,6,8}}`,
`cr(pi) = 3`) is reproduced by both conventions; the checks below reproduce it.

---

## 2. The smallest counterexample: `(n, ell) = (4, 2)`

`Set(4,2)` has exactly 7 partitions. Their `Des` and `cr` (hand-checkable;
also printed by every script below):

| partition            | Des(pi)  | cr(pi) | arcs (length >= 2) |
|----------------------|----------|--------|--------------------|
| {{1,2,3},{4}}        | {1,2}    | 0      | none               |
| {{1,2,4},{3}}        | {1}      | 0      | (1,4)              |
| {{1,2},{3,4}}        | {1,3}    | 0      | none               |
| {{1,3,4},{2}}        | {3}      | 0      | (1,4)              |
| {{1,3},{2,4}}        | {}       | **1**  | (1,3), (2,4) cross |
| {{1,4},{2,3}}        | {2}      | 0      | (1,4)              |
| {{1},{2,3,4}}        | {2,3}    | 0      | (2,4)              |

Only the crossing partition `{{1,3},{2,4}}` contributes a factor `q`.

The left-hand side of (1) at `(4,2)` is therefore the symmetric function

    Q_4(Set(4,2)) = q * s_(4) + s_(3,1) + s_(2,1,1)

(exact Schur expansion; derivation below and in both scripts).

The right-hand side of (1) at `(4,2)` has three terms:

    k=0:  T_{3,2}(q) s_(4)  =  C(3, 4) T_4(q) s_(4)  =  0        <-- the vanishing term
    k=1:  T_{2,1}(q) s_(3,1) = C(2, 2) T_2(q) s_(3,1) = 1 * s_(3,1)
    k=2:  T_{1,0}(q) s_(2,1,1) = C(1, 0) T_0(q) s_(2,1,1) = 1 * s_(2,1,1)

so the paper's RHS is `s_(3,1) + s_(2,1,1)` and is **missing the entire
`q * s_(4)` term**:

    LHS - paper RHS = q * s_(4) = q * h_4.

Equivalently, in the monomial basis (4 variables), the paper's RHS equals
`h_2 e_2`, while the true LHS equals `q h_4 + h_2 e_2`.

### How the exact expansion is obtained (so the verifier can redo it)

Build `Q_4(Set(4,2))` from the definitions:

1. List the 7 partitions with their `Des` and `cr` (table above).
2. For each exponent tuple `alpha = (a_1,a_2,a_3,a_4)` with `sum a_i = 4`,
   the corresponding weakly increasing sequence is the sorted expansion of
   `alpha`; its rise set is `R(alpha) = { j : sorted(seq)_{j-1} < sorted(seq)_j }`.
   The coefficient of `x_1^{a_1} x_2^{a_2} x_3^{a_3} x_4^{a_4}` in
   `F_{4,D}` is 1 iff `D` is a subset of `R(alpha)`, else 0.
3. So the coefficient of that monomial in `Q_4(Set(4,2))` at q-degree `d` is
   the number of partitions `pi` with `cr(pi) = d` and `Des(pi) subset of R(alpha)`.

The resulting monomial coefficients (constant on permutations of the
exponents; the function is symmetric, which the scripts assert directly):

| weight mu (as partition) | coeff at q^0 | coeff at q^1 |
|--------------------------|--------------|--------------|
| (4)                      | 0            | 1            |
| (3,1)                    | 1            | 1            |
| (2,2)                    | 1            | 1            |
| (2,1,1)                  | 3            | 1            |
| (1,1,1,1)                | 6            | 1            |

Convert monomial coefficients to Schur coefficients with the Kostka matrix
`s_lambda = sum_mu K_{lambda,mu} m_mu` (K computed by enumerating semistandard
Young tableaux, so nothing is trusted from literature):

    K (rows = lambda, cols = mu in the order (4), (3,1), (2,2), (2,1,1), (1^4)):
    (4):      [1, 1, 1, 1, 1]
    (3,1):    [0, 1, 1, 2, 3]
    (2,2):    [0, 0, 1, 1, 2]
    (2,1,1):  [0, 0, 0, 1, 3]
    (1^4):    [0, 0, 0, 0, 1]

Solving `K^T c = f` exactly (fractions) gives:

    q^0: c_(4) = 0, c_(3,1) = 1, c_(2,2) = 0, c_(2,1,1) = 1, c_(1^4) = 0
    q^1: c_(4) = 1, all others 0

i.e. `Q_4(Set(4,2)) = q * s_(4) + s_(3,1) + s_(2,1,1)`. Hand-check of the
nonzero q^0 coefficients from the Kostka rows:
`s_(3,1) = m_(3,1) + m_(2,2) + 2 m_(2,1,1) + 3 m_(1^4)` and
`s_(2,1,1) = m_(2,1,1) + 3 m_(1^4)`; their sum has monomial coefficients
`(4): 0, (3,1): 1, (2,2): 1, (2,1,1): 3, (1^4): 6`, exactly the q^0 column of
the table above; the scripts verify this algebraically instead of trusting
the hand expansion.

---

## 3. Reproduction commands and expected output

Run all commands from the repository root
(`https://github.com/zuwasi/mantiq-touchard-riordan`). Requirements: Python 3.8+
(standard library only) for checks 1 and 2; a Mathematica installation with
`wolframscript` for check 3; Lean 4 with `elan`/`lake` for check 4.

### Check 1 (pure Python, ~2 seconds, no dependencies): the counterexample

    python prototype\standalone_check_4_2.py

The script rebuilds `Q_4(Set(4,2))` from scratch (7 partitions, weakly
increasing sequences, arcs, crossings), asserts symmetry, computes the Kostka
matrix by enumerating SSYTs, inverts it exactly with Fractions, and asserts
the result. Expected final output:

    Schur coefficients at q^0: {(3, 1): 1, (2, 1, 1): 1}
    Schur coefficients at q^1: {(4,): 1}

    VERIFIED: Q_4(Set(4,2)) = q*s_(4) + s_(3,1) + s_(2,1,1).
    Paper Theorem 3.1 RHS for (4,2) = s_(3,1) + s_(2,1,1) (T_{3,2} = 0),
    so the theorem is FALSE: it misses the q*s_(4) term.
    Equivalently, LHS - paper RHS = q*h_4 = q*s_(4).

It also prints the 7-partition table and the per-monomial coefficients; the
table must match Section 2 above.

### Check 2 (pure Python, seconds): the full census

    python prototype\verify3.py 7

Expected key lines (full output also saved to `verify3_results.json`):

    V2 hook SYT descent sets = all k-subsets: PASS (n<=7)
    V3 collapse c_D = A(n-|D|,ell) for every D: PASS (n<=7)
    V1 Gessel Thm 2.4 in M-basis (#{T: Des<=pd(beta)} = K_(lam,mu)): PASS (n<=7)
    V4 corrected theorem f = sum_k A(n-k,ell,q) s_hook: PASS (n<=7)
    V5 paper Theorem 3.1 (T-coeffs): 10 of 28 (n,ell) pairs FAIL (n<=7)
    V5 first differing (n,ell,k): true A vs paper T:
        (n=4, ell=2, k=0): A(4,2,q) = 1q^1   vs   T_{3,2}(q) = 0
        (n=5, ell=2, k=0): A(5,2,q) = 1q^2   vs   T_{4,3}(q) = 0
        (n=5, ell=2, k=1): A(4,2,q) = 1q^1   vs   T_{3,2}(q) = 0
        (n=5, ell=3, k=0): A(5,3,q) = 2 + 5q^1   vs   T_{4,2}(q) = 2 + 1q^1
        (n=6, ell=2, k=0): A(6,2,q) = 1q^3   vs   T_{5,4}(q) = 0
        (n=6, ell=2, k=1): A(5,2,q) = 1q^2   vs   T_{4,3}(q) = 0
        (n=6, ell=2, k=2): A(4,2,q) = 1q^1   vs   T_{3,2}(q) = 0
        (n=6, ell=3, k=0): A(6,3,q) = 5q^1 + 9q^2 + 1q^3   vs   T_{5,3}(q) = 0
    V6 Kostka-solve Schur expansion == corrected hook formula: PASS (n<=6)
    V7a NC collapse c^NC_D = A^NC(n-|D|,ell) for every D: PASS (n<=7)
    V7b NC total = sum_k C(n-1,k) A^NC(n-k,ell) = Narayana(n,ell): PASS (n<=7)
    V7c A^NC(m,ell) = #NC matchings on m-1 pts with m-ell chords = C(m-1,2b)Cat_b: PASS (n<=7)
    V7d sum_ell A^NC(m,ell) = Motzkin(m-1): PASS (n<=7)
    V8 inclusion-exclusion A(m,ell,q) = sum_j (-1)^j C(m-1,j) H(m-j,ell,q): PASS (n<=7)
    V9 a(m,ell)=sum_j (-1)^j C(m-1,j) S(m-j,ell) and sum_ell a(m,ell)=Bell(m-1): PASS

Interpretation: the paper's Theorem 3.1 agrees exactly on the pairs with
`ell = 1`, `ell = n`, or `n - ell = 1`, and fails for every pair with
`ell >= 2` and `n - ell >= 2` (10 of the 28 pairs with n <= 7). V6
independently re-derives the Schur expansion through a second path (Kostka
solve) inside the same script.

### Check 3 (Wolfram Mathematica, about 20 seconds): symbolic algebra

    wolframscript -file prototype\wl_check_4_2.wls

Requires the Wolfram package `mantiq-project/src/SetPartitionQSym.wl` (loaded
by the script). Expected key lines:

    LHS == CorrectedRHS: True
    LHS - PaperRHS = q*(...)        (exactly q times the complete homogeneous
                                     function h_4: every monomial of degree 4
                                     appears with coefficient 1)
    pi = {{1, 3}, {2, 4}}  Des = {}  cr = 1     (the only crossing partition)
    A(4,2,q) = q
    T_(3,2)(q) = 0

The printed `Q_4(Set(4,2))` polynomial is the direct symbolic sum over the 7
partitions; its non-q part has no `x1^4, x2^4, x3^4, x4^4` terms (the (4)
weight orbit has q-degree-0 coefficient 0), and its q^0 coefficients on the
remaining weight orbits (3,1)/(2,2)/(2,1,1)/(1^4) of degree 4 are
1 / 1 / 3 / 6, matching the table in Section 2.

A broader Wolfram audit (11 checks A1-A11, including the corrected theorem to
n = 6 and the exact failure census to n = 6):

    wolframscript -file mantiq-project\proof_audit.wls

Expected summary line: `== audit summary: 10 passed, 0 failed, 1 warnings ==`
(the single WARN is check A3, which flags the counterexample itself; it is a
finding, not a failure of the audit).

### Check 4 (Lean 4, under 2 seconds): formal definitions and census to n = 5

    cd lean
    lake build

Expected: `Build completed successfully` (4 jobs). What `lean/Verify.lean`
formalizes and proves (as finite exhaustive checks): set partitions as
restricted growth strings with descent sets and arc crossings; the corrected
coefficient `A(m, ell, q)` and the paper's `T_{a,b}(q)` (the latter by
enumerating all matchings); the counterexample (`A(4,2) = q` nonzero while
`T_{3,2} = 0`); the failure census, exactly `(4,2)`, `(5,2)`, `(5,3)` within
n <= 5; the collapse lemma for every descent set at `(4,2)` and `(5,3)` and
for every `(n, ell)` with n <= 4; and the Bell, inclusion-exclusion, Motzkin
triangle, Motzkin row-sum, and Narayana identities. It does not formalize the
fundamental quasisymmetric functions or the Schur expansion itself, so it
checks the census and the coefficient identities, not a re-derivation of the
Schur expansion. Proof style: 11 of the 12 theorems are `native_decide`
(trusted compiler evaluation of an exhaustive finite enumeration), and
`T_3_2` is `rfl`; this is a finite cross-check of the same definitions, not a
kernel-level proof about all n. To force recompilation from scratch:

    cd lean
    lake env lean Verify.lean

(exits 0 with no output when everything checks).

---

## 4. What is true instead (the corrected theorem)

For all `(n, ell)` with `1 <= ell <= n` (corrected identity and full failure
census: Python to n <= 7, Wolfram to n <= 6; Lean 4 verifies the census of
failing pairs and the collapse and coefficient identities in the same finite
range, to n <= 5):

    Q_n(Set(n, ell)) = sum over k of A(n-k, ell, q) * s_(n-k, 1^k)

where `A(m, ell, q)` is the crossing generating polynomial of the
**adjacency-free** partitions of `[m]` with `ell` blocks (partitions with
empty descent set). The engine behind it is the collapse lemma

    c_D(n, ell) = A(n - |D|, ell, q)   for every descent set D,

where `c_D(n, ell)` is the crossing-number distribution of the partitions of
`[n]` with `ell` blocks whose descent set is **exactly** `D` (for each `cr`
value, the number of such partitions with that crossing number). The lemma
says this distribution depends only on `|D|`: it equals the distribution of
the adjacency-free partitions of `[n - |D|]` with `ell` blocks. The intuition
is that each maximal run of descents can be collapsed by deleting the upper
endpoint of each descent adjacency and relabeling, which preserves the arc
structure and hence `cr`. This is Check V3 (Python, every descent set for
every `(n, ell)` with n <= 7); Wolfram check A10 is the spot check at
`(n, ell) = (4,2)` for `D = {}, {1}, {2}, {3}`, and Lean proves it for every
descent set at `(4,2)` and `(5,3)` and for every `(n, ell)` with n <= 4.

Consequences, all verified:

- The corrected identity is manifestly Schur-positive (every `A` is a
  polynomial with nonnegative integer coefficients).
- Corollary 3.3 (3) and (4) of the paper remain TRUE: the noncrossing
  coefficients are the Motzkin triangle numbers (V7c: `A^NC(m, ell) =
  C(m-1, 2b) * Cat_b` where `b = m - ell`), row sums are Motzkin numbers
  (V7d), and the hook totals are the Narayana numbers (V7b).
- New identities at q = 1 (V8, V9):
  `a(m, ell) = sum_j (-1)^j C(m-1, j) S(m-j, ell)` and
  `sum over ell of a(m, ell) = Bell(m-1)`, i.e. adjacency-free partitions of
  `[m]` (all block sizes) are counted by the Bell number `Bell(m-1)`.

The reason the paper's proof breaks: the bijection it constructs from pairs
(SYT, chord diagram) to partitions is injective but **not surjective**: chord
diagrams whose chords would force forbidden adjacencies (consecutive integers
in the same block arising at positions that are not descents of the SYT) are
never hit. The fibers are exactly the descent runs, which is what the collapse
lemma measures.

---

## 5. Verifier checklist

Work through these in order; each is independently checkable.

1. **Definitions.** Confirm from the paper (arXiv:2606.13149v1, PDF included
   in this repository as `2606.13149v1.pdf`): Definition 2.3 (F, rises
   required at D, extra rises allowed), the linear-graph construction and
   Definition 2.6 (cr), Section 2.2 (T_{n,j}(q) = C(n,2j) T_{2j}(q)),
   Section 3 (Des(pi) = same-block adjacency set; Theorem 3.1, Equation (1)).
2. **Straight segments cannot cross** (Section 1 convention note). Convince
   yourself this is airtight; it removes the only ambiguity in `cr`.
3. **Hand table.** Verify the 7 partitions of `Set(4,2)`, their `Des` sets,
   and their `cr` values (only `{{1,3},{2,4}}` crosses, with `cr = 1`).
4. **Hand monomial table.** For a couple of weights, e.g. `(3,1)` and `(2,1,1)`,
   count by hand which partitions have `Des subset of R(alpha)` and check the
   coefficients 1/1 (weight (3,1)) and 3/1 (weight (2,1,1)) at q^0/q^1.
   For example, for weight `(2,1,1)` (monomials like `x1^2 x2 x3`): the rise
   set is `{2,3}` for the sequence (1,1,2,3), so the contributing partitions
   are those with `Des subset of {2,3}`: `{{1,4},{2,3}}`, `{{1,3,4},{2}}`,
   `{{1},{2,3,4}}` (cr = 0, three of them) and `{{1,3},{2,4}}` (cr = 1):
   coefficients 3 at q^0 and 1 at q^1.
5. **Run Check 1** (`standalone_check_4_2.py`) and read its source: it is
   ~150 lines, uses only the standard library, and every definition is
   written out literally. The verdict asserts fail loudly on any mismatch.
6. **Run Check 2** (`verify3.py 7`): the census and the corrected theorem.
   Note that V5 first differs exactly at `(4,2)` with `A(4,2,q) = q` vs
   `T_{3,2}(q) = 0`.
7. **Run Check 3** (`wl_check_4_2.wls`) if you have Mathematica: the symbolic
   difference `LHS - PaperRHS` is exactly `q * h_4`.
8. **Run Check 4** (Lean) if you have Lean 4: `lake build` compiles the
   formalization.
9. **Sanity-check the escape hatches.** (a) Could `F` have been meant with
   rises *forbidden* outside D? No: Definition 2.3 says `i_j < i_{j+1}` if
   `j in D`, and the paper cites Gessel; the check scripts implement exactly
   the quoted condition. (b) Could `T_{n,j}` be meant differently? Section
   2.2 of the paper defines `T_{n,j}(q) = C(n,2j) T_{2j}(q)` explicitly; with
   that definition `T_{3,2}(q) = C(3,4) T_4(q) = 0` regardless of `T_4`. (c)
   Could `cr` count something else? Definition 2.6 is quoted above, and the
   paper's own Example 2.5 (`cr = 3`) is reproduced by the scripts.

**Bottom line.** The coefficient of `s_(4)` in `Q_4(Set(4,2))` is `q`, but
Theorem 3.1 assigns it `T_{3,2}(q) = 0`. Theorem 3.1 is false as stated; the
corrected statement, the census of all failing parameter pairs, and three
new positive identities are documented above and verified by machine.
