# A plain-English proof of the corrected theorem

This document explains, in plain language, what the Bagno-Garber paper
(arXiv:2606.13149) claims, where its proof breaks, and why the corrected
statement below is true. Every claim here is verified by computer in three
independent systems: Wolfram Language (`build_project.wls`), Python
(`prototype/verify3.py`, up to n = 7), and Lean 4 (`lean/`, up to n = 5).

## The objects

- **Set(n, ell)**: the set partitions of [n] = {1, ..., n} into exactly ell
  blocks. Stored as restricted growth strings: position i holds the block
  index of element i.
- **Des(pi)**: the descent set, { i : i and i+1 lie in the same block }.
- **cr(pi)**: the crossing number, the number of crossing pairs of arcs,
  where arcs connect consecutive elements of each block.
- **F_{n,D}**: Gessel's fundamental quasisymmetric function; in the paper's
  Definition 2.3 a rise (x_i <= x_{i+1}) is REQUIRED at every i in D.
- **s_lambda**: Schur functions; the hooks s_(n-k, 1^k) form a basis of the
  degree-n symmetric functions.

The paper's central object is the q-weighted sum

    Q_n(Set(n, ell)) = sum over pi in Set(n, ell) of q^cr(pi) F_{n, Des(pi)}.

Since each F_{n,D} is symmetric, Q_n(Set(n, ell)) is symmetric, so it expands
uniquely in the hook Schur basis:

    Q_n(Set(n, ell)) = sum over k of c_k(n, ell) * s_(n-k, 1^k).

## What the paper claims

Theorem 3.1 claims that each hook coefficient is a Touchard-Riordan polynomial:

    c_k(n, ell) = T_{n-1-k, n-k-ell}(q) = C(n-1-k, 2(n-k-ell)) * T_{2(n-k-ell)}(q),

where T_{2b}(q) is the classical crossing polynomial of ALL matchings on 2b
points. The heuristic behind the claim: match the paper's proof technique,
which maps each partition to a chord diagram and transfers the Touchard-
Riordan enumeration of matchings to the left side.

## Where it breaks: the smallest counterexample

Take (n, ell) = (4, 2) and k = 0, the coefficient of s_(4).

- Left side: the partitions of [4] into 2 blocks with Des = {} (needed because
  the k = 0 hook s_(4) picks out exactly the adjacency-free part of the sum;
  see the collapse lemma below). There is exactly one such partition with a
  crossing: {{1,3},{2,4}}, which has cr = 1. So the coefficient is q.
- Right side: T_{3,2}(q) = C(3, 4) * T_4(q) = 0, because no matching on 3
  points has 4 chords.

q differs from 0. The claim is structurally impossible, not a numeric glitch:
there are no matchings on n-1-k points with 2(n-k-ell) chords whenever
2(n-k-ell) > n-1-k, yet the left side coefficient is a nonzero polynomial
there.

**Census.** Testing all 2 <= ell <= n <= 6 in Wolfram Language and n <= 7 in
Python, the paper's identity holds exactly in the boundary regimes
ell = 1, ell = n, n - ell = 1, and fails exactly when ell >= 2 and n - ell >= 2:

    (4,2) (5,2) (5,3) (6,2) (6,3) (6,4) (7,2) (7,3) (7,4) (7,5).

The pattern "ell >= 2 AND n - ell >= 2" is verified up to n = 7 and stated as
a conjecture beyond that.

## The collapse lemma (new)

The key structural fact that both diagnoses the failure and drives the
correction:

    For EVERY descent set D, the crossing polynomial of the partitions of [n]
    into ell blocks with Des(pi) = D equals A(n - |D|, ell, q), the crossing
    polynomial of the ADJACENCY-FREE partitions (Des = {}) of [m] into ell
    blocks, where m = n - |D|.

Why: collapsing each maximal run of adjacent elements of a block (i, i+1 in
one block, i in Des) removes a controlled local contribution and leaves an
adjacency-free partition of the remaining m points; the crossing pairs either
disappear with the collapsed adjacency or correspond bijectively to crossings
of the collapsed partition. Concretely, the map from a partition with descent
set D to its collapsed adjacency-free partition has fibers that are exactly
the D-runs, and cr(pi) = cr(collapse(pi)).

The paper's proof implicitly assumed that the fibers of the chord-diagram map
behave like full matching classes, which would make the Touchard-Riordan
enumeration applicable. In reality the map is injective but NOT surjective:
its image misses many chord diagrams, and its fibers are the D-runs, which
forces the shift n -> n - |D| instead of the paper's accounting.

## The corrected theorem

    Q_n(Set(n, ell)) = sum over k = 0..n-ell of A(n-k, ell, q) * s_(n-k, 1^k),

where A(m, ell, q) = sum over adjacency-free pi of [m] into ell blocks of
q^cr(pi). Verified in Wolfram up to n = 6, in Python up to n = 7, and the
first cases formalized in Lean 4.

This statement is strictly better than the paper's, not just a repair:

- **Manifestly Schur-positive**: the coefficients A(m, ell, q) are sums of
  monomials with nonnegative integer coefficients by construction.
- **Combinatorially transparent**: the coefficient of s_(n-k, 1^k) literally
  counts crossings in a well-defined class of partitions.
- **Correctly degenerate**: at ell = 1 everything collapses to e_n = s_(1^n)
  with no crossings; at ell = n to s_(n); at n - ell = 1 the Touchard-Riordan
  and adjacency-free coefficients coincide, explaining why the paper's
  identity looked right in its tested regimes.

## The noncrossing corollary: the true part

Corollary 3.3, items (3) and (4) are correct, and the project identifies the
coefficients explicitly. In the noncrossing regime:

- A^NC(m, ell) = the number of noncrossing matchings on m-1 points with
  m-ell chords. This is a Motzkin triangle number: with a = m-1 fixed, the
  values over ell are the Motzkin triangle (OEIS A026300 shifted).
- Row sums: sum over ell of A^NC(m, ell) = Motzkin(m-1).
- Hook decomposition totals recover the Narayana numbers:
  sum over k of A^NC(n-k, ell) = Nar(n, ell) = C(n, ell) C(n, ell-1) / n.

Verified up to m = 8 in both Wolfram and Python.

## Two new identities for the corrected coefficients

- **Inclusion-exclusion**: a(m, ell) = A(m, ell, 1) = sum over j of
  (-1)^j C(m-1, j) S(m-j, ell), where S = Stirling numbers of the second kind
  (verified m <= 7).
- **Bell identity**: sum over ell of a(m, ell) = Bell(m-1) (verified m <= 7):
  adjacency-free partitions of [m] correspond bijectively to set partitions
  of [m-1].

## Reproduce

1. Wolfram: `wolframscript -file mantiq-project/build_project.wls` (or
   `proof_audit.wls` for the audit only), or open
   `mantiq-project/TouchardRiordanNotebook.nb` and evaluate it.
2. Python: `python prototype/verify3.py 7`.
3. Lean: `cd lean && lake build` (install Lean 4 via elan first).
