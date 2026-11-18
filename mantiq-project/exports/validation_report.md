# Validation report -- MantiQ Factor reproduction of Bagno-Garber (arXiv:2606.13149)

Wolfram version: 15.0.1 for Microsoft Windows (64-bit) (July 2, 2026)

* CorrectedTheorem: PASS
* PaperTheorem31: FAIL -- {{4, 2}, {5, 2}, {5, 3}, {6, 2}, {6, 3}, {6, 4}}
* CollapseLemma: PASS
* HookSYTDescentSets: PASS
* NCMotzkinTriangle: PASS
* NCTotals: PASS
* InclusionExclusionBell: PASS

## Headline results

- Theorem 3.1 of the paper is false from (n,l) = (4,2) on: the coefficient
  of s_(n-k,1^k) is the crossing polynomial A(n-k,l,q) of adjacency-free
  partitions, not the Touchard-Riordan value T_(n-1-k,n-k-l)(q).
- Corollary 3.3 (3)-(4) (noncrossing) is correct: the coefficients are the
  Motzkin triangle numbers NCMatchCount(m-1, m-l), and the block-count
  totals are Narayana(n,l) and Motzkin numbers.
- New identities verified: collapse lemma c_D = A(n-|D|,l,q) for every
  descent set D; inclusion-exclusion for A(m,l,q); Sum_l a(m,l) = Bell(m-1).
