# Email to Dr. Eli Bagno

**To:** bagnoe@g.jct.ac.il
**Subject:** An independent computational audit of your Touchard-Riordan paper, with a corrected theorem and a full reproducible workflow

Dear Dr. Bagno,

I am Daniel Liezrowice, from the Engineering Software Lab (ESL). We develop an automated mathematical workflow we call the MantiQ Factor: Mathematica for exact computation, Python for an independent cross-check, and Lean 4 for a formal, machine-checked guard. We ran it end to end on your recent paper with Dr. Maria Garber, "Touchard-Riordan Polynomials and Schur-positivity of Set Partitions" (arXiv:2606.13149), and I would like to share the results with you, because they concern Theorem 3.1 directly.

**1. What we found**

Theorem 3.1, as stated, is false. The smallest counterexample is (n, ell) = (4, 2): the true coefficient of the Schur function s_(4) in the q-weighted partition sum Q_4(Set(4,2)) is the monomial q (it comes from exactly one of the seven two-block partitions, the crossing one), while the paper's formula T_(3,2)(q) evaluates to 0 because a binomial coefficient vanishes.

We then swept the whole parameter space. The identity holds exactly when ell = 1, when ell = n, and when n - ell = 1. It fails for every (n, ell) with ell >= 2 and n - ell >= 2, starting at (4, 2). Three independent systems agree on every overlapping case: Wolfram Mathematica verified n <= 6, a from-scratch Python verifier verified n <= 7, and a Lean 4 library verified n <= 5.

**2. What we believe the correct statement is**

The proof idea of the paper (mapping partitions to chord diagrams and borrowing the Touchard-Riordan count) is not salvageable as is, because the map is injective but not surjective: chord diagrams whose chords would force forbidden adjacencies are never hit. The fibers are exactly the descent runs.

The statement that survives, and that we verified to n = 7 in all three systems, is:

Q_n(Set(n, ell)) = Sum over k of A(n-k, ell, q) * s_(n-k, 1^k)

where A(m, ell, q) is the crossing polynomial of adjacency-free partitions (partitions whose descent set is empty). The collapse lemma behind it, c_D(n, ell) = A(n - |D|, ell, q) for every descent set D, is a new structural identity and, we think, of independent interest. The corrected theorem is manifestly Schur-positive: every coefficient is a sum of monomials with nonnegative integer coefficients, visible by inspection.

Your Corollary 3.3, parts (3) and (4), is correct and quite beautiful: the noncrossing coefficients are exactly the Motzkin triangle numbers (noncrossing matchings with unmatched points allowed), the row sums are the Motzkin numbers, and the hook totals recover the Narayana numbers.

**3. New identities the audit surfaced**

Reading the corrected table at q = 1 gives a new identity, verified in all three systems: a(m, ell) = Sum_j (-1)^j C(m-1, j) S(m-j, ell), and Sum over ell of a(m, ell) = Bell(m-1): a bijection between adjacency-free partitions of m points and set partitions of m-1.

**4. Everything is public and reproducible**

- Full workflow repository (Wolfram package, teaching notebook, audit scripts, Python verifier, Lean 4 library, figures, documentation): https://github.com/zuwasi/mantiq-touchard-riordan
- A technical slide deck describing the audit: https://zuwasi.github.io/Public-html-pages/mantiq-touchard-riordan/technical/
- A short demo video (about 5 minutes) walking through the audit; I am happy to send it or walk you through it on a call.

Everything reruns in minutes on a laptop: the Wolfram build takes about 19 seconds, the Python verifier a few seconds, the Lean build under two seconds.

**5. Why I am writing to you**

Beyond this specific paper: we believe this workflow (Mathematica as the engine of exact experimentation, with Lean 4 as the formal guard) is genuinely useful for applied combinatorics research and teaching. It found the gap, produced the corrected statement, and now protects it through every future revision. If you find this interesting, I would be glad to set up a call, run the workflow live with you on a question of your choosing, or help you apply it to your current work. The teaching notebook we built from your paper's definitions could also serve directly as course material on Schur positivity and set partitions.

With respect for your work, and with thanks for the stimulating paper,

Daniel Liezrowice
Engineering Software Lab (ESL)
www.eswlab.com | sales@eswlab.com
