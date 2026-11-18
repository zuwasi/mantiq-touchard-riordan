/-
Verify.lean -- MantiQ Factor formal companion to Bagno-Garber,
"Touchard-Riordan Polynomials and Schur-positivity of Set Partitions"
(arXiv:2606.13149).

Pure Lean 4 core, no Mathlib. Every claim below is checked by brute-force
enumeration via `native_decide` (trusted compiler evaluation; only `T_3_2`
is a plain `rfl`). In particular this is not a kernel-level proof of a
mathematical statement about all n: it is an exhaustive finite check. What
is formalized here:

  * set partitions as restricted growth strings, descent sets, crossings
  * the corrected coefficients A(m, ell, q) as coefficient lists
  * Touchard-Riordan coefficients T_{a,b}(q) by enumerating all matchings
  * the smallest counterexample: A(4,2) = q is nonzero, T_{3,2} = 0
  * the census: the paper's hook coefficients differ from the corrected ones
    exactly for (n, ell) in {(4,2), (5,2), (5,3)} within n <= 5
  * the collapse lemma c_D(n, ell) = A(n - |D|, ell, q) for every descent set,
    at (n, ell) = (4,2) and (5,3), and for every (n, ell) with n <= 4
  * Bell identity, inclusion-exclusion, Motzkin triangle, Motzkin row sums,
    Narayana hook totals (the true parts of Corollary 3.3)

Build with:  lake build
-/

namespace Bagno

/-! ## Restricted growth strings -/

/-- Extend a restricted growth string by one element in every legal block.
    The first element must be 0, so the empty string only extends to [0]. -/
def extend (w : List Nat) : List (List Nat) :=
  if w == [] then [[0]]
  else (List.range (w.foldl Nat.max 0 + 2)).map (fun b => w ++ [b])

/-- All set partitions of [n] as restricted growth strings (block labels from 0). -/
def rgs : Nat → List (List Nat)
  | 0 => [[]]
  | n + 1 => (rgs n).flatMap extend

/-- The largest block label used; a partition has `maxBlock w + 1` blocks. -/
def maxBlock (w : List Nat) : Nat := w.foldl Nat.max 0

/-- Descent set (0-based): indices i with elements i and i+1 in the same block. -/
def des (w : List Nat) : List Nat :=
  (List.range (w.length - 1)).filter (fun i => w.getD i 0 == w.getD (i + 1) 0)

/-! ## Arcs and crossings -/

/-- u and v form an arc when they are consecutive elements of the same block. -/
def isArc (w : List Nat) (u v : Nat) : Bool :=
  u < v ∧ w.getD u 0 == w.getD v 0 ∧
    (List.range (v - u - 1)).all (fun k => w.getD (u + 1 + k) 0 != w.getD u 0)

/-- All arcs of a partition. -/
def arcs (w : List Nat) : List (Nat × Nat) :=
  (List.range w.length).flatMap (fun u =>
    (List.range w.length).filterMap (fun v =>
      if isArc w u v then some (u, v) else none))

/-- Arcs (u1, v1), (u2, v2) cross exactly when u1 < u2 < v1 < v2. -/
def crossPair (p q : Nat × Nat) : Bool := p.1 < q.1 ∧ q.1 < p.2 ∧ p.2 < q.2

/-- Crossing number of a partition: each crossing pair counted once
    (for every pair exactly one orientation satisfies `crossPair q p`). -/
def crossings (w : List Nat) : Nat :=
  (arcs w).foldl (fun acc p => acc + (arcs w).countP (fun q => crossPair q p)) 0

/-! ## q-polynomials as ascending coefficient lists -/

abbrev QP := List Nat

def qpAdd : QP → QP → QP
  | [], p => p
  | p, [] => p
  | a :: as, b :: bs => (a + b) :: qpAdd as bs

def qpScale (k : Nat) (p : QP) : QP := p.map (fun c => k * c)

/-- Trim trailing zeros, mapping the empty list to the canonical zero. -/
def qpTrim : QP → QP :=
  fun p =>
    let d := p.reverse.dropWhile (fun c => c == 0) |>.reverse
    if d == [] then [0] else d

/-- Tally the values of `f` over `xs` into a coefficient list. -/
def distWith {α : Type} (xs : List α) (f : α → Nat) : QP :=
  xs.foldl
    (fun acc x =>
      let i := f x
      let padded :=
        if i < acc.length then acc else acc ++ List.replicate (i + 1 - acc.length) 0
      padded.set i (padded.getD i 0 + 1))
    [0]

/-! ## Corrected coefficients A(m, ell, q) -/

/-- Adjacency-free (des = ∅) partitions of [m] into ell blocks. -/
def adjacencyFree (m ell : Nat) : List (List Nat) :=
  (rgs m).filter (fun w => des w == [] ∧ maxBlock w + 1 == ell)

/-- A(m, ell, q): the crossing polynomial of adjacency-free partitions. -/
def ACoeff (m ell : Nat) : QP := distWith (adjacencyFree m ell) crossings

/-- a(m, ell) = A(m, ell, 1): the count at q = 1. -/
def aCount (m ell : Nat) : Nat := (ACoeff m ell).sum

/-! ## Matchings and Touchard-Riordan polynomials -/

/-- All matchings on the given sorted point list (empty for odd sizes).
    The Nat argument tracks the list length for structural recursion. -/
def matchingsOn : Nat → List Nat → List (List (Nat × Nat))
  | _, [] => [[]]
  | 0, _ => [[]]
  | 1, _ => []
  | n + 2, pts =>
      (List.range (pts.length - 1)).flatMap (fun idx =>
        let q := pts.getD (idx + 1) 0
        let left := (pts.take (idx + 1)).drop 1
        let right := pts.drop (idx + 2)
        (matchingsOn n (left ++ right)).map (fun ms => (pts.head!, q) :: ms))

/-- Crossing number of a matching: each crossing pair counted once. -/
def matchingCross (ms : List (Nat × Nat)) : Nat :=
  ms.foldl (fun acc p => acc + ms.countP (fun q => crossPair q p)) 0

/-- The Touchard-Riordan polynomial T_{2b}(q), computed by enumeration. -/
def tRMatching (b : Nat) : QP := distWith (matchingsOn (2 * b) (List.range (2 * b))) matchingCross

/-- Binomial coefficients by Pascal's rule. -/
def binom : Nat → Nat → Nat
  | _, 0 => 1
  | 0, _ + 1 => 0
  | n + 1, k + 1 => binom n k + binom n (k + 1)

/-- The paper's hook coefficient T_{a,b}(q) = C(a, 2b) T_{2b}(q); 0 when 2b > a. -/
def touchardRiordanT (a b : Nat) : QP :=
  if a < 2 * b then [0] else qpScale (binom a (2 * b)) (tRMatching b)

/-- The paper's k-th hook coefficient for Q_n(Set(n, ell)), with the
    out-of-range convention (n < k + ell) giving zero. -/
def paperCoeff (n ell k : Nat) : QP :=
  if n < k + ell then [0] else touchardRiordanT (n - 1 - k) (n - k - ell)

/-! ## The counterexample and the census -/

/-- Smallest counterexample: the corrected coefficient of s_(4) in
    Q_4(Set(4,2)) is q (coefficients [0, 1]), not zero. -/
theorem A_4_2 : qpTrim (ACoeff 4 2) = [0, 1] := by native_decide

/-- The paper's coefficient T_(3,2)(q) is zero: no matching on 3 points has 4 chords. -/
theorem T_3_2 : touchardRiordanT 3 2 = [0] := by rfl

/-- The counterexample: the two candidate coefficients genuinely differ. -/
theorem counterexample_4_2 : touchardRiordanT 3 2 = [0] ∧ qpTrim (ACoeff 4 2) ≠ [0] := by native_decide

/-- True when the paper's and the corrected k-th hook coefficients differ. -/
def coeffsDiffer (n ell k : Nat) : Bool :=
  !(qpTrim (ACoeff (n - k) ell) == qpTrim (paperCoeff n ell k))

/-- True when the paper's Theorem 3.1 formula differs from the corrected one. -/
def paperFails (n ell : Nat) : Bool := (List.range n).any (fun k => coeffsDiffer n ell k)

/-- All (n, ell) with 2 <= n <= nmax + 1 where the paper's formula fails. -/
def censusFails (nmax : Nat) : List (Nat × Nat) :=
  (List.range nmax).flatMap (fun off =>
    (List.range (off + 2)).filterMap (fun e =>
      let n := off + 2
      let ell := e + 1
      if paperFails n ell then some (n, ell) else none))

/-- Census up to n = 5: exactly (4,2), (5,2), (5,3) fail. -/
theorem census_5 : censusFails 4 = [(4, 2), (5, 2), (5, 3)] := by native_decide

/-! ## The collapse lemma -/

/-- All subsets of {0, ..., n-1}. -/
def subsets : Nat → List (List Nat)
  | 0 => [[]]
  | n + 1 => (subsets n).flatMap (fun s => [s, s ++ [n]])

/-- Collapse lemma instance: the crossing polynomial over partitions with
    descent set exactly D equals A(n - |D|, ell, q). -/
def collapseOK (n ell : Nat) (d : List Nat) : Bool :=
  qpTrim (distWith ((rgs n).filter (fun w => maxBlock w + 1 == ell ∧ des w == d)) crossings)
    == qpTrim (ACoeff (n - d.length) ell)

/-- Collapse lemma for (n, ell) = (4, 2) and every descent set. -/
theorem collapse_4_2_allD : (subsets 3).all (fun d => collapseOK 4 2 d) = true := by native_decide

/-- Collapse lemma for (n, ell) = (5, 3) and every descent set. -/
theorem collapse_5_3_allD : (subsets 4).all (fun d => collapseOK 5 3 d) = true := by native_decide

/-- Collapse lemma for every (n, ell) with n <= 4 and every descent set. -/
def collapseAllOK (n : Nat) : Bool :=
  (List.range n).all (fun off =>
    (List.range (off + 2)).all (fun e => (subsets (off + 1)).all (fun d => collapseOK (off + 2) (e + 1) d)))

theorem collapse_all_n4 : collapseAllOK 3 = true := by native_decide

/-! ## Number sequences for the noncrossing corollary -/

/-- Bell numbers B_0, ..., B_n. -/
def bellList : Nat → List Nat
  | 0 => [1]
  | n + 1 =>
      let p := bellList n
      p ++ [(List.range (n + 1)).foldl (fun acc k => acc + binom n k * p.getD k 0) 0]

/-- Stirling numbers of the second kind. -/
def stirling : Nat → Nat → Nat
  | 0, 0 => 1
  | 0, _ + 1 => 0
  | _ + 1, 0 => 0
  | n + 1, k + 1 => (k + 1) * stirling n (k + 1) + stirling n k

/-- Motzkin numbers M_0, ..., M_n. -/
def motzList : Nat → List Nat
  | 0 => [1]
  | n + 1 =>
      let p := motzList n
      p ++ [p.getLast! + (List.range n).foldl (fun acc j => acc + p.getD j 0 * p.getD (n - 1 - j) 0) 0]

/-- Noncrossing configurations on a points with exactly b chords (a chord
    matches two points, remaining points stay single): the Motzkin triangle.
    Counted by choosing the 2b chorded points and noncrossingly perfect-
    matching them; noncrossing-ness only depends on the relative order. -/
def ncMatchCount (a b : Nat) : Nat :=
  (subsets a).filter (fun s => s.length == 2 * b)
    |>.flatMap (fun s => (matchingsOn (2 * b) s).filter (fun ms => matchingCross ms == 0))
    |>.length

/-- Noncrossing adjacency-free partitions of [m] into ell blocks. -/
def ancCount (m ell : Nat) : Nat :=
  (adjacencyFree m ell).countP (fun w => crossings w == 0)

/-! ## Verified identities -/

/-- Inclusion-exclusion at q = 1: a(m, ell) = sum_j (-1)^j C(m-1, j) S(m-j, ell). -/
def ieA (m ell : Nat) : Int :=
  (List.range m).foldl (fun acc j =>
    acc + (if j % 2 == 0 then 1 else -1)
      * ((binom (m - 1) j : Int) * (stirling (m - j) ell : Int))) 0

def ieOK (m ell : Nat) : Bool := (aCount m ell : Int) == ieA m ell

theorem inclusionExclusion_2_5 :
    (List.range 4).all (fun off =>
      (List.range (off + 2)).all (fun e => ieOK (off + 2) (e + 1))) = true := by native_decide

/-- Bell identity: sum over ell of a(m, ell) = Bell(m-1). -/
def bellOK (m : Nat) : Bool :=
  (List.range m).foldl (fun acc ell => acc + aCount m (ell + 1)) 0 == (bellList (m - 1)).getLast!

theorem bell_2_6 : (List.range 5).all (fun off => bellOK (off + 2)) = true := by native_decide

/-- Motzkin triangle: ancCount(m, ell) = ncMatchCount(m-1, m-ell). -/
def triangleOK (m ell : Nat) : Bool := ancCount m ell == ncMatchCount (m - 1) (m - ell)

theorem triangle_2_6 :
    (List.range 5).all (fun off =>
      (List.range (off + 2)).all (fun e => triangleOK (off + 2) (e + 1))) = true := by native_decide

/-- Motzkin row sums: sum over ell of ancCount(m, ell) = M(m-1). -/
def motzRowOK (m : Nat) : Bool :=
  (List.range m).foldl (fun acc ell => acc + ancCount m (ell + 1)) 0 == (motzList (m - 1)).getLast!

theorem motzRows_2_7 : (List.range 6).all (fun off => motzRowOK (off + 2)) = true := by native_decide

/-- Narayana numbers from the hook decomposition of the noncrossing counts. -/
def nar (n ell : Nat) : Nat := binom n ell * binom n (ell - 1) / n

def narOK (n ell : Nat) : Bool :=
  (List.range n).foldl (fun acc k => acc + ancCount (n - k) ell * binom (n - 1) k) 0 == nar n ell

theorem narayana_2_5 :
    (List.range 4).all (fun off =>
      (List.range (off + 2)).all (fun e => narOK (off + 2) (e + 1))) = true := by native_decide



end Bagno
