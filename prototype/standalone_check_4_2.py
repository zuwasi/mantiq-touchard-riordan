"""Standalone re-derivation of Q_4(Set(4,2)) with no dependencies.

Definitions follow the paper (arXiv 2606.13149):
  Set(4,2)   : set partitions of {1,2,3,4} into 2 blocks.
  Des(pi)    : {i : i and i+1 lie in the same block}            (Section 3 usage)
  cr(pi)     : crossings of the linear graph: arcs join the elements of a block
               that are separated by other blocks (consecutive block elements at
               distance >= 2); adjacent integers joined by straight segments
               contribute no crossings. Arcs (a,b),(c,d) cross iff a<c<b<d.
  F_{n,D}(x) : sum over weakly increasing i1 <= ... <= i4 in [4] with a strict
               rise at every j in D, of x_{i1}...x_{i4}.
  Q_n(S)     : sum over pi in S of q^cr(pi) F_{n,Des(pi)}.

The paper's Theorem 3.1 RHS for (n,ell) = (4,2) is
  T_{3,2}(q) s_(4) + T_{2,1}(q) s_(3,1) + T_{1,0}(q) s_(2,1,1)
with T_{a,b}(q) = C(a,2b) T_{2b}(q), so T_{3,2} = C(3,4) T_4 = 0,
T_{2,1} = C(2,2) T_2 = 1, T_{1,0} = C(1,0) T_0 = 1:
  paper RHS = s_(3,1) + s_(2,1,1), i.e. coefficient 0 for s_(4).
"""
from itertools import permutations
from collections import defaultdict
from fractions import Fraction

# The 7 set partitions of {1,2,3,4} into 2 blocks.
parts = []
for s in range(1, 5):                              # singleton {s}
    others = frozenset(i for i in range(1, 5) if i != s)
    parts.append(frozenset((frozenset({s}), others)))
for a in range(2, 5):                              # block {1,a}, complement
    others = frozenset(i for i in range(1, 5) if i not in (1, a))
    parts.append(frozenset((frozenset({1, a}), others)))
parts = sorted(set(parts))
assert len(parts) == 7, len(parts)

def des(p):
    d = set()
    for b in p:
        s = sorted(b)
        d |= {u for u, v in zip(s, s[1:]) if v == u + 1}
    return d

def arcs(p):
    out = []
    for b in p:
        s = sorted(b)
        out += [(u, v) for u, v in zip(s, s[1:]) if v - u >= 2]
    return out

def cr(p):
    a = sorted(arcs(p))
    return sum(1 for i in range(len(a)) for j in range(i + 1, len(a))
               if a[i][0] < a[j][0] < a[i][1] < a[j][1])

print("The 7 partitions of Set(4,2) with Des(pi) and cr(pi):")
for p in parts:
    print("  ", sorted(sorted(b) for b in p), "Des =", sorted(des(p)), "cr =", cr(p))

# Q_4(Set(4,2)): one term per weakly increasing sequence, i.e. per monomial.
poly = defaultdict(lambda: defaultdict(int))
for p in parts:
    D, c = des(p), cr(p)
    def gen(prefix, remaining, first):
        if remaining == 0:
            yield list(prefix)
            return
        for v in range(first, 5):
            prefix.append(v)
            yield from gen(prefix, remaining - 1, v)
            prefix.pop()
    for seq in gen([], 4, 1):
        alpha = tuple(seq.count(v) for v in range(1, 5))
        rises = {j for j in range(1, 4) if seq[j - 1] < seq[j]}
        if D <= rises:
            poly[alpha][c] += 1

print("Q_4(Set(4,2)) monomial coefficients (alpha = exponents of x1..x4):")
for e in sorted(poly, key=lambda t: tuple(-x for x in t)):
    print("  ", e, "->", dict(poly[e]))

# Symmetry check: the coefficient of x^alpha depends only on the multiset.
for alpha in poly:
    key = tuple(sorted(alpha, reverse=True))
    orbit = set(permutations(key))
    for qe in poly[alpha]:
        vals = {poly.get(a, {}).get(qe, 0) for a in orbit}
        assert len(vals) == 1, (alpha, qe, vals)
print("Symmetry check passed (coefficient constant on variable permutations).")

# Schur expansion via inverse Kostka (exact fractions).
def ssyt_count(shape, weight):
    cells = [(r, c) for r, rowlen in enumerate(shape) for c in range(rowlen)]
    fill = {}
    def rec(i, rem):
        if i == len(cells):
            return 1
        r, c = cells[i]
        t = 0
        for v in range(1, 5):
            if rem[v - 1] == 0: continue
            if r > 0 and v <= fill[(r - 1, c)]: continue
            if c > 0 and v < fill[(r, c - 1)]: continue
            fill[(r, c)] = v; rem[v - 1] -= 1
            t += rec(i + 1, rem)
            rem[v - 1] += 1; del fill[(r, c)]
        return t
    return rec(0, list(weight) + [0] * (4 - len(weight)))

shapes = [(4,), (3, 1), (2, 2), (2, 1, 1), (1, 1, 1, 1)]
# Distinct weight partitions; the coefficient is constant on each orbit.
weights = sorted({tuple(sorted(alpha, reverse=True)) for alpha in poly},
                 key=lambda t: tuple(-x for x in t))
K = [[ssyt_count(l, list(mu)) for mu in weights] for l in shapes]
print("Kostka matrix K[lambda][mu] (rows = Schur shapes, cols = monomial weights):")
for i, l in enumerate(shapes):
    print("  ", l, K[i])

n = len(shapes)
nw = len(weights)
for qe in sorted({q for e in poly for q in poly[e]}):
    # s_lambda = sum_mu K[lamba][mu] m_mu  =>  f_mu = sum_lambda c_lambda K[lambda][mu],
    # i.e. K^T c = f; solve that system exactly.
    f = [poly[weights[j]].get(qe, 0) for j in range(nw)]
    a = [[K[i][j] for i in range(n)] + [f[j]] for j in range(nw)]
    for col in range(n):
        piv = next(r for r in range(col, nw) if a[r][col] != 0)
        a[col], a[piv] = a[piv], a[col]
        pv = a[col][col]
        a[col] = [Fraction(x, pv) for x in a[col]]
        for r in range(nw):
            if r != col and a[r][col]:
                fac = a[r][col]
                a[r] = [x - fac * y for x, y in zip(a[r], a[col])]
    coeffs = [int(a[i][nw]) for i in range(n)]
    print(f"Schur coefficients at q^{qe}:",
          {shapes[i]: coeffs[i] for i in range(n) if coeffs[i]})
    if qe == 0:
        q0 = {shapes[i]: coeffs[i] for i in range(n) if coeffs[i]}
    else:
        q1 = {shapes[i]: coeffs[i] for i in range(n) if coeffs[i]}

# Verdict: the true expansion is q*s_(4) + s_(3,1) + s_(2,1,1); the paper's
# Theorem 3.1 RHS is s_(3,1) + s_(2,1,1) and misses the q*s_(4) term.
assert q0 == {(3, 1): 1, (2, 1, 1): 1}, q0
assert q1 == {(4,): 1}, q1
print()
print("VERIFIED: Q_4(Set(4,2)) = q*s_(4) + s_(3,1) + s_(2,1,1).")
print("Paper Theorem 3.1 RHS for (4,2) = s_(3,1) + s_(2,1,1) (T_{3,2} = 0),")
print("so the theorem is FALSE: it misses the q*s_(4) term.")
print("Equivalently, LHS - paper RHS = q*h_4 = q*s_(4).")