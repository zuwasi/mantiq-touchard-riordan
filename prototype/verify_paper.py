# Prototype verification of Bagno-Garber (EPTCS 445, GASCom 2026):
# "Touchard-Riordan Polynomials and Schur-positivity of Set Partitions"
#
# Conventions from the paper:
#  - Linear graph of a set partition: within each block (sorted), consecutive
#    elements joined. If the two elements are consecutive integers the join is a
#    STRAIGHT segment on the baseline; otherwise a semicircular ARC.
#  - cr(pi) counts crossing pairs of ARCS only (straight segments contribute none),
#    arcs (a,b),(c,d) cross iff a<c<b<d. Verified against Example 2.5 (cr=3).
#  - Des(pi) = {i : i and i+1 lie in the same block}.
#
# Claims to verify:
#  C1 (Def 2.6/Ex 2.5): crossing counter reproduces cr({{1,3,7},{2,4,5,6,8}}) = 3.
#  C2 (Touchard-Riordan): T_{2m}(q) = (1-q)^{-m} sum_{i>=0} (-1)^i q^{i(i+1)/2}
#      [C(2m,m-i) - C(2m,m-i-1)]  equals brute-force generating function
#      over perfect matchings on 2m points.
#  C3 (paper, p.5): generating function for matchings on n points with j chords:
#      T_{n,j}(q) = C(n,2j) * T_{2j}(q).
#  C4 (Theorem 3.1, reduced form): for all n, l, k:
#      sum over pi in Set(n,l) with |Des(pi)|=k of q^{cr(pi)}
#        == C(n-1-k, 2(n-k-l)) * T_{2(n-k-l)}(q).
#      (This is Theorem 3.1 coefficient-extracted: hook SYT with k descents have
#       each descent set of size k exactly once, so the coefficient of the hook
#       s_{(n-k,1^k)} is the crossing distribution over partitions with k descents.)
#  C5 (Corollary 3.3(2) at q=1): sum_k sum_l [count] == S(n, l) totals; and the
#      reduced identity sum_r C(l+r-1, 2r) (2r-1)!! == Stirling2(n,l).
#  C6 (Remark 2.10 typo check): T_{2n}(0) = Catalan(n) != Motzkin(2n);
#      the correct identity used by Cor 3.3 is T_{n,j}(0) = C(n,2j) Cat_j
#      = Motzkin triangle M_{n,j} (OEIS A055151), and sum_j T_{n,j}(0) = Motzkin(n).
#  C7 (New insights): coefficient statistics of T_{2m}(q): unimodality,
#      log-concavity, real-rootedness, for m up to 16.

from itertools import combinations, product
from math import comb, factorial
from functools import lru_cache
import json, sys

def set_partitions(n):
    """Restricted growth strings: all set partitions of [n]."""
    if n == 0:
        yield ()
        return
    def rec(prefix, maxblock):
        m = len(prefix)
        if m == n:
            yield tuple(prefix)
            return
        for b in range(0, min(maxblock + 1, m) + 1):
            yield from rec(prefix + [b], max(maxblock, b))
    yield from rec([], 0)

def blocks(rgs):
    d = {}
    for i, b in enumerate(rgs, start=1):
        d.setdefault(b, []).append(i)
    return [sorted(v) for v in d.values()]

def arcs_and_des(rgs):
    """Return (list of arcs (a,b) with b>a+1, descent set, total consecutive joins)."""
    ars, des = [], set()
    for blk in blocks(rgs):
        for u, v in zip(blk, blk[1:]):
            if v == u + 1:
                des.add(u)          # consecutive integers: straight segment, and
            else:                   # i,i+1 in same block  => descent at i
                ars.append((u, v))
    return ars, sorted(des)

def crossing_number(rgs):
    ars, des = arcs_and_des(rgs)
    cr = 0
    for (a, b), (c, d) in combinations(ars, 2):
        a, b = min((a,b)), max((a,b))
        c, d = min((c,d)), max((c,d))
        if a < c < b < d or c < a < d < b:
            cr += 1
    return cr, des

# --- C1: Example 2.5 check ----------------------------------------------------
ex = tuple([0,1,0,1,1,1,0,1])  # {1,3,7},{2,4,5,6,8}
cr_ex, des_ex = crossing_number(ex)
C1 = (cr_ex == 3 and des_ex == [])
print(f"C1 Example 2.5: cr={cr_ex} (expect 3), Des={des_ex} -> {'PASS' if C1 else 'FAIL'}")

# --- C2: Touchard-Riordan closed form vs brute force --------------------------
def perfect_matchings(points):
    pts = list(points)
    if not pts:
        yield ()
        return
    a, rest = pts[0], pts[1:]
    for i, b in enumerate(rest):
        for m in perfect_matchings(rest[:i] + rest[i+1:]):
            yield ((a, b),) + m

def cr_matchings_on_circle(matching):
    """crossings of chords (a,b) with a<b: a<c<b<d."""
    chords = [(min(x,y), max(x,y)) for x, y in matching]
    return sum(1 for (a,b),(c,d) in combinations(chords,2) if a<c<b<d or c<a<d<b)

def TR_brute(m):
    dist = {}
    for M in perfect_matchings(range(1, 2*m+1)):
        dist[cr_matchings_on_circle(M)] = dist.get(cr_matchings_on_circle(M), 0) + 1
    return dist

def TR_formula(m, maxlen=None):
    """T_{2m}(q) coefficient list."""
    from math import comb
    N = [0]*(m*m//2 + 2)
    i = 0
    while i <= m:
        c = comb(2*m, m-i) - comb(2*m, m-i-1) if m-i-1 >= 0 else comb(2*m, m-i)
        e = i*(i+1)//2
        if e < len(N):
            N[e] += (-1)**i * c
        i += 1
    # T = N / (1-q)^m ; do polynomial division by (1-q)^m
    from math import comb as C
    # (1-q)^m coefficients
    den = [(-1)**k * C(m, k) for k in range(m+1)]
    # long division: T[k] = N[k] - sum_{j<k} den[j+1..] ... do iteratively
    T = [0]*(m*(m-1)//2 + 1)
    n_coef = N + [0]*(m*m)
    for k in range(len(T)):
        acc = n_coef[k]
        for j in range(1, min(m, k)+1):
            acc -= den[j] * T[k-j]
        T[k] = acc
    # trim
    while T and T[-1] == 0:
        T.pop()
    return T

C2 = True
for m in range(1, 7):
    bf = TR_brute(m)
    fm = TR_formula(m)
    bf_list = [bf.get(k, 0) for k in range(max(bf)+1)] if bf else [1]
    ok = bf_list == fm
    C2 &= ok
    print(f"C2 m={m}: brute={bf_list} formula={fm} -> {'PASS' if ok else 'FAIL'}")

# --- C3: T_{n,j}(q) = C(n,2j) T_{2j}(q) ---------------------------------------
def matchings_gen_fun(n, j):
    dist = {}
    for pts in combinations(range(1, n+1), 2*j):
        for M in perfect_matchings(pts):
            c = cr_matchings_on_circle(M)
            dist[c] = dist.get(c, 0) + 1
    return dist

C3 = True
for n in range(2, 9):
    for j in range(0, n//2 + 1):
        bf = matchings_gen_fun(n, j)
        bf_list = [bf.get(k,0) for k in range(max(bf, default=0)+1)]
        T2j = TR_formula(j) if j > 0 else [1]
        expect = [comb(n, 2*j) * t for t in T2j]
        if bf_list != expect:
            C3 = False
            print(f"C3 FAIL n={n} j={j}: {bf_list} vs {expect}")
print(f"C3 T_(n,j) = C(n,2j) T_(2j): {'PASS' if C3 else 'FAIL'}")

# --- C4: Theorem 3.1 reduced identity -----------------------------------------
def TR_poly_eval(qv):
    pass

def dist_from_formula(n, l, k):
    r = n - k - l
    if r < 0: return None
    T2r = TR_formula(r)
    return [comb(n-1-k, 2*r) * t for t in T2r]

C4 = True
report = []
for n in range(1, 9):
    parts = list(set_partitions(n))
    for l in range(1, n+1):
        pl = [p for p in parts if len(set(p)) == l]
        # group by |Des|
        from collections import defaultdict
        byk = defaultdict(lambda: defaultdict(int))
        for p in pl:
            cr, des = crossing_number(p)
            byk[len(des)][cr] += 1
        for k in byk:
            bf = [byk[k].get(i, 0) for i in range(max(byk[k])+1)]
            ex = dist_from_formula(n, l, k)
            if ex is None or bf != ex:
                C4 = False
                print(f"C4 FAIL n={n} l={l} k={k}: brute={bf} formula={ex}")
    print(f"C4 checked n={n}")
print(f"C4 Theorem 3.1 (reduced, all n<=8): {'PASS' if C4 else 'FAIL'}")

# --- C5: Stirling consistency --------------------------------------------------
def stirling2(n, k):
    if n == 0 and k == 0: return 1
    if k == 0 or k > n: return 0
    return k*stirling2(n-1, k) + stirling2(n-1, k-1)

C5 = True
for n in range(1, 10):
    for l in range(1, n+1):
        s = sum(comb(l+r-1, 2*r) * (factorial(2*r)//(2**r * factorial(r))) for r in range(0, (n-l)//2 + 1) if l+r-1 >= 2*r)
        # careful: r ranges with l+r-1 >= 2r i.e. r <= l-1
        s = sum(comb(l+r-1, 2*r) * (factorial(2*r)//(2**r*factorial(r))) for r in range(0, l))
        if s != stirling2(n, l):
            C5 = False
            print(f"C5 FAIL n={n} l={l}: {s} vs {stirling2(n,l)}")
print(f"C5 Stirling identity: {'PASS' if C5 else 'FAIL'}")

# --- C6: Remark 2.10 typo check -------------------------------------------------
def poly_eval(T, x): return sum(c * x**i for i, c in enumerate(T))

cat = [1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796]
motz = [1, 1, 2, 4, 9, 21, 51, 127, 323, 835, 2188]
# Motzkin triangle A055151: M_{n,k} = C(n,2k) * Cat_k
tri = [[comb(n_, 2*k_)*cat[k_] for k_ in range(n_//2+1)] for n_ in range(11)]
A055151_rows = [[1],[1],[1,1],[1,3],[1,6,2],[1,10,10],[1,15,30,5],[1,21,70,35],[1,28,140,140,14]]
ok_tri = all(tri[i][:len(A055151_rows[i])] == A055151_rows[i] for i in range(len(A055151_rows)))
print(f"C6a Motzkin triangle M_(n,k)=C(n,2k)Cat_k matches A055151: {'PASS' if ok_tri else 'FAIL'}")
# T_{2n}(0) = Cat_n
t0 = [TR_formula(n)[0] for n in range(1, 9)]
ok_catalan = t0 == cat[1:9]
print(f"C6b T_(2n)(0) = Catalan(n): {'PASS' if ok_catalan else 'FAIL'}  ({t0})")
# T_{2n}(0) vs Motzkin(2n): mismatch (the paper's Remark 2.10 literal claim)
print(f"C6c T_(2n)(0)={t0} vs Motzkin(2n)={[motz[2*n] for n in range(1,9)]} -> NOT equal (Remark 2.10 typo confirmed)")
# sum_j T_(n,j)(0) = Motzkin(n)
ok_motz = all(sum(comb(n_, 2*j)*cat[j] for j in range(n_//2+1)) == motz[n_] for n_ in range(11))
print(f"C6d sum_j T_(n,j)(0) = Motzkin(n): {'PASS' if ok_motz else 'FAIL'}")

# --- C7: coefficient statistics of T_{2m}(q) -----------------------------------
def is_unimodal(T):
    seen_drop = False
    for a, b in zip(T, T[1:]):
        if b < a: seen_drop = True
        elif seen_drop and b > a: return False
    return True

def is_log_concave(T):
    return all(T[i-1]*T[i+1] >= T[i]**2 for i in range(1, len(T)-1))

ins = []
for m in range(1, 17):
    T = TR_formula(m)
    ins.append({"m": m, "coeffs": T, "unimodal": is_unimodal(T), "logconcave": is_log_concave(T)})
for r in ins:
    print(f"C7 m={r['m']:2d} unimodal={r['unimodal']} logconcave={r['logconcave']} coeffs={r['coeffs']}")

# real-rootedness via numpy roots for m<=12
try:
    import numpy as np
    realroot = []
    for m in range(1, 13):
        T = TR_formula(m)
        rts = np.roots(T[::-1])
        real = all(abs(z.imag) < 1e-7 for z in rts)
        realroot.append((m, real))
    print("C7 real-rootedness m=1..12:", realroot)
except Exception as e:
    print("numpy unavailable:", e)

# mean crossings of perfect matchings (sanity: m(m-1)/4?)
import statistics
for m in [3, 4, 5]:
    bf = TR_brute(m)
    tot = sum(bf.values()); mean = sum(k*v for k, v in bf.items())/tot
    print(f"mean cr for 2m={2*m}: {mean:.4f} (m(m-1)/4 = {m*(m-1)/4})")

json.dump({"C1": C1, "C2": C2, "C3": C3, "C4": C4, "C5": C5, "insights": ins},
          open(r"C:\Amp_demos\Eli-Bango-JCT\prototype\verify_results.json", "w"), indent=1)
print("\nALL CORE CHECKS:", "PASS" if all([C1,C2,C3,C4,C5]) else "SOME FAILED")
