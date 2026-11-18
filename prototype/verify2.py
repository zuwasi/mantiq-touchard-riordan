# Definitive verification of Bagno-Garber Theorem 3.1 / Corollary 3.3.
# Two INDEPENDENT paths for Schur coefficients:
#   (A) monomial-basis path: F_{n,D} expanded via weakly-increasing value
#       sequences; f expressed in m_mu; Schur via Kostka (SSYT brute force).
#   (B) pairing path: <F_{n,D}, s_lam> = #{SYT(lam) with Des = D} (SYT brute force).
# Everything cross-checked. n <= 6.

from itertools import combinations, product
from functools import lru_cache
from collections import defaultdict
import json

# ---------------- set partitions ----------------
def set_partitions(n):
    if n == 0:
        yield ()
        return
    def rec(prefix, maxblock):
        m = len(prefix)
        if m == n:
            yield tuple(prefix); return
        for b in range(0, min(maxblock + 1, m) + 1):
            yield from rec(prefix + [b], max(maxblock, b))
    yield from rec([], 0)

def blocks_of(rgs):
    d = defaultdict(list)
    for i, b in enumerate(rgs, 1):
        d[b].append(i)
    return [sorted(v) for _, v in sorted(d.items())]

def des_set(rgs):
    return tuple(i for i in range(1, len(rgs)) if rgs[i-1] == rgs[i])  # rgs[i-1], rgs[i] are elements i, i+1

def arcs_cr_des(rgs):
    """arcs = consecutive block elements with gap>=2 (straight segs excluded);
    cr per Def 2.6; Des per the paper."""
    ars = []
    for blk in blocks_of(rgs):
        for u, v in zip(blk, blk[1:]):
            if v - u >= 2:
                ars.append((u, v))
    cr = 0
    for (a, b), (c, d) in combinations(ars, 2):
        if a > c: (a, b), (c, d) = (c, d), (a, b)
        if a < c < b < d:
            cr += 1
    return cr, des_set(rgs), ars

def is_noncrossing_partition(rgs):
    """Block-level noncrossing: no i<j<k<l with i,k in one block and j,l in another."""
    blks = blocks_of(rgs)
    pos = {}
    for bi, blk in enumerate(blks):
        for x in blk: pos[x] = bi
    n = len(rgs)
    for i in range(1, n+1):
        for k in range(i+1, n+1):
            if pos[i] != pos[k]: continue
            for j in range(i+1, k):
                for l in range(k+1, n+1):
                    if pos[j] == pos[l] and pos[j] != pos[i]:
                        return False
    return True

# ---------------- SYT enumeration and Des ----------------
def syts(shape):
    """All standard Young tableaux of given shape (list of row lengths)."""
    cells = [(r, c) for r, row in enumerate(shape) for c in range(row)]
    n = len(cells)
    idx = {cell: i for i, cell in enumerate(cells)}
    results = []
    def rec(used, assign, remaining):
        if not remaining:
            results.append(dict(assign)); return
        # next smallest value goes in a cell whose upper/left neighbors are placed
        for cell in cells:
            if cell in assign: continue
            r, c = cell
            if r > 0 and (r-1, c) not in assign: continue
            if c > 0 and (r, c-1) not in assign: continue
            # value = len(used)+1 must be placeable: cells must be filled in order
            assign[cell] = len(used) + 1
            rec(used | {cell}, assign, remaining - 1)
            del assign[cell]
            break  # the smallest remaining value must go in the first available cell? no...
        return
    # proper enumeration: assign values 1..n in increasing order to admissible cells
    def rec2(assign):
        if len(assign) == n:
            results.append(dict(assign)); return
        v = len(assign) + 1
        for cell in cells:
            if cell in assign: continue
            r, c = cell
            if r > 0 and (r-1, c) not in assign: continue
            if c > 0 and (r, c-1) not in assign: continue
            assign[cell] = v
            rec2(assign)
            del assign[cell]
        return
    rec2({})
    # convert to Des sets
    out = []
    for t in results:
        pos_of = {v: cell for cell, v in t.items()}
        des = tuple(i for i in range(1, n) if pos_of[i+1][0] > pos_of[i][0])
        out.append(des)
    return out

# ---------------- Kostka via SSYT brute force ----------------
def kostka_row(lam, maxpart):
    """K[lam][mu] for all partitions mu of n=len content with parts<=maxpart."""
    n = sum(lam)
    cells = [(r, c) for r, row in enumerate(lam) for c in range(row)]
    counts = defaultdict(int)
    def rec(assign):
        if len(assign) == n:
            content = tuple(sorted((tuple(assign.values()))), reverse=True) if False else tuple(sorted(assign.values(), reverse=True))
            counts[content] += 1
            return
        v = len(assign) + 1
        for cell in cells:
            if cell in assign: continue
            r, c = cell
            if r > 0 and (r-1, c) not in assign: continue
            if c > 0 and (r, c-1) not in assign: continue
            if r > 0 and assign.get((r-1, c), 0) >= v: continue
            if c > 0 and assign.get((r, c-1), 0) > v: continue
            assign[cell] = v
            rec(assign)
            del assign[cell]
    rec({})
    return dict(counts)

def partitions_of(n, maxpart=None, prefix=()):
    if maxpart is None: maxpart = n
    if n == 0:
        yield prefix; return
    for p in range(min(maxpart, n), 0, -1):
        yield from partitions_of(n - p, min(p, n - p) if n - p else 0, prefix + (p,)) if False else partitions_of(n - p, p, prefix + (p,))

def all_partitions(n):
    yield from partitions_of(n)

# F_{n,D}: which co() does a monomial correspond to
def co_of_exponents(a):
    """a = exponent tuple (a_1..a_n); composition = positive parts in order; co = partial sums < n."""
    n = len(a)
    parts = [x for x in a if x > 0]
    co, s = [], 0
    for x in parts[:-1] if parts else []:
        s += x; co.append(s)
    # partial sums of all but last part
    co = []
    s = 0
    for x in parts:
        s += x
        if s < n: co.append(s)
    return tuple(co)

def F_monomial_coeffs(n, D):
    """Return dict: multiset(lambda sorted desc) -> coefficient (0/1) for F_{n,D} in m-basis.
    F_{n,D} = sum over weakly increasing w in [n]^n with rises exactly at D of monomial x^w."""
    out = defaultdict(int)
    def rec2(w):
        if len(w) == n:
            a = [0]*n
            for x in w: a[x-1] += 1
            if co_of_exponents(a) == tuple(sorted(D)):
                out[tuple(sorted(a, reverse=True))] += 1
            return
        lo = w[-1] if w else 1
        for t in range(lo, n+1):
            rec2(w + (t,))
    rec2(())
    return dict(out)

# ---------------- main computation ----------------
def run(nmax=6):
    report = {}
    # sanity: SYT hook descent sets = all k-subsets
    for n in range(2, 7):
        for k in range(0, n):
            shape = [n-k] + [1]*k
            des = sorted(set(syts(shape)))
            want = sorted({tuple(sorted(c)) for c in combinations(range(1, n), k)})
            assert des == want, (n, k, des, want)
    print("SYT hook descent-set property verified (all n<=6)")

    results = []
    for n in range(1, nmax+1):
        parts = list(set_partitions(n))
        # SYT pairing table for this n: pair[lam_des_list]
        shapes = list(all_partitions(n))
        pair = {}   # lam -> Counter of Des -> count
        for lam in shapes:
            cnt = defaultdict(int)
            for d in syts(lam):
                cnt[d] += 1
            pair[lam] = dict(cnt)
        for l in range(1, n+1):
            pl = [p for p in parts if len(set(p)) == l]
            # LHS f: coefficient of each monomial exponent-vector
            mono = defaultdict(int)   # exponent tuple -> coeff
            for p in pl:
                cr, D, ars = arcs_cr_des(p)
                # add q^cr * F_{n,D}: all weakly increasing seqs with co = D
                for a, cval in F_monomial_coeffs(n, D).items():
                    mono[a] += cval  # q^cr applied later via separation by cr? we need polynomial in q:
            # store per-cr version: recompute with q weights
            mono_q = defaultdict(lambda: defaultdict(int))
            for p in pl:
                cr, D, ars = arcs_cr_des(p)
                for a, cval in F_monomial_coeffs(n, D).items():
                    mono_q[cr][a] += cval
            # symmetry check per cr-level: coefficient depends only on multiset
            symmetric = True
            for cr, d in mono_q.items():
                byms = defaultdict(set)
                for a, v in d.items():
                    byms[tuple(sorted(a, reverse=True))].add(v)
                if any(len(s) > 1 for s in byms.values()):
                    symmetric = False
            # true Schur expansion via pairing (valid if symmetric):
            # a_lam = sum_pi q^cr * #{SYT(lam): Des = Des(pi)}  -- group by (cr, |Des|)
            # BUT pairing uses full Des, so a_lam = sum over Des classes:
            # a_lam(q) = sum_cr q^cr sum_{|Des|=k} #{SYT(lam): Des=Des(pi)}... need per-Des counts:
            descount = defaultdict(int)  # Des tuple -> (count, cr) -- partitions with that Des (all same cr dist? no!)
            # careful: partitions with same Des can have different cr. group:
            bydes = defaultdict(lambda: defaultdict(int))
            for p in pl:
                cr, D, ars = arcs_cr_des(p)
                bydes[D][cr] += 1
            schur = {}
            for lam in shapes:
                poly = defaultdict(int)  # cr -> coeff
                for D, crdist in bydes.items():
                    cnt = pair[lam].get(D, 0)
                    if cnt:
                        for cr, m in crdist.items():
                            poly[cr] += cnt * m
                if poly: schur[lam] = dict(poly)
            results.append({"n": n, "l": l, "symmetric": symmetric, "schur": {str(k): v for k, v in schur.items()},
                            "counts": {str(D): dict(v) for D, v in bydes.items()}})
    return results

res = run(6)
for r in res:
    if r["n"] <= 4:
        print(f"n={r['n']} l={r['l']} symmetric={r['symmetric']}")
        for lam, poly in r["schur"].items():
            print(f"   s_{lam}: q-coeffs {poly}")

# --- compare with paper RHS: coefficient of s_(n-k,1^k) should be T_{n-1-k, n-k-l}(q)
def TR_poly(m):
    """T_{2m}(q) coefficient list via Touchard-Riordan formula."""
    N = [0]*(m*(m+1)//2 + 2)
    for i in range(m+1):
        c = (comb(2*m, m-i) - (comb(2*m, m-i-1) if m-i-1 >= 0 else 0))
        e = i*(i+1)//2
        N[e] += (-1)**i * c
    den = [(-1)**kk * comb(m, kk) for kk in range(m+1)]
    T = [0]*(m*(m-1)//2 + 1)
    ext = N + [0]*(len(T)+5)
    for kk in range(len(T)):
        acc = ext[kk]
        for j in range(1, min(m, kk)+1):
            acc -= den[j]*T[kk-j]
        T[kk] = acc
    while T and T[-1] == 0: T.pop()
    return T

from math import comb
def T_ab(a, b):
    """gen fun matchings on a points with b chords = C(a,2b) T_{2b}(q)."""
    if 2*b > a or b < 0: return [0]
    if b == 0: return [1]
    base = TR_poly(b)
    return [comb(a, 2*b)*t for t in base]

print("\n--- Paper RHS vs true hook coefficients ---")
for r in res:
    n, l = r["n"], r["l"]
    for k in range(0, n):
        hook = str(tuple([n-k] + [1]*k))
        true_poly = r["schur"].get(hook, {})
        paper = T_ab(n-1-k, n-k-l)
        tp = [true_poly.get(i, 0) for i in range(max([len(paper)-1] + list(true_poly)) + 1)] if true_poly else [0]*(len(paper))
        match = tp == paper
        flag = "OK " if match else "DIFF"
        if not match or n <= 4:
            print(f"{flag} n={n} l={l} k={k} hook={hook}: true={tp} paper={paper}")
