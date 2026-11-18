# Definitive verification of Bagno-Garber (arXiv 2606.13149) Theorem 3.1 and
# Corollary 3.3, under the paper's own Definition 2.3:
#   F_{n,D} = sum over weakly increasing i_1<=...<=i_n with i_j < i_{j+1} IF j in D
# (rises REQUIRED at D, extra rises allowed) -- Gessel's classical fundamental
# quasisymmetric function, index set co(D) (partial sums of D-runs), with
# co(empty) = (n).
#
# Verified here:
#   V1  Gessel Theorem 2.4: [M_beta] s_lambda = #{T in SYT(lambda): Des(T) <= pd(beta)}
#   V2  hook property: SYT(n-k,1^k) descent sets = all k-subsets of [n-1]
#   V3  collapse claim: c_D(n,ell) = A(n-|D|,ell,q) for EVERY D (cr preserved)
#   V4  corrected theorem: f(n,ell) = sum_k A(n-k,ell,q) s_(n-k,1^k)   [M-basis check]
#   V5  paper theorem 3.1 with T_{n-1-k,n-k-ell}(q): table of agreements/failures
#   V6  Kostka-solve cross-check of the Schur expansion (independent path)
#   V7  NC: c^NC_D = A^NC(n-|D|,ell); Q_n(NC) = sum_k A^NC(n-k,ell) s_hook;
#       total count = Narayana; comparison with paper corollary
#   V8  inclusion-exclusion: A(m,ell,q) = sum_j (-1)^j C(m-1,j) H(m-j,ell,q)
#   V9  a(m,ell) = sum_j (-1)^j C(m-1,j) S(m-j,ell); a(m) = Bell(m-1)
#   V10 A vs T first-difference table

from itertools import combinations
from collections import defaultdict
from fractions import Fraction
from math import comb
import json, sys

# ---------------- set partitions of [n] (restricted growth strings) -------------
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

def nblocks(rgs):
    return max(rgs) + 1 if rgs else 0

def des_set(rgs):
    return tuple(i for i in range(1, len(rgs)) if rgs[i-1] == rgs[i])

def arcs_cr(rgs):
    """arcs = consecutive block elements with gap >= 2 (adjacent integers are
    straight lines, never counted); cr = number of crossing arc pairs (Def 2.6)."""
    d = defaultdict(list)
    for i, b in enumerate(rgs, 1):
        d[b].append(i)
    ars = []
    for _, blk in sorted(d.items()):
        for u, v in zip(blk, blk[1:]):
            if v - u >= 2:
                ars.append((u, v))
    cr = 0
    for (a, b), (c, e) in combinations(ars, 2):
        if a > c: (a, b), (c, e) = (c, e), (a, b)
        if a < c < b < e:
            cr += 1
    return cr, ars

def is_noncrossing(rgs):
    d = defaultdict(list)
    for i, b in enumerate(rgs, 1):
        d[b].append(i)
    pos = {}
    for bi, (_, blk) in enumerate(sorted(d.items())):
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

# ---------------- polynomial helpers (dict degree -> coeff) ---------------------
def padd(p, q):
    r = dict(p)
    for k, v in q.items(): r[k] = r.get(k, 0) + v
    return {k: v for k, v in r.items() if v}

def pscale(p, c):
    return {k: v*c for k, v in p.items()} if c else {}

def pmul(p, q):
    r = defaultdict(int)
    for k1, v1 in p.items():
        for k2, v2 in q.items(): r[k1+k2] += v1*v2
    return dict(r)

def peq(p, q):
    return {k: v for k, v in p.items() if v} == {k: v for k, v in q.items() if v}

def pstr(p):
    if not p: return "0"
    return " + ".join(f"{v}q^{k}" if k else f"{v}" for k, v in sorted(p.items()))

# ---------------- Touchard-Riordan: T_{a,b}(q) = C(a,2b) T_{2b}(q) --------------
def TR_perfect(m):
    """T_{2m}(q) as coeff list (Touchard-Riordan formula)."""
    size = m*(m-1)//2 + 1
    N = [0]*(size + m + 2)
    for i in range(m+1):
        hi = comb(2*m, m-i)
        lo = comb(2*m, m-i-1) if m-i-1 >= 0 else 0
        e = i*(i+1)//2
        N[e] += (-1)**i * (hi - lo)
    den = [(-1)**k * comb(m, k) for k in range(m+1)]
    T = []
    for kk in range(size):
        acc = N[kk]
        for j in range(1, min(m, kk)+1):
            acc -= den[j]*T[kk-j]
        T.append(acc)
    while T and T[-1] == 0: T.pop()
    return T or [0]

def T_ab(a, b):
    if b < 0 or 2*b > a: return {0: 0}
    if b == 0: return {0: 1}
    base = TR_perfect(b)
    return {k: comb(a, 2*b)*v for k, v in enumerate(base)}

# ---------------- M-basis machinery ----------------------------------------------
def pd_of(beta):
    """partial sums of composition beta, excluding the total."""
    s, out = 0, []
    for x in beta[:-1]:
        s += x; out.append(s)
    return tuple(out)

def mono_key(beta):
    """composition beta -> partition (sorted desc) for m_mu keying."""
    return tuple(sorted(beta, reverse=True))

def kostka_row(lam, n):
    """K[lam][beta] for all compositions beta of n, via SSYT brute force.
    Keyed by the EXACT content vector (c_1,...,c_k): beta_i = number of entries
    equal to i. Fill cells in row-major order; rows weakly increase, columns
    strictly increase; content must use every value 1..max (no gaps)."""
    cells = [(r, c) for r, row in enumerate(lam) for c in range(row)]
    N = len(cells)
    counts = defaultdict(int)
    def rec(idx, assign):
        if idx == N:
            mult = defaultdict(int)
            for v in assign: mult[v] += 1
            if all(mult[v] > 0 for v in range(1, max(mult) + 1)):
                counts[tuple(mult[v] for v in range(1, max(mult) + 1))] += 1
            return
        r, c = cells[idx]
        lo = 1
        if c > 0: lo = max(lo, assign[idx - 1])             # left neighbor
        if r > 0 and c < lam[r - 1]:
            lo = max(lo, assign[idx - lam[r - 1]] + 1)      # top neighbor
        for v in range(lo, n + 1):
            assign.append(v)
            rec(idx + 1, assign)
            assign.pop()
    rec(0, [])
    return dict(counts)

def kostka_merged(krow):
    """Merge exact-content Kostka over compositions with the same sort."""
    m = defaultdict(int)
    for beta, cnt in krow.items():
        m[mono_key(beta)] += cnt
    return dict(m)

def syt_des_sets(shape):
    cells = [(r, c) for r, row in enumerate(shape) for c in range(row)]
    N = len(cells)
    out = []
    def rec(assign):
        if len(assign) == N:
            pos = {v: cell for cell, v in assign.items()}
            out.append(tuple(i for i in range(1, N) if pos[i+1][0] > pos[i][0]))
            return
        v = len(assign) + 1
        for cell in cells:
            if cell in assign: continue
            r, c = cell
            if r > 0 and (r-1, c) not in assign: continue
            if c > 0 and (r, c-1) not in assign: continue
            assign[cell] = v
            rec(assign)
            del assign[cell]
    rec({})
    return out

# ---------------- main ------------------------------------------------------------
def main(nmax=7, json_out="verify3_results.json"):
    log = []
    def say(s):
        print(s); log.append(s)

    say(f"=== verify3: definitive check, paper convention (F: rises required at D) ===")

    # ---------- per-n data ----------
    data = {}   # n -> dict(pl=..., bydes=D->crdist, ncbydes, ...)
    for n in range(1, nmax+1):
        pl = list(set_partitions(n))
        bydes = defaultdict(lambda: defaultdict(int))       # D -> cr -> count
        ncbydes = defaultdict(lambda: defaultdict(int))
        for p in pl:
            cr, _ = arcs_cr(p)
            D = des_set(p)
            bydes[D][cr] += 1
            if is_noncrossing(p):
                ncbydes[D][cr] += 1
        data[n] = dict(pl=pl, bydes={k: dict(v) for k, v in bydes.items()},
                       ncbydes={k: dict(v) for k, v in ncbydes.items()})
    say(f"partition tables built for n<={nmax}")

    # ---------- V2: hook SYT descent sets ----------
    ok2 = True
    for n in range(2, nmax+1):
        for k in range(0, n):
            got = sorted(set(syt_des_sets([n-k] + [1]*k)))
            want = sorted({tuple(sorted(c)) for c in combinations(range(1, n), k)})
            if got != want:
                ok2 = False; say(f"V2 FAIL n={n} k={k}")
    say(f"V2 hook SYT descent sets = all k-subsets: {'PASS' if ok2 else 'FAIL'} (n<={nmax})")

    # ---------- V3: collapse claim c_D = A(n-|D|,ell) ----------
    # A(m,ell,q) = c_emptyset(m,ell): adjacency-free partitions of [m] (Des = ()).
    A = {}   # (m,ell) -> crdist
    for m in range(1, nmax+1):
        for ell in range(1, m+1):
            d = defaultdict(int)
            for p in data[m]["pl"]:
                if nblocks(p) == ell and des_set(p) == ():
                    cr, _ = arcs_cr(p)
                    d[cr] += 1
            A[(m, ell)] = dict(d)
    ok3 = True; bad3 = []
    # redo properly: bydes must be per (n,ell)
    bydes_le = {}   # (n,ell) -> D -> crdist
    ncbydes_le = {}
    for n in range(1, nmax+1):
        for ell in range(1, n+1):
            bd = defaultdict(lambda: defaultdict(int))
            nbd = defaultdict(lambda: defaultdict(int))
            for p in data[n]["pl"]:
                if nblocks(p) != ell: continue
                cr, _ = arcs_cr(p)
                D = des_set(p)
                bd[D][cr] += 1
                if is_noncrossing(p):
                    nbd[D][cr] += 1
            bydes_le[(n, ell)] = {k: dict(v) for k, v in bd.items()}
            ncbydes_le[(n, ell)] = {k: dict(v) for k, v in nbd.items()}
    for (n, ell), bd in bydes_le.items():
        for D, crdist in bd.items():
            if crdist != A.get((n - len(D), ell), {}):
                ok3 = False; bad3.append((n, ell, D, crdist, A.get((n-len(D), ell))))
    say(f"V3 collapse c_D = A(n-|D|,ell) for every D: {'PASS' if ok3 else 'FAIL'} (n<={nmax})"
        + ("" if ok3 else f"  first failures: {bad3[:3]}"))

    # ---------- V1: Gessel theorem in M-basis ----------
    ok1 = True; bad1 = []
    kostka = {}   # (n,lam) -> mu -> K
    for n in range(1, nmax+1):
        lams = sorted(_partitions(n))
        for lam in lams:
            K = kostka_row(list(lam), n)
            kostka[(n, lam)] = K
            des_list = syt_des_sets(list(lam))
            # for every composition beta of n: #{T: Des(T) subseteq pd(beta)} == K[sort(beta)]
            for beta in _compositions(n):
                lhs = sum(1 for D in des_list if set(D) <= set(pd_of(beta)))
                rhs = K.get(beta, 0)
                if lhs != rhs:
                    ok1 = False; bad1.append((n, lam, beta, lhs, rhs))
    say(f"V1 Gessel Thm 2.4 in M-basis (#{'{'}T: Des<=pd(beta){'}'} = K_(lam,mu)): "
        f"{'PASS' if ok1 else 'FAIL'} (n<={nmax})" + ("" if ok1 else f"  first: {bad1[:3]}"))

    # ---------- V4/V5: f(n,ell): direct M-coeffs vs corrected-hook vs paper-hook ----
    # direct: [M_beta] f = sum_{D subseteq pd(beta)} c_D    (monomial x^a in F_{n,D} iff pd(beta) superseteq D)
    # corrected hook RHS: sum_k A(n-k,ell,q) * s_(n-k,1^k):  [M_beta] = sum_k A(n-k,ell,q) * K[hook_k][mu(beta)]
    # paper RHS: same with T_(n-1-k, n-k-ell)(q).
    results = []
    ok4 = True; bad4 = []
    diff_table = []
    for n in range(1, nmax+1):
        for ell in range(1, n+1):
            bd = bydes_le[(n, ell)]
            cD = {D: {cr: v for cr, v in d.items()} for D, d in bd.items()}   # cr -> count (poly in q)
            # direct M-coeffs, keyed by the EXACT composition beta:
            #   [M_beta] f = sum_{D subseteq pd(beta)} c_D
            comps = list(_compositions(n))
            direct = {}
            for beta in comps:
                S = set(pd_of(beta))
                acc = defaultdict(int)
                for D, crdist in cD.items():
                    if set(D) <= S:
                        for cr, v in crdist.items(): acc[cr] += v
                direct[beta] = dict(acc)
            # corrected hook RHS: [M_beta] = sum_k A(n-k,ell,q) * K_{hook_k, beta}
            hookM = {}
            for beta in comps:
                acc = defaultdict(int)
                for k in range(0, n):
                    Av = A.get((n - k, ell), {})
                    kv = kostka[(n, tuple([n - k] + [1] * k))].get(beta, 0)
                    if kv:
                        for cr, v in Av.items(): acc[cr] += kv * v
                hookM[beta] = dict(acc)
            match4 = all(peq(direct.get(b, {}), hookM.get(b, {})) for b in comps)
            if not match4:
                ok4 = False; bad4.append((n, ell))
            # paper RHS M-coeffs: same hooks with T_(n-1-k, n-k-ell)(q)
            paperM = {}
            for beta in comps:
                acc = defaultdict(int)
                for k in range(0, n):
                    Tv = T_ab(n - 1 - k, n - k - ell)
                    kv = kostka[(n, tuple([n - k] + [1] * k))].get(beta, 0)
                    if kv:
                        for cr, v in Tv.items(): acc[cr] += kv * v
                paperM[beta] = dict(acc)
            match5 = all(peq(direct.get(b, {}), paperM.get(b, {})) for b in comps)
            # hook-coefficient comparison table: true A(n-k,ell,q) vs paper T_{n-1-k,n-k-ell}(q)
            row = []
            for k in range(0, n):
                Av = A.get((n-k, ell), {})
                Tv = T_ab(n-1-k, n-k-ell)
                row.append((k, Av, Tv, peq(Av, Tv)))
            results.append(dict(n=n, ell=ell, direct=direct, hookM=hookM, paperM=paperM,
                                match_corrected=match4, match_paper=match5,
                                hooks=[dict(k=k, A=A1, T=T1, eq=eq) for k, A1, T1, eq in row]))
            for k, A1, T1, eq in row:
                if not eq:
                    diff_table.append((n, ell, k, A1, T1))
    say(f"V4 corrected theorem f = sum_k A(n-k,ell,q) s_hook: {'PASS' if ok4 else 'FAIL'} (n<={nmax})"
        + ("" if ok4 else f"  failures: {bad4}"))
    n_paper_fail = sum(1 for r in results if not r["match_paper"])
    say(f"V5 paper Theorem 3.1 (T-coeffs): {n_paper_fail} of {len(results)} (n,ell) pairs FAIL (n<={nmax})")
    say("V5 first differing (n,ell,k): true A vs paper T:")
    shown = 0
    for n, ell, k, A1, T1 in diff_table:
        say(f"    (n={n}, ell={ell}, k={k}): A({n-k},{ell},q) = {pstr(A1)}   vs   T_{{{n-1-k},{n-k-ell}}}(q) = {pstr(T1)}")
        shown += 1
        if shown >= 8: break

    # ---------- V6: Kostka-solve cross-check ----------
    ok6 = True
    for n in range(1, nmax+1):
        lams = sorted(_partitions(n))
        kostkaM = {lam: kostka_merged(kostka[(n, lam)]) for lam in lams}
        for ell in range(1, n+1):
            # merged direct M-coeffs over partitions mu (exact-composition -> sort)
            directC = next(r for r in results if r["n"] == n and r["ell"] == ell)["direct"]
            direct = defaultdict(lambda: defaultdict(int))
            for beta, d in directC.items():
                for cr, v in d.items(): direct[mono_key(beta)][cr] += v
            direct = {k: dict(v) for k, v in direct.items()}
            # solve K a = direct-m using all shapes, triangular in dominance order
            # dominance order: mu dominates nu iff partial sums >=; solve from dominant end.
            def dom_key(mu): return tuple(-s for s in _cumsum(mu))
            lams_sorted = sorted(lams, key=dom_key)   # most dominant first
            coeff = {}
            for mu in lams_sorted:
                # coeff[mu] * K[mu][mu] + sum_{nu dominates mu, nu != mu} coeff[nu]*K[nu][mu] = direct[mu]
                target = dict(direct.get(mu, {}))
                for nu in lams_sorted:
                    if nu == mu: continue
                    kv = kostkaM[nu].get(mu, 0)
                    if kv and nu in coeff:
                        # subtract coeff[nu]*kv  (polynomial times int)
                        for cr, v in coeff[nu].items():
                            target[cr] = target.get(cr, 0) - kv*v
                target = {k: v for k, v in target.items() if v}
                # diagonal of the merged Kostka matrix: K~[mu][mu] = #{comps sorting to mu}
                diag = kostkaM[mu].get(mu, 0)
                if diag == 0:
                    ok6 = False; say(f"V6 FAIL n={n} ell={ell} lam={lam}: zero diagonal")
                    coeff[mu] = {}
                else:
                    coeff[mu] = {k: Fraction(v, diag) for k, v in target.items()}
                    if any(c.denominator != 1 for c in coeff[mu].values()):
                        ok6 = False; say(f"V6 FAIL n={n} ell={ell} lam={lam}: non-integer solve {coeff[mu]}")
            # compare with corrected hook formula restricted to hooks + zeros elsewhere
            hookcoeff = {}
            for k in range(0, n):
                hookcoeff[tuple([n-k]+[1]*k)] = A.get((n-k, ell), {})
            for lam in lams_sorted:
                want = hookcoeff.get(lam, {})
                if not peq(coeff.get(lam, {}), want):
                    ok6 = False
                    say(f"V6 FAIL n={n} ell={ell} lam={lam}: solve={coeff.get(lam)} want={want}")
    say(f"V6 Kostka-solve Schur expansion == corrected hook formula: {'PASS' if ok6 else 'FAIL'} (n<=6)")

    # ---------- V7: NC ----------
    # A^NC(m,ell) = # NC partitions of [m], ell blocks, Des = () (adjacency-free)
    ANC = {}
    for m in range(1, nmax+1):
        for ell in range(1, m+1):
            cnt = sum(1 for p in data[m]["pl"] if nblocks(p) == ell and des_set(p) == () and is_noncrossing(p))
            ANC[(m, ell)] = cnt
    ok7 = True; bad7 = []
    for (n, ell), nbd in ncbydes_le.items():
        for D, crdist in nbd.items():
            if crdist.get(0, 0) != ANC.get((n - len(D), ell), 0):
                ok7 = False; bad7.append((n, ell, D, crdist, ANC.get((n-len(D), ell))))
    say(f"V7a NC collapse c^NC_D = A^NC(n-|D|,ell) for every D: {'PASS' if ok7 else 'FAIL'} (n<={nmax})"
        + ("" if ok7 else f"  first: {bad7[:3]}"))
    # Narayana count check
    ok7b = True
    for n in range(1, nmax+1):
        for ell in range(1, n+1):
            total = sum(sum(d.values()) for d in ncbydes_le[(n, ell)].values())
            formula = sum(comb(n-1, k)*ANC.get((n-k, ell), 0) for k in range(0, n))
            nar = _narayana(n, ell)
            if total != nar or formula != nar:
                ok7b = False; say(f"V7b FAIL n={n} ell={ell}: total={total} formula={formula} narayana={nar}")
    say(f"V7b NC total = sum_k C(n-1,k) A^NC(n-k,ell) = Narayana(n,ell): {'PASS' if ok7b else 'FAIL'} (n<={nmax})")

    # ---------- V7c: Motzkin identification ----------
    # Paper Cor 3.3(3) claims Q_n(NC(n,ell)) = sum_k M_{n-1-k, n-ell-k} s_hook with
    # M_{a,b} = # noncrossing matchings on a points with b chords (Motzkin triangle).
    # Our collapse gives Q_n(NC(n,ell)) = sum_k A^NC(n-k, ell) s_hook, so item 3 is
    # TRUE iff A^NC(m, ell) = nc_match_count(m-1, m-ell) for all m, ell.
    from functools import lru_cache
    def nc_match_count(a, b):
        """# noncrossing matchings (partial, disjoint chords) on [a] with exactly b chords."""
        @lru_cache(maxsize=None)
        def go(S):
            if not S: return {0: 1}
            out = defaultdict(int)
            for k, v in go(S[1:]).items(): out[k] += v          # min point unmatched
            for qi in range(1, len(S)):                          # min point paired with S[qi]
                inside, outside = S[1:qi], S[qi+1:]
                for k1, v1 in go(inside).items():
                    for k2, v2 in go(outside).items():
                        out[k1 + k2 + 1] += v1 * v2
            return dict(out)
        return go(tuple(range(1, a + 1))).get(b, 0)
    def _catalan(b):
        return comb(2*b, b) // (b + 1)
    ok7c = True; bad7c = []
    for m in range(1, nmax+1):
        for ell in range(1, m+1):
            got = ANC.get((m, ell), 0)
            want = nc_match_count(m - 1, m - ell)
            if got != want:
                ok7c = False; bad7c.append((m, ell, got, want))
            # structural identity: nc_match_count(a,b) = C(a,2b) Catalan(b)
            a_, b_ = m - 1, m - ell
            if nc_match_count(a_, b_) != comb(a_, 2*b_) * _catalan(b_):
                ok7c = False; bad7c.append(("cat", a_, b_))
    say(f"V7c A^NC(m,ell) = #NC matchings on m-1 pts with m-ell chords = C(m-1,2b)Cat_b: "
        f"{'PASS' if ok7c else 'FAIL'} (n<={nmax})" + ("" if ok7c else f"  first: {bad7c[:3]}"))
    # V7d: sum_ell A^NC(m,ell) = Motzkin(m-1)  [Cor 3.3(4) coefficient totals]
    def _motzkin(n):
        memo = [0]*(n+2); memo[0] = 1
        for i in range(1, n+1):
            memo[i] = memo[i-1] + sum(memo[j]*memo[i-2-j] for j in range(0, i-1))
        return memo[n]
    ok7d = True; bad7d = []
    for m in range(1, nmax+1):
        tot = sum(ANC.get((m, ell), 0) for ell in range(1, m+1))
        if tot != _motzkin(m - 1):
            ok7d = False; bad7d.append((m, tot, _motzkin(m-1)))
    say(f"V7d sum_ell A^NC(m,ell) = Motzkin(m-1): {'PASS' if ok7d else 'FAIL'} (n<={nmax})"
        + ("" if ok7d else f"  first: {bad7d[:3]}"))

    # ---------- V8/V9 ----------
    ok8 = True
    for m in range(1, nmax+1):
        for ell in range(1, m+1):
            # H(c,ell,q) = total crossing poly of Set(c,ell)
            ie = defaultdict(int)
            for j in range(0, m):
                c = m - j
                if c < ell: continue
                Hc = defaultdict(int)
                for p in data[c]["pl"]:
                    if nblocks(p) == ell:
                        cr, _ = arcs_cr(p); Hc[cr] += 1
                sgn = (-1)**j * comb(m-1, j)
                for cr, v in Hc.items(): ie[cr] += sgn*v
            ie = {k: v for k, v in ie.items() if v}
            if not peq(ie, A.get((m, ell), {})):
                ok8 = False; say(f"V8 FAIL m={m} ell={ell}: IE={ie} A={A.get((m,ell))}")
    say(f"V8 inclusion-exclusion A(m,ell,q) = sum_j (-1)^j C(m-1,j) H(m-j,ell,q): {'PASS' if ok8 else 'FAIL'} (n<={nmax})")
    ok9 = True
    for m in range(1, nmax+1):
        tot = sum(ANC_and_A_counts(m, ell, A) for ell in range(1, m+1))
        bell = _bell(m-1)
        if tot != bell:
            ok9 = False; say(f"V9 FAIL m={m}: sum_ell a(m,ell)={tot} Bell(m-1)={bell}")
        for ell in range(1, m+1):
            ie = sum((-1)**j * comb(m-1, j) * _stirling(m-j, ell) for j in range(0, m))
            a_tot = sum(A.get((m, ell), {}).values())
            if ie != a_tot:
                ok9 = False; say(f"V9b FAIL m={m} ell={ell}: IE={ie} a={a_tot}")
    say(f"V9 a(m,ell)=sum_j (-1)^j C(m-1,j) S(m-j,ell) and sum_ell a(m,ell)=Bell(m-1): {'PASS' if ok9 else 'FAIL'}")

    # ---------- save ----------
    out = dict(log=log, results=[{k: v for k, v in r.items() if k in ("n", "ell", "match_corrected", "match_paper", "hooks")} for r in results],
               A={f"{m},{ell}": v for (m, ell), v in A.items()},
               ANC={f"{m},{ell}": v for (m, ell), v in ANC.items()})
    with open(json_out, "w") as f:
        json.dump(out, f, indent=1, default=str)
    say(f"saved {json_out}")
    return results, A, ANC

def _partitions(n):
    if n == 0:
        yield (); return
    for p in range(n, 0, -1):
        for rest in _partitions(n - p):
            if not rest or p >= rest[0]:
                yield (p,) + rest

def _compositions(n):
    if n == 0:
        yield (); return
    for p in range(1, n+1):
        for rest in _compositions(n - p):
            yield (p,) + rest

def _cumsum(t):
    out, s = [], 0
    for x in t: s += x; out.append(s)
    return tuple(out)

def _narayana(n, k):
    if k < 1 or k > n: return 0
    return comb(n, k)*comb(n, k-1)//n

def _stirling(n, k):
    if n == 0 and k == 0: return 1
    if n == 0 or k == 0: return 0
    return k*_stirling(n-1, k) + _stirling(n-1, k-1)

def _bell(n):
    return sum(_stirling(n, k) for k in range(0, n+1))

def ANC_and_A_counts(m, ell, A):
    return sum(A.get((m, ell), {}).values())

if __name__ == "__main__":
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    main(nmax)
