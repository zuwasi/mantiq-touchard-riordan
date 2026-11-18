(* ::Package:: *)

(* SetPartitionQSym.wl -- MantiQ Factor companion package for
   Bagno-Garber, "Touchard-Riordan Polynomials and Schur-positivity of Set
   Partitions" (arXiv:2606.13149).

   Conventions follow the paper:
   - Set partitions of [n] are encoded as restricted growth strings (rgs):
     rgs[[j]] = block index of j, blocks ordered by smallest element.
   - Des[pi] = { i in [n-1] : i and i+1 lie in the same block }  (Sec. 3)
   - cr(pi) = number of crossing pairs of arcs, where arcs are consecutive
     elements of a block with gap >= 2 (adjacent integers give straight lines,
     Def. 2.6).
   - F[n, D] = Gessel's fundamental quasisymmetric function under the paper's
     Definition 2.3 (rises REQUIRED at positions in D, extra rises allowed):
     sum over weakly increasing i_1<=...<=i_n with a strict rise at every j in D,
     returned as a polynomial in x_1..x_n.

   Main result implemented here (correcting Theorem 3.1 of the paper):
     Sum_{pi in Set(n,l)} q^cr(pi) F_{n,Des(pi)}
       = Sum_{k=0}^{n-1} A(n-k,l,q) s_(n-k,1^k),
   where A(m,l,q) is the crossing polynomial of adjacency-free (Des = {})
   partitions of [m] into l blocks. The paper's coefficient T_{n-1-k,n-k-l}(q)
   agrees only when every adjacency-free partition comes from the chord-diagram
   map, which fails from n = 4, l = 2 on. The noncrossing analogue of the paper
   (Corollary 3.3 (3)-(4)) is correct: A^NC(m,l) equals the Motzkin triangle
   number counting noncrossing matchings on m-1 points with m-l chords. *)

BeginPackage["SetPartitionQSym`"]

SetPartitionsRGS::usage = "SetPartitionsRGS[n] returns all set partitions of [n] as restricted growth strings (lists of block indices)."
BlocksOfPartition::usage = "BlocksOfPartition[rgs] returns the blocks of the partition as sorted lists."
DescentSetOfPartition::usage = "DescentSetOfPartition[rgs] returns Des(pi) = {i : i and i+1 lie in the same block}."
ArcsOfPartition::usage = "ArcsOfPartition[rgs] returns the arcs {u,v}, the consecutive block elements with gap >= 2."
CrossingNumberPartition::usage = "CrossingNumberPartition[rgs] returns cr(pi), the number of crossing arc pairs."
NoncrossingQ::usage = "NoncrossingQ[rgs] tests whether the partition is noncrossing."
FundamentalQSym::usage = "FundamentalQSym[n, D, vars] returns F_{n,D} as a polynomial in vars (paper Def. 2.3: rises required at D)."
QSet::usage = "QSet[n, ell, q, vars] returns Sum q^cr(pi) F_{n,Des(pi)} over pi in Set(n, ell)."
QSetNC::usage = "QSetNC[n, ell, q, vars] is the noncrossing analogue over NC(n, ell)."
HookSchur::usage = "HookSchur[a, b, vars] returns the Schur polynomial s_(a,1^b) via s_(a,1^b) = Sum_j (-1)^j h_{a+j} e_{b-j}."
ACoeff::usage = "ACoeff[m, ell, q] returns the corrected coefficient: the crossing polynomial of adjacency-free (Des = {}) partitions of [m] into ell blocks."
ANCCount::usage = "ANCCount[m, ell] returns the number of noncrossing adjacency-free partitions of [m] into ell blocks."
TouchardRiordanT::usage = "TouchardRiordanT[a, b, q] returns the paper's coefficient T_{a,b}(q) = C(a,2b) T_{2b}(q), or 0 when 2b > a."
CorrectedRHS::usage = "CorrectedRHS[n, ell, q, vars] returns Sum_k ACoeff[n-k, ell, q] s_(n-k,1^k): the corrected right side of Theorem 3.1."
PaperRHS::usage = "PaperRHS[n, ell, q, vars] returns Sum_k TouchardRiordanT[n-1-k, n-k-ell, q] s_(n-k,1^k): the paper's right side of Theorem 3.1."
SYTDescentSets::usage = "SYTDescentSets[shape] returns the descent sets of all standard Young tableaux of the given shape."
NCMatchCount::usage = "NCMatchCount[a, b] returns the number of noncrossing matchings on a points with exactly b chords (a Motzkin triangle number)."
MotzkinNumber::usage = "MotzkinNumber[n] returns the n-th Motzkin number."
PartitionArcDiagram::usage = "PartitionArcDiagram[rgs] draws the arc diagram of a set partition."
RunValidationSuite::usage = "RunValidationSuite[nmax:6] runs the full verification suite and returns an association of results."
RunProofAudit::usage = "RunProofAudit[] runs the paper-specific audit (corrected theorem, counterexample census, boundary cases, Motzkin identification) and returns {lines, passed, failed, warned}."
polyDiffZero::usage = "polyDiffZero[a, b] tests whether the polynomial a - b simplifies to 0."

polyDiffZero[a_, b_] := Simplify[Expand[a - b]] === 0

Begin["`Private`"]

(* ---------- set partitions as restricted growth strings ---------- *)

SetPartitionsRGS[n_] := SetPartitionsRGS[n] = Module[{step},
  step[w_] := Table[Append[w, b], {b, 0, If[w === {}, 0, Max[w] + 1]}];
  Nest[Join @@ (step /@ #) &, {{}}, n]
]

BlocksOfPartition[rgs_List] := With[
  {d = GroupBy[Range[Length[rgs]], rgs[[#]] &]},
  Table[Sort[d[b]], {b, 0, Max[rgs]}]
]

DescentSetOfPartition[rgs_List] := Select[Range[Length[rgs] - 1],
  rgs[[#]] == rgs[[# + 1]] &]

ArcsOfPartition[rgs_List] := Catenate[
  Select[Partition[#, 2, 1], #[[2]] - #[[1]] >= 2 &] & /@
    BlocksOfPartition[rgs]]

CrossingNumberPartition[rgs_List] := Module[
  {arcs = ArcsOfPartition[rgs], pairs, cnt = 0, ab, ce},
  pairs = Subsets[arcs, {2}];
  Do[
    {ab, ce} = SortBy[pairs[[i]], First];
    If[ab[[1]] < ce[[1]] < ab[[2]] < ce[[2]], cnt++],
    {i, Length[pairs]}];
  cnt
]

NoncrossingQ[rgs_List] := Catch[Module[{n = Length[rgs], i, j, k, l},
  For[i = 1, i <= n - 2, i++,
    For[k = i + 1, k <= n - 1, k++,
      If[rgs[[i]] === rgs[[k]],
        For[j = i + 1, j <= k - 1, j++,
          For[l = k + 1, l <= n, l++,
            If[rgs[[j]] === rgs[[l]] && rgs[[j]] =!= rgs[[i]],
              Throw[False]]]]]]];
  True
]]

(* ---------- fundamental quasisymmetric function, paper Def. 2.3 ---------- *)

(* F_{n,D} = Sum over compositions beta of n whose partial sums (excluding the
   total) contain D, of the monomial quasisymmetric function M_beta restricted
   to the given variables: monomials x_{i1}^b1 ... x_{ik}^bk over strictly
   increasing index sets. This is exactly the paper's sum over weakly
   increasing i_1<=...<=i_n with a strict rise at every j in D. *)

partialSums[comp_List] := Most[Accumulate[comp]]

riseOK[comp_List, D_List] := SubsetQ[partialSums[comp], D]

monomialSym[beta_List, vars_List] := Sum[
  Times @@ Table[vars[[idx[[i]]]]^beta[[i]], {i, Length[beta]}],
  {idx, Subsets[Range[Length[vars]], {Length[beta]}]}]

FundamentalQSym[n_Integer, D_List, vars_List] :=
  Total[monomialSym[#, vars] & /@
    Select[Catenate[Permutations /@ IntegerPartitions[n]], riseOK[#, D] &]]

(* ---------- the main object and the two right sides ---------- *)

QSet[n_Integer, ell_Integer, q_Symbol, vars_List] := Module[
  {parts = Select[SetPartitionsRGS[n], Max[#] + 1 == ell &]},
  Total[(q^CrossingNumberPartition[#])*
      FundamentalQSym[n, DescentSetOfPartition[#], vars] & /@ parts]
]

QSetNC[n_Integer, ell_Integer, q_Symbol, vars_List] := Module[
  {parts = Select[SetPartitionsRGS[n], Max[#] + 1 == ell && NoncrossingQ[#] &]},
  Total[(q^CrossingNumberPartition[#])*
      FundamentalQSym[n, DescentSetOfPartition[#], vars] & /@ parts]
]

(* ---------- symmetric-function helpers ---------- *)

(* h_r as the sum of all monomials of degree r: group integer partitions of r,
   pad with zeros, and sum over distinct exponent arrangements. *)
monoExps[exp_List, vars_List] :=
  Total[Times @@ (vars^#) & /@ DeleteDuplicates[Permutations[exp]]]

homogeneousH[r_Integer, vars_List] := If[r == 0, 1,
  Total[monoExps[Join[#, ConstantArray[0, Length[vars] - Length[#]]], vars] & /@
    IntegerPartitions[r]]]

elementaryE[r_Integer, vars_List] := If[r == 0, 1, SymmetricPolynomial[r, vars]]

HookSchur[a_Integer, b_Integer, vars_List] :=
  Expand[Sum[(-1)^j homogeneousH[a + j, vars] elementaryE[b - j, vars], {j, 0, b}]]

(* ---------- corrected and paper coefficients ---------- *)

adjacencyFreePartitions[m_Integer, ell_Integer] :=
  Select[SetPartitionsRGS[m],
    Max[#] + 1 == ell && DescentSetOfPartition[#] === {} &]

ACoeff[m_Integer, ell_Integer, q_Symbol] := Module[{d},
  If[ell > m, 0,
    d = Tally[CrossingNumberPartition /@ adjacencyFreePartitions[m, ell]];
    If[d === {}, 0, Total[q^#[[1]] #[[2]] & /@ d]]]
]

ANCCount[m_Integer, ell_Integer] :=
  If[ell > m, 0,
    Length[Select[adjacencyFreePartitions[m, ell], NoncrossingQ]]]

TouchardRiordanT[a_Integer, b_Integer, q_Symbol] := Module[{m = b, num},
  If[b < 0 || 2 b > a, 0,
    If[b == 0, 1,
      num = Sum[(-1)^i (Binomial[2 m, m - i] - Binomial[2 m, m - i - 1])*
          q^Binomial[i + 1, 2], {i, 0, m}];
      Binomial[a, 2 b] Factor[num/(1 - q)^m]]]
]

CorrectedRHS[n_Integer, ell_Integer, q_Symbol, vars_List] :=
  Sum[ACoeff[n - k, ell, q] HookSchur[n - k, k, vars], {k, 0, n - 1}]

PaperRHS[n_Integer, ell_Integer, q_Symbol, vars_List] :=
  Sum[TouchardRiordanT[n - 1 - k, n - k - ell, q]*
      HookSchur[n - k, k, vars], {k, 0, n - 1}]

(* ---------- standard Young tableaux (descent sets) ---------- *)

SYTDescentSets[shape_List] := Module[
  {cells, n, out, rec, cellOf},
  cells = Catenate[Table[{r, c}, {r, Length[shape]}, {c, shape[[r]]}]];
  n = Length[cells];
  out = {};
  rec[assign_Association, v_] := If[v > n,
    cellOf = Association[Reverse /@ Normal[assign]];
    AppendTo[out, Select[Range[n - 1],
      cellOf[# + 1][[1]] > cellOf[#][[1]] &]],
    Do[
      With[{cell = cells[[i]], r = cells[[i, 1]], c = cells[[i, 2]]},
        If[FreeQ[Keys[assign], cell] &&
          (r == 1 || (KeyExistsQ[assign, {r - 1, c}] && assign[{r - 1, c}] < v)) &&
          (c == 1 || (KeyExistsQ[assign, {r, c - 1}] && assign[{r, c - 1}] <= v)),
          rec[Append[assign, cell -> v], v + 1]]],
      {i, Length[cells]}]];
  rec[<||>, 1];
  out
]

(* ---------- noncrossing matchings and Motzkin numbers ---------- *)

NCMatchCount[a_Integer, b_Integer] := Module[{rec},
  rec[S_List] := rec[S] = If[S === {}, <|0 -> 1|>,
    Module[{out = <||>, kv1, kv2},
      Do[ (* smallest point unmatched *)
        out[kv1[[1]]] = Lookup[out, kv1[[1]], 0] + kv1[[2]],
        {kv1, List @@@ Normal[rec[Rest[S]]]}];
      Do[ (* smallest point paired: chords split into inside and outside *)
        With[{sub1 = rec[S[[2 ;; qi - 1]]], sub2 = rec[S[[qi + 1 ;; -1]]]},
          Do[out[kv1[[1]] + kv2[[1]] + 1] =
              Lookup[out, kv1[[1]] + kv2[[1]] + 1, 0] + kv1[[2]] kv2[[2]],
            {kv1, List @@@ Normal[sub1]}, {kv2, List @@@ Normal[sub2]}]],
        {qi, 2, Length[S]}];
      out]];
  Lookup[rec[Range[a]], b, 0]
]

MotzkinNumber[n_Integer] := Module[{m},
  m = ConstantArray[0, Max[n + 1, 2]];
  m[[1]] = 1;
  If[n >= 1,
    Do[m[[i + 1]] = m[[i]] + Sum[m[[j + 1]] m[[i - 1 - j]], {j, 0, i - 2}], {i, n}]];
  m[[n + 1]]
]

(* ---------- arc diagram drawing ---------- *)

PartitionArcDiagram[rgs_List] := Module[
  {n = Length[rgs], arcs, cols},
  arcs = ArcsOfPartition[rgs];
  cols = ColorData[97] /@ (rgs + 1);
  Graphics[{
    Table[{cols[[j]], PointSize[0.02], Point[{j, 0}],
        Text[Style[j, 10], {j, 0.12}]}, {j, n}],
    {GrayLevel[0.5], Line[{{1, 0}, {n, 0}}]},
    Table[With[{u = a[[1]], v = a[[2]]},
      {cols[[rgs[[u]] + 1]], Circle[{(u + v)/2, 0}, (v - u)/2, {0, Pi}]}],
      {a, arcs}]},
    ImageSize -> 340, PlotRange -> {{0.4, n + 0.6}, {-0.4, 1.4}}]
]

(* ---------- validation suite ---------- *)

RunValidationSuite[nmax_:6] := Module[
  {res = <||>, vars, q = Global`q, ok, fails, bad, n, ell, m},
  (* C1: corrected theorem QSet == CorrectedRHS *)
  fails = {};
  Do[
    vars = Table[Symbol["x" <> ToString[i]], {i, n}];
    If[! polyDiffZero[QSet[n, ell, q, vars], CorrectedRHS[n, ell, q, vars]],
      AppendTo[fails, {n, ell}]],
    {n, 2, Min[nmax, 6]}, {ell, n}];
  res["CorrectedTheorem"] = <|"pass" -> fails === {}, "failures" -> fails|>;
  (* C2: paper Theorem 3.1 *)
  fails = {};
  Do[
    vars = Table[Symbol["x" <> ToString[i]], {i, n}];
    If[! polyDiffZero[QSet[n, ell, q, vars], PaperRHS[n, ell, q, vars]],
      AppendTo[fails, {n, ell}]],
    {n, 2, Min[nmax, 6]}, {ell, n}];
  res["PaperTheorem31"] = <|"pass" -> fails === {}, "failures" -> fails,
    "firstCounterexample" -> If[fails === {}, None, First[fails]]|>;
  (* C3: collapse lemma c_D = A(n-|D|, ell) for every descent set *)
  bad = {};
  Do[
    With[{byD = GroupBy[
        Select[SetPartitionsRGS[n], Max[#] + 1 == ell &],
        DescentSetOfPartition]},
      Do[
        With[{dist = Tally[CrossingNumberPartition /@ byD[D]]},
          If[! polyDiffZero[Total[q^#[[1]] #[[2]] & /@ dist],
              ACoeff[n - Length[D], ell, q]],
            AppendTo[bad, {n, ell, D}]]],
        {D, Keys[byD]}]],
    {n, 2, Min[nmax, 6]}, {ell, n}];
  res["CollapseLemma"] = <|"pass" -> bad === {}, "failures" -> bad|>;
  (* C4: hook SYT descent sets are all k-subsets (Gessel duality, hook case) *)
  bad = {};
  Do[
    With[{got = Sort[SYTDescentSets[Join[{n - k}, Table[1, {k}]]]]},
      With[{want = Sort[Subsets[Range[n - 1], {k}]]},
        If[got =!= want, AppendTo[bad, {n, k}]]]],
    {n, 2, Min[nmax, 6]}, {k, 0, n - 1}];
  res["HookSYTDescentSets"] = <|"pass" -> bad === {}, "failures" -> bad|>;
  (* C5: NC coefficients = Motzkin triangle (paper Cor. 3.3 (3)) *)
  bad = {};
  Do[
    If[ANCCount[n, ell] =!= NCMatchCount[n - 1, n - ell],
      AppendTo[bad, {n, ell}]],
    {n, 2, Min[nmax, 7]}, {ell, n}];
  res["NCMotzkinTriangle"] = <|"pass" -> bad === {}, "failures" -> bad|>;
  (* C6: NC totals: Motzkin (Cor. 3.3 (4)) and Narayana per block count *)
  bad = {};
  Do[
    If[Sum[ANCCount[n, ell], {ell, n}] =!= MotzkinNumber[n - 1],
      AppendTo[bad, {n, "motzkin"}]];
    With[{nar = Binomial[n, ell] Binomial[n, ell - 1]/n},
      If[Sum[ANCCount[n - k, ell] Binomial[n - 1, k], {k, 0, n - 1}] =!= nar,
        AppendTo[bad, {n, ell}]]],
    {n, 2, Min[nmax, 7]}, {ell, n}];
  res["NCTotals"] = <|"pass" -> bad === {}, "failures" -> bad|>;
  (* C7: inclusion-exclusion identity for the corrected coefficients at q = 1 *)
  bad = {};
  Do[
    With[{a1 = Total[CoefficientList[ACoeff[m, ell, q], q]],
        ie = Sum[(-1)^j Binomial[m - 1, j] StirlingS2[m - j, ell], {j, 0, m - 1}]},
      If[a1 =!= ie, AppendTo[bad, {m, ell, a1, ie}]]],
    {m, 2, Min[nmax, 6]}, {ell, m}];
  res["InclusionExclusionBell"] = <|"pass" -> bad === {}, "failures" -> bad|>;
  res
]

(* ---------- proof audit: paper-specific checks A1-A11 ---------- *)

RunProofAudit[] := Module[
  {lines = {}, passed = 0, failed = 0, warned = 0,
   q = Global`q, vars, okAll, badAll, paperFail, okNC, badNC, okMotz, okBell,
   log, audit, warn, n, ell, m, k},
  log[s_] := AppendTo[lines, ToString[s, InputForm]];
  audit[name_, test_, detail_: ""] :=
    If[TrueQ[test], passed++; log["PASS: " <> name],
      failed++; log["FAIL: " <> name <> " -- " <> ToString[detail, InputForm]]];
  warn[name_, detail_] := (warned++;
    log["WARN: " <> name <> " -- " <> ToString[detail, InputForm]]);
  log["== MantiQ proof audit: Bagno-Garber arXiv:2606.13149 =="];
  log["Wolfram: " <> $Version];
  audit["A1 corrected coefficient A(4,2,q) = q",
    ACoeff[4, 2, q] === q, ACoeff[4, 2, q]];
  audit["A2 paper coefficient T_(3,2)(q) = 0",
    TouchardRiordanT[3, 2, q] === 0, TouchardRiordanT[3, 2, q]];
  warn["A3 Theorem 3.1 first counterexample (n,l) = (4,2): coefficient of s_(4)",
    "true q vs paper 0"];
  okAll = True; badAll = {};
  Do[
    vars = Table[Symbol["x" <> ToString[i]], {i, n}];
    If[! polyDiffZero[QSet[n, ell, q, vars], CorrectedRHS[n, ell, q, vars]],
      okAll = False; AppendTo[badAll, {n, ell}]],
    {n, 2, 6}, {ell, n}];
  audit["A4 corrected theorem QSet == Sum_k A(n-k,l,q) s_(n-k,1^k) for n <= 6",
    okAll, badAll];
  paperFail = {};
  Do[
    vars = Table[Symbol["x" <> ToString[i]], {i, n}];
    If[! polyDiffZero[QSet[n, ell, q, vars], PaperRHS[n, ell, q, vars]],
      AppendTo[paperFail, {n, ell}]],
    {n, 2, 6}, {ell, n}];
  audit["A5 paper Theorem 3.1 fails for exactly the recorded (n,l) pairs (n <= 6)",
    paperFail === {{4, 2}, {5, 2}, {5, 3}, {6, 2}, {6, 3}, {6, 4}}, paperFail];
  audit["A6 ell = 1: QSet(n,1) = e_n = s_(1,...,1) (all descents, no crossings)",
    TrueQ[And @@ Table[
      vars = Table[Symbol["x" <> ToString[i]], {i, n}];
      polyDiffZero[QSet[n, 1, q, vars], HookSchur[1, n - 1, vars]],
      {n, 2, 5}]]];
  audit["A7 ell = n: QSet(n,n) = s_(n) (no descents, no crossings)",
    TrueQ[And @@ Table[
      vars = Table[Symbol["x" <> ToString[i]], {i, n}];
      polyDiffZero[QSet[n, n, q, vars], HookSchur[n, 0, vars]],
      {n, 2, 5}]], "ell=n mismatch"];
  okNC = True; badNC = {};
  Do[
    If[ANCCount[m, ell] =!= NCMatchCount[m - 1, m - ell],
      okNC = False; AppendTo[badNC, {m, ell}]],
    {m, 2, 7}, {ell, m}];
  audit["A8 A^NC(m,l) = #NC matchings on m-1 points with m-l chords (Motzkin triangle), m <= 7",
    okNC, badNC];
  okMotz = And @@ Table[
    Sum[ANCCount[m, ell], {ell, m}] === MotzkinNumber[m - 1], {m, 2, 7}];
  audit["A9 Sum_l A^NC(m,l) = Motzkin(m-1), m <= 7", okMotz];
  audit["A10 collapse c_D(4,2) = A(4-|D|,2,q) for D = {}, {1}, {2}, {3}",
    TrueQ[And @@ (Function[{d},
        With[{parts = Select[SetPartitionsRGS[4],
            Max[#] + 1 == 2 && DescentSetOfPartition[#] === d &]},
          polyDiffZero[Total[q^CrossingNumberPartition[#] & /@ parts],
            ACoeff[4 - Length[d], 2, q]]]] /@ {{}, {1}, {2}, {3}})]];
  okBell = And @@ Table[
    Sum[Total[CoefficientList[ACoeff[m, ell, q], q]], {ell, m}] === BellB[m - 1],
    {m, 2, 6}];
  audit["A11 Sum_l a(m,l) = Bell(m-1), m <= 6", okBell];
  log["== audit summary: " <> ToString[passed] <> " passed, " <>
     ToString[failed] <> " failed, " <> ToString[warned] <> " warnings =="];
  {lines, passed, failed, warned}
]

End[]
EndPackage[]
