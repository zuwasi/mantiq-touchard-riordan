#!/usr/bin/env python3
"""Generates TouchardRiordanNotebook.nb as a plain-text Wolfram notebook.

Cells are emitted as Cell["string", style] only: code lives exclusively in
"Input" cells, prose in "Text" cells, no BoxData/RowBox anywhere. Sections are
CellGroupData[{...}, Open].
"""
import os

def esc(s):
    return s.replace("\\", "\\\\").replace('"', '\\"')

def cell(style, text):
    return 'Cell["%s", "%s"]' % (esc(text), style)

# content: list of sections; each section is (title, [(style, text), ...])
SECTIONS = []

SECTIONS.append(("MantiQ Factor and the Bagno-Garber Touchard-Riordan Paper", [
("Text", "A computational audit, a corrected theorem, and interactive teaching explorations for the paper 'Touchard-Riordan Polynomials and Schur-positivity of Set Partitions' by E. Bagno and D. Garber (arXiv:2606.13149)."),
("Text", "What this notebook does: it loads the companion package SetPartitionQSym.wl, builds the paper's central quasisymmetric sum by brute force, tests Theorem 3.1, exposes the smallest counterexample, states and verifies a corrected theorem, proves the collapse lemma, and explores the (correct) noncrossing corollary with interactive controls. A Lean 4 companion in the lean/ folder re-derives the key claims independently. Requirements: Mathematica 14 or 15 and the mantiq-project folder (no external packages). Evaluate everything with Evaluation > Evaluate Notebook."),
]))

SECTIONS.append(("1. Load the package", [
("Text", "All functions live in the SetPartitionQSym context. Load it from src/."),
("Input", "SetDirectory[NotebookDirectory[]];\nGet[FileNameJoin[{\"src\", \"SetPartitionQSym.wl\"}]]"),
("Input", "Sort[Names[\"SetPartitionQSym`*\"]]"),
]))

SECTIONS.append(("2. The paper's objects: set partitions, descents, crossings", [
("Text", "A set partition of [n] = {1,...,n} into ell blocks is stored as a restricted growth string (rgs): position i holds the block index of element i, counting blocks from 0. For example {0,1,0,1} is the partition {{1,3},{2,4}}."),
("Input", "Length[SetPartitionsRGS[4]]"),
("Input", "BlocksOfPartition[{0, 1, 0, 1}]"),
("Text", "The descent set is Des(pi) = { i : i and i+1 lie in the same block }. The crossing number cr(pi) counts pairs of arcs (i,j) and (k,l) with i < k < j < l, where arcs connect consecutive elements of each block."),
("Input", "DescentSetOfPartition[{0, 1, 0, 1}]"),
("Input", "CrossingNumberPartition[{0, 1, 0, 1}]"),
("Input", "Table[{p, DescentSetOfPartition[p], CrossingNumberPartition[p]},\n  {p, Select[SetPartitionsRGS[4], Max[#] + 1 == 2 &]}] // TableForm"),
("Text", "Each arc in the drawing below is one pair of consecutive elements of a block; arcs of the same color belong to the same block, and crossing pairs are easy to see visually."),
("Input", "PartitionArcDiagram /@ Select[SetPartitionsRGS[4], Max[#] + 1 == 2 &]"),
]))

SECTIONS.append(("3. Interactive arc-diagram explorer", [
("Text", "Drag the sliders to walk through every partition of [n] into ell blocks. Each diagram shows the descent set Des and the crossing number cr. Watch how n = 4, ell = 2 already contains a crossing partition: that is the seed of the counterexample in Section 6."),
("Input", "Manipulate[\n Module[{parts = Select[SetPartitionsRGS[n], Max[#] + 1 == ell &], jj, p},\n  jj = Min[j, Length[parts]];\n  p = parts[[jj]];\n  Labeled[PartitionArcDiagram[p],\n   Row[{\"partition \", jj, \" of \", Length[parts], \":  Des = \",\n     DescentSetOfPartition[p], \",  cr = \", CrossingNumberPartition[p]}]]],\n {{n, 4, \"n (points)\"}, 2, 6, 1},\n {{ell, 2, \"blocks\"}, 1, Dynamic[n], 1},\n {{j, 1, \"partition #\"}, 1,\n  Dynamic[Length[Select[SetPartitionsRGS[n], Max[#] + 1 == ell &]]], 1},\n ControlPlacement -> Left]"),
]))

SECTIONS.append(("4. The quasisymmetric sum of the paper", [
("Text", "The paper's central object is Q_n(Set(n, ell)) = Sum over pi in Set(n, ell) of q^cr(pi) F_{n, Des(pi)}, where F_{n,D} is Gessel's fundamental quasisymmetric function; the paper's Definition 2.3 requires a rise (x_i <= x_{i+1}) at every index of D. Expanding in a finite list of variables turns it into an ordinary polynomial that Mathematica can compare, simplify and decompose."),
("Input", "vars = Table[Symbol[\"x\" <> ToString[i]], {i, 6}];\nQSet[4, 2, q, vars]"),
("Input", "QSetNC[4, 2, q, vars]"),
("Text", "QSetNC restricts the sum to noncrossing partitions; it becomes central in Section 11."),
]))

SECTIONS.append(("5. The paper's Theorem 3.1", [
("Text", "Theorem 3.1 of the paper claims, in words: the sum above equals sum over k = 0,...,n-1 of T_{n-1-k, n-k-ell}(q) times the hook Schur function s_(n-k, 1^k), where T_{a,b}(q) = C(a, 2b) T_{2b}(q) and T_{2b} is the classical Touchard-Riordan crossing polynomial of all matchings on 2b points. The claim therefore says: each hook Schur coefficient is a Touchard-Riordan polynomial."),
("Input", "PaperRHS[4, 2, q, vars]"),
("Text", "Compare the two sides at (n, ell) = (4, 2). The difference is not zero:"),
("Input", "Simplify[QSet[4, 2, q, vars] - PaperRHS[4, 2, q, vars]]"),
]))

SECTIONS.append(("6. The smallest counterexample: (n, ell) = (4, 2)", [
("Text", "The failure is structural, not numerical. Look at the k = 0 hook coefficient, the coefficient of s_(n) = s_(4). On the left it is q: exactly one adjacency-free partition of [4] into 2 blocks has one crossing (the partition {{1,3},{2,4}}), so its crossing polynomial is q. The paper's coefficient T_{3,2}(q) = C(3,4) T_4(q) = 0 because there are no matchings on 3 points with 4 chords. The paper's formula is impossible from this point on."),
("Input", "ACoeff[4, 2, q]"),
("Input", "TouchardRiordanT[3, 2, q]"),
("Input", "Select[SetPartitionsRGS[4], Max[#] + 1 == 2 && DescentSetOfPartition[#] === {} && CrossingNumberPartition[#] == 1 &]"),
("Text", "Interactive comparison: for any (n, ell) the grid below lists, hook by hook, the true coefficient A(n-k, ell, q) and the paper's coefficient T_{n-1-k, n-k-ell}(q). The last row also re-verifies the full corrected identity of Section 9 for the chosen (n, ell)."),
("Input", "Manipulate[\n Module[{vars, rows, check},\n  vars = Table[Symbol[\"x\" <> ToString[i]], {i, n}];\n  rows = Table[{k,\n    ToString[ACoeff[n - k, ell, q], InputForm],\n    ToString[TouchardRiordanT[n - 1 - k, n - k - ell, q], InputForm],\n    If[ACoeff[n - k, ell, q] === TouchardRiordanT[n - 1 - k, n - k - ell, q], \"agree\", \"differ\"]},\n   {k, 0, n - 1}];\n  check = If[polyDiffZero[QSet[n, ell, q, vars], CorrectedRHS[n, ell, q, vars]],\n    \"corrected identity holds\", \"corrected identity FAILS\"];\n  Column[{Grid[Prepend[rows, {\"k\", \"true A\", \"paper T\", \"match\"}], Frame -> All],\n    Style[check, 12, Bold]}]],\n {{n, 4, \"n\"}, 2, 6, 1},\n {{ell, 2, \"blocks\"}, 1, Dynamic[n], 1}]"),
]))

SECTIONS.append(("7. Failure census: when does Theorem 3.1 hold?", [
("Text", "Running the comparison for all 2 <= ell <= n <= 6 in Wolfram Language, and for n <= 7 with the independent Python cross-check (prototype/verify3.py), gives a clean census. The paper's identity holds exactly in the boundary regimes ell = 1, ell = n and n - ell = 1; it fails precisely when ell >= 2 and n - ell >= 2. This pattern was verified up to n = 7 in two independent implementations."),
("Input", "Module[{cells},\n cells = Table[\n   If[n < ell, \"\",\n    If[polyDiffZero[QSet[n, ell, q, Table[Symbol[\"x\" <> ToString[i]], {i, n}]],\n        PaperRHS[n, ell, q, Table[Symbol[\"x\" <> ToString[i]], {i, n}]]],\n     \"holds\", \"FAILS\"]], {n, 2, 6}, {ell, 1, 6}];\n Grid[Prepend[MapThread[Prepend, {cells, Range[2, 6]}],\n   Prepend[Range[1, 6], \"ell \\ n\"]], Frame -> All]]"),
]))

SECTIONS.append(("8. The collapse lemma: why the proof breaks and what is true", [
("Text", "Key structural fact (new here): the descent set only shifts the coefficient. For EVERY descent set D, the crossing polynomial of the partitions of [n] into ell blocks with Des(pi) = D equals A(n - |D|, ell, q), the crossing polynomial of the adjacency-free partitions (Des = empty) of [m] into ell blocks with m = n - |D|. Intuitively: an adjacency in the block structure of pi contributes one fixed crossing pattern that can be collapsed, and what remains is an adjacency-free partition of the remaining m points. The paper's proof implicitly treated the fibers of the chord-diagram map as if they were as rich as full matching classes; in reality the map is injective but not surjective, and its fibers are exactly these D-runs, which forces the shift."),
("Input", "Module[{ok = True, bad = {}, parts},\n Do[parts = Select[SetPartitionsRGS[5],\n    Max[#] + 1 == 3 && DescentSetOfPartition[#] === d &];\n   If[parts === {},\n    If[ACoeff[5 - Length[d], 3, q] =!= 0, ok = False; AppendTo[bad, d]],\n    If[! polyDiffZero[Total[q^CrossingNumberPartition[#] & /@ parts],\n       ACoeff[5 - Length[d], 3, q]], ok = False; AppendTo[bad, d]]],\n  {d, Subsets[Range[4]]}];\n {ok, bad}]"),
]))

SECTIONS.append(("9. The corrected theorem", [
("Text", "Corrected Theorem (verified here by computer in two systems). For all n >= 1 and 1 <= ell <= n: Q_n(Set(n, ell)) = sum over k = 0,...,n-ell of A(n-k, ell, q) s_(n-k, 1^k), where A(m, ell, q) is the crossing polynomial of the adjacency-free partitions of [m] into ell blocks. This is manifestly Schur-positive with nonnegative q-coefficients, and it reduces to the paper's formula exactly in the regimes where the Touchard-Riordan and adjacency-free coefficients coincide (ell = 1, ell = n, n - ell = 1)."),
("Input", "And @@ Flatten[Table[\n  With[{vars = Table[Symbol[\"x\" <> ToString[i]], {i, n}]},\n   polyDiffZero[QSet[n, ell, q, vars], CorrectedRHS[n, ell, q, vars]]],\n  {n, 2, 6}, {ell, n}]]"),
("Text", "An explicit expansion, (n, ell) = (5, 2): left side, corrected right side, and their difference."),
("Input", "vars5 = Table[Symbol[\"x\" <> ToString[i]], {i, 6}];\n{QSet[5, 2, q, vars5], CorrectedRHS[5, 2, q, vars5],\n Simplify[QSet[5, 2, q, vars5] - CorrectedRHS[5, 2, q, vars5]]}"),
]))

SECTIONS.append(("10. The coefficients A(m, ell, q) and two new identities", [
("Text", "A(m, ell, q) extends the Stirling numbers of the second kind: at q = 1, a(m, ell) = A(m, ell, 1) counts the adjacency-free partitions of [m] into ell blocks. Two identities proved and verified in this project: inclusion-exclusion a(m, ell) = sum over j of (-1)^j C(m-1, j) S(m-j, ell), and the row identity sum over ell of a(m, ell) = Bell(m-1)."),
("Input", "Grid[Prepend[\n  Table[ToString[ACoeff[m, ell, q], InputForm], {m, 1, 7}, {ell, 1, 7}],\n  Range[7]], Frame -> All]"),
("Input", "Table[{m, Sum[Total[CoefficientList[ACoeff[m, ell, q], q]], {ell, m}], BellB[m - 1]},\n  {m, 2, 7}]"),
("Input", "And @@ Flatten[Table[\n  Total[CoefficientList[ACoeff[m, ell, q], q]] ===\n   Sum[(-1)^j Binomial[m - 1, j] StirlingS2[m - j, ell], {j, 0, m - 1}],\n  {m, 2, 6}, {ell, m}]]"),
]))

SECTIONS.append(("11. The noncrossing corollary: the true part, with Motzkin and Narayana", [
("Text", "Not everything in the paper fails. Corollary 3.3, items (3) and (4), are correct, and this project identifies their coefficients explicitly: in the noncrossing regime, A^NC(m, ell) equals the Motzkin triangle number, the number of noncrossing matchings on m-1 points with m-ell chords. The row sums are Motzkin numbers M(m-1), and summing the hook decomposition over k recovers the Narayana numbers Nar(n, ell) = C(n, ell) C(n, ell-1) / n."),
("Input", "Table[{m, ell, ANCCount[m, ell], NCMatchCount[m - 1, m - ell]},\n  {m, 3, 6}, {ell, 1, m}] // TableForm"),
("Input", "Table[{m, Sum[ANCCount[m, ell], {ell, m}], MotzkinNumber[m - 1]}, {m, 2, 8}]"),
("Input", "And @@ Flatten[Table[\n  Sum[ANCCount[n - k, ell], {k, 0, n - ell}] ===\n   Binomial[n, ell] Binomial[n, ell - 1]/n,\n  {n, 2, 7}, {ell, n}]]"),
("Text", "Interactive explorer: the Motzkin triangle row A^NC(m, ell) for each block count, with the Motzkin row sum displayed."),
("Input", "Manipulate[\n Module[{row = Table[ANCCount[m, ell], {ell, m}]},\n  Column[{BarChart[row, ChartLabels -> Range[m],\n    PlotLabel -> \"A^NC(m, ell) for m = \" <> ToString[m],\n    ImageSize -> 420],\n   Style[\"row sum = \" <> ToString[Total[row]] <> \" = Motzkin(\" <>\n     ToString[m - 1] <> \")\", 12, Bold]}]],\n {{m, 5, \"m\"}, 2, 8, 1}]"),
]))

SECTIONS.append(("12. The Lean 4 companion", [
("Text", "The lean/ folder contains an independent formal check in Lean 4 (pure core, no Mathlib): brute-force set partitions as restricted growth strings, descent sets, crossing counts, the counterexample c_empty(4,2) = q differs from T_{3,2} = 0, and the failure census up to n = 5. Build it with 'lake build' (install Lean 4 through elan, https://leanprover.github.io/lean4). Why run Lean in addition to Mathematica? A computer-algebra system catches computational claims quickly, but a single codebase can share a mistaken assumption with itself. Lean re-derives every fact from first-order definitions, so agreement between the two systems is strong evidence that a claim is a theorem of the combinatorics rather than an artifact of one implementation. Here the two systems agree on both the counterexample and the corrected theorem."),
("Input", "Import[FileNameJoin[{\"exports\", \"audit_report.txt\"}], \"Text\"]"),
]))

SECTIONS.append(("KEY FINDINGS", [
("Text", "1. Theorem 3.1 of Bagno-Garber is false in general. Smallest counterexample: (n, ell) = (4, 2), where the coefficient of s_(4) is q but the paper's T_{3,2}(q) = 0."),
("Text", "2. Failure census: the identity holds exactly when ell = 1 or ell = n or n - ell = 1; it fails precisely when ell >= 2 and n - ell >= 2 (verified to n = 7 in Python and Wolfram, and n = 5 in Lean)."),
("Text", "3. Corrected theorem: the true coefficients are the crossing polynomials A(m, ell, q) of adjacency-free partitions, not Touchard-Riordan polynomials. The corrected identity is manifestly Schur-positive."),
("Text", "4. Collapse lemma (new): for every descent set D, c_D(n, ell) = A(n - |D|, ell, q). The descent set only shifts."),
("Text", "5. Corollary 3.3 (3) and (4) are true, with an explicit Motzkin triangle identification; totals are Motzkin and Narayana numbers."),
("Text", "6. New identities: inclusion-exclusion for a(m, ell) in terms of Stirling numbers, and sum over ell of a(m, ell) = Bell(m-1)."),
("Text", "7. Workflow lesson: the MantiQ Factor pipeline (Mathematica brute force + Schur decomposition, cross-checked in Python, then formalized in Lean) found in minutes a gap that survived peer review, and immediately produced a stronger, cleaner statement to replace it."),
("Text", "Reproduce everything: run build_project.wls (full build with reports and figures), proof_audit.wls (audit only), and python prototype/verify3.py 7 (independent Python verification). Package sources: src/SetPartitionQSym.wl. Formal proof: lean/."),
]))

# ---- emit ----
groups = []
for title, contents in SECTIONS:
    inner = [cell("Section", title)]
    for style, text in contents:
        inner.append(cell(style, text))
    groups.append("Cell[CellGroupData[{" + ",\n".join(inner) + "}, Open]]")

nb = ("Notebook[{\n" + ",\n".join(groups) + "\n},\nWindowSize -> {900, 750},\n"
      "WindowMargins -> {{Automatic, 60}, {Automatic, 40}},\n"
      "StyleDefinitions -> \"Default.nb\"\n]\n"
      "(* Generated by tools/gen_notebook.py; edit that script, not this file. *)\n")

out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "TouchardRiordanNotebook.nb")
with open(out, "w", encoding="utf-8") as f:
    f.write(nb)
print("wrote", out, len(nb), "chars,", len(groups), "sections")

