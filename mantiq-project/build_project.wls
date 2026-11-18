#!/usr/bin/env wolframscript
(* build_project.wls -- MantiQ Factor build for the Bagno-Garber project.
   Loads the package, runs the validation suite, exports figures, tables and
   reports into exports/. *)

projectDir = DirectoryName[$InputFileName];
SetDirectory[projectDir];

Do[If[! DirectoryQ[FileNameJoin[{projectDir, d}]],
   CreateDirectory[FileNameJoin[{projectDir, d}]]],
  {d, {"exports", "docs"}}];

Get[FileNameJoin[{projectDir, "src", "SetPartitionQSym.wl"}]];
Print["package loaded, ", $Version];

q = Global`q;
t0 = AbsoluteTime[];

(* ---------- 1. validation suite ---------- *)
suite = RunValidationSuite[6];
Export[FileNameJoin[{projectDir, "exports", "validation_results.json"}],
  suite, "JSON"];
Print["validation suite done"];

(* ---------- 2. audit-style checks beyond the suite ---------- *)
{lines, passed, failed, warned} = RunProofAudit[];
Scan[Print, lines];

Export[FileNameJoin[{projectDir, "exports", "audit_report.txt"}],
  StringJoin[Riffle[lines, "\n"]], "Text"];

(* ---------- 3. figures and data exports ---------- *)

(* figure 1: corrected vs paper coefficients for (n,l) grid, n <= 7 *)
data = Flatten[Table[
    With[{a = ACoeff[n - k, ell, q], t = TouchardRiordanT[n - 1 - k, n - k - ell, q]},
      <|"n" -> n, "ell" -> ell, "k" -> k,
        "A_string" -> ToString[a, InputForm],
        "T_string" -> ToString[t, InputForm],
        "agree" -> Simplify[a - t] === 0|>],
    {n, 2, 7}, {ell, n}, {k, 0, n - 1}], 2];
Export[FileNameJoin[{projectDir, "exports", "coefficient_census.csv"}], data];

(* figure: first counterexample bar chart, n = 4, ell = 2 *)
chart = BarChart[
  {Total@CoefficientRules[ACoeff[4, 2, q], q][[All, 2]],
   Total@CoefficientRules[TouchardRiordanT[3, 2, q], q][[All, 2]]},
  ChartLabels -> {"corrected A(4,2,q) = q", "paper T_(3,2)(q) = 0"},
  ChartStyle -> {ColorData[97, 1], ColorData[97, 4]},
  PlotLabel -> "Coefficient of s_(4) in Q_4(Set(4,2))",
  ImageSize -> 480];
Export[FileNameJoin[{projectDir, "exports", "figure_counterexample_4_2.png"}], chart];

(* figure: A(m,l,q) table as grid for m <= 7 *)
atable = Table[
  If[ell > m, "", ToString[ACoeff[m, ell, q], InputForm]],
  {m, 1, 7}, {ell, 1, 7}];
Export[FileNameJoin[{projectDir, "exports", "A_coefficients_table.csv"}], atable];

(* figure: Narayana / Motzkin distribution plot for n = 8 *)
narPlot = ListLinePlot[
  Table[Binomial[8, ell] Binomial[8, ell - 1]/8, {ell, 1, 8}],
  PlotLabel -> "Narayana(8,ell): total NC adjacency-free counts",
  Mesh -> All, ImageSize -> 480];
Export[FileNameJoin[{projectDir, "exports", "figure_narayana_8.png"}], narPlot];

(* arc-diagram gallery of Set(4,2) with Des and cr annotated *)
gallery = GraphicsGrid[
  Partition[
    With[{parts = Select[SetPartitionsRGS[4], Max[#] + 1 == 2 &]},
      Labeled[PartitionArcDiagram[#],
        "Des=" <> ToString[DescentSetOfPartition[#], InputForm] <>
          ", cr=" <> ToString[CrossingNumberPartition[#]]] & /@ parts],
    UpTo[4]],
  ImageSize -> 720];
Export[FileNameJoin[{projectDir, "exports", "figure_set42_gallery.png"}], gallery];

(* ---------- 4. validation report (human readable) ---------- *)
reportLines = Map[
  "* " <> ToString[#] <> ": " <>
     If[TrueQ[suite[#]["pass"]], "PASS", "FAIL -- " <>
        ToString[suite[#]["failures"], InputForm]] &,
  Select[Keys[suite], StringQ[#] && AssociationQ[suite[#]] &]];
reportLines = Select[reportLines, StringQ];
report2 = StringJoin[
  "# Validation report -- MantiQ Factor reproduction of Bagno-Garber (arXiv:2606.13149)\n\n",
  "Wolfram version: ", $Version, "\n\n",
  StringJoin[Riffle[reportLines, "\n"]], "\n\n",
  "## Headline results\n\n",
  "- Theorem 3.1 of the paper is false from (n,l) = (4,2) on: the coefficient\n",
  "  of s_(n-k,1^k) is the crossing polynomial A(n-k,l,q) of adjacency-free\n",
  "  partitions, not the Touchard-Riordan value T_(n-1-k,n-k-l)(q).\n",
  "- Corollary 3.3 (3)-(4) (noncrossing) is correct: the coefficients are the\n",
  "  Motzkin triangle numbers NCMatchCount(m-1, m-l), and the block-count\n",
  "  totals are Narayana(n,l) and Motzkin numbers.\n",
  "- New identities verified: collapse lemma c_D = A(n-|D|,l,q) for every\n",
  "  descent set D; inclusion-exclusion for A(m,l,q); Sum_l a(m,l) = Bell(m-1).\n"];
Export[FileNameJoin[{projectDir, "exports", "validation_report.md"}], report2, "Text"];

(* ---------- 5. assert exports exist ---------- *)
required = {"audit_report.txt", "validation_report.md", "validation_results.json",
   "coefficient_census.csv", "figure_counterexample_4_2.png",
   "A_coefficients_table.csv", "figure_narayana_8.png", "figure_set42_gallery.png"};
missing = Select[required, FileType[FileNameJoin[{projectDir, "exports", #}]] =!= File &];
If[missing =!= {},
  Print["MISSING EXPORTS: ", missing]; Exit[1],
  Print["build OK: all exports present in exports/"]];
Print["elapsed seconds: ", N[AbsoluteTime[] - t0]];
