(* Wolfram Language re-derivation of the (4,2) counterexample to Theorem 3.1
   of Bagno-Garber, arXiv:2606.13149. Run from the repository root:
     wolframscript -file prototype\wl_check_4_2.wls
   Requires Wolfram Mathematica (any recent version with wolframscript). *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "mantiq-project", "src", "SetPartitionQSym.wl"}]];
vars = {x1, x2, x3, x4};
q = Global`q;
Q42 = QSet[4, 2, q, vars];
Print["Q_4(Set(4,2)) = ", Q42];
(* Schur expansion in the hook basis used by the paper *)
corr = CorrectedRHS[4, 2, q, vars];
Print["CorrectedRHS(4,2) = ", corr];
Print["LHS == CorrectedRHS: ", Simplify[Expand[Q42 - corr]] === 0];
paper = PaperRHS[4, 2, q, vars];
Print["PaperRHS(4,2) = ", paper];
Print["LHS - PaperRHS = ", Simplify[Expand[Q42 - paper]]];
(* per-partition table *)
Do[
  Module[{rgs = p, des = DescentSetOfPartition[p], cr = CrossingNumberPartition[p],
     F = FundamentalQSym[4, DescentSetOfPartition[p], vars], blocks},
    blocks = BlocksOfPartition[p];
    Print["pi = ", blocks, "  Des = ", des, "  cr = ", cr]],
  {p, Select[SetPartitionsRGS[4], Max[#] + 1 == 2 &]}];
(* coefficient of s_(4) on the LHS: the corrected theorem says it is A(4,2,q) *)
Print["A(4,2,q) = ", ACoeff[4, 2, q]];
Print["T_(3,2)(q) = ", TouchardRiordanT[3, 2, q]];
