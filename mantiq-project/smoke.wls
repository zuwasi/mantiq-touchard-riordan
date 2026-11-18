(* smoke test for SetPartitionQSym *)
SetDirectory[$InputFileName // DirectoryName];
Get[FileNameJoin[{Directory[], "src", "SetPartitionQSym.wl"}]];

Print["version: ", $Version];

(* basic counts *)
Print["partitions of 4: ", Length[SetPartitionsRGS[4]]];  (* expect 15 *)
Print["Des of {{1,3},{2,4}} rgs {0,1,0,1}: ", DescentSetOfPartition[{0, 1, 0, 1}]];
Print["cr of {0,1,0,1}: ", CrossingNumberPartition[{0, 1, 0, 1}]];  (* expect 1 *)
Print["NoncrossingQ {0,1,0,1}: ", NoncrossingQ[{0, 1, 0, 1}]];  (* False *)
Print["NoncrossingQ {0,1,1,0}: ", NoncrossingQ[{0, 1, 1, 0}]];  (* True *)

(* F test: F[2, {}, {x1,x2}] = M_(2)+M_(1,1) part -> x1^2 + x1 x2 + x2^2 ... *)
Print["F[2,{},]: ", FundamentalQSym[2, {}, {x1, x2}]];
Print["F[2,{1},]: ", FundamentalQSym[2, {1}, {x1, x2}]];

(* hook Schur: s_(2,1) in 3 vars *)
Print["HookSchur[2,1]: ", HookSchur[2, 1, {a, b, c}]];

(* ACoeff table *)
Print["A(4,2,q) = ", ACoeff[4, 2, q]];  (* expect q *)
Print["A(5,3,q) = ", ACoeff[5, 3, q]];  (* expect 2+5q *)
Print["T(3,2,q) = ", TouchardRiordanT[3, 2, q]];  (* expect 0 *)
Print["T(4,2,q) = ", TouchardRiordanT[4, 2, q]];  (* expect 2+q *)

(* counterexample: n=4, ell=2, k=0 *)
vars = Table[x[i], {i, 4}];
Print["QSet(4,2) == CorrectedRHS(4,2): ",
  Simplify[Expand[QSet[4, 2, q, vars] - CorrectedRHS[4, 2, q, vars]]] === 0];
Print["QSet(4,2) == PaperRHS(4,2): ",
  Simplify[Expand[QSet[4, 2, q, vars] - PaperRHS[4, 2, q, vars]]] === 0];

(* SYT descent sets for hook (3,1) *)
Print["SYT Des sets of (3,1): ", Sort[SYTDescentSets[{3, 1}]]];
Print["NCMatchCount[3,1] = ", NCMatchCount[3, 1], "  (expect 3)"];
Print["MotzkinNumber[5] = ", MotzkinNumber[5], "  (expect 21)"];

(* full suite, small *)
res = RunValidationSuite[5];
Print["suite: ", res]
