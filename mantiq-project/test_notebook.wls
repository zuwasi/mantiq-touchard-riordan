#!/usr/bin/env wolframscript
(* test_notebook.wls -- validates the generated teaching notebook:
   1. the file parses as a Notebook expression,
   2. no cell contains BoxData (plain-string cells only),
   3. every Input cell parses,
   4. all Input cells evaluate in order without error. *)

projectDir = DirectoryName[$InputFileName];
SetDirectory[projectDir];
nbPath = FileNameJoin[{projectDir, "TouchardRiordanNotebook.nb"}];
nb = Get[nbPath];

(* load the package up front so unqualified names in later cells resolve to
   SetPartitionQSym` symbols during the parse loop, exactly as they would in
   the frontend after the load cell has run *)
Get[FileNameJoin[{projectDir, "src", "SetPartitionQSym.wl"}]];

If[Head[nb] =!= Notebook, Print["FAIL: not a Notebook expression"]; Exit[1]];
boxData = Cases[nb, BoxData[___], Infinity];
If[boxData =!= {}, Print["FAIL: found BoxData cells"]; Exit[1]];
Print["no BoxData: OK"];

inputCells = Cases[nb, Cell[code_String, "Input", ___] :> code, Infinity];
Print["input cells: ", Length[inputCells]];
Do[
  With[{c = inputCells[[i]]},
   If[ToExpression[c, InputForm, Hold] === $Failed,
     Print["PARSE FAIL in cell ", i, ": ", StringTake[c, UpTo[80]]];
     Exit[1]]],
  {i, Length[inputCells]}];
Print["all input cells parse: OK"];

(* evaluate in order, with NotebookDirectory[] replaced by the project dir *)
ok = True;
Do[
  With[{c = inputCells[[i]],
    code = StringReplace[inputCells[[i]], "SetDirectory[NotebookDirectory[]]" ->
      "SetDirectory[" <> ToString[projectDir, InputForm] <> "]"]},
    res = Check[ToExpression[code, InputForm]; Null, $Failed];
    If[res === $Failed || Head[res] === $Failed,
      ok = False; Print["EVAL FAIL in cell ", i, ": ",
        StringTake[c, UpTo[80]]]]],
  {i, Length[inputCells]}];
If[ok, Print["all input cells evaluate: OK"], Exit[1]];
Print["notebook validation PASSED"];
