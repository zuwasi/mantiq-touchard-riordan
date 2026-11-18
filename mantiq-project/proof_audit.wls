#!/usr/bin/env wolframscript
(* proof_audit.wls -- standalone run of the paper-specific proof audit.
   Prints the audit lines and exports exports/audit_report.txt. *)

projectDir = DirectoryName[$InputFileName];
SetDirectory[projectDir];
Get[FileNameJoin[{projectDir, "src", "SetPartitionQSym.wl"}]];
Print["package loaded, ", $Version];

{lines, passed, failed, warned} = RunProofAudit[];
Scan[Print, lines];

Export[FileNameJoin[{projectDir, "exports", "audit_report.txt"}],
  StringJoin[Riffle[lines, "\n"]], "Text"];

If[failed > 0, Exit[1], Print["proof audit OK"]];
