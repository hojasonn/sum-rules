(* ::Package:: *)

BeginPackage["sumrules`"]
Unprotect @@ Names["sumrules`*"];
ClearAll @@ Names["sumrules`*"];

s0opt::usage""

Begin["`Private`"]
opts0[s0init_, lowerbound\[Tau]init_] :=
 
 Module[{s0initial, s0opt1, s0opt2, lowerbound\[Tau]initial, 
   lowerbound\[Tau]1, lowerbound\[Tau]2},
  s0initial = s0init;
  s0opt1 = s0initial + 10;
  s0opt2 = s0initial - 10;
  lowerbound\[Tau]initial = lowerbound\[Tau]init;
  lowerbound\[Tau]1 = lowerbound\[Tau]initial + 0.1;
  lowerbound\[Tau]2 = lowerbound\[Tau]initial + 0.05;
  While[Abs[s0opt1 - s0opt2] > 0.05 || 
    Abs[lowerbound\[Tau]1 - lowerbound\[Tau]2] > 0.005,
   Module[{avgM, s0, \[Tau], root1, root2, tauInterval},
     (*sample 20 points in the tau range*)
     
     tauInterval = Abs[upperbound\[Tau] - lowerbound\[Tau]initial]/20;
      avgM = Sum[
       mx[\[Tau], s0]/
        Length[Table[\[Tau], {\[Tau], (lowerbound\[Tau]initial), \
(upperbound\[Tau]), tauInterval}]], {\[Tau], 
        lowerbound\[Tau]initial, (upperbound\[Tau]), tauInterval}];
     s0opt1 = 
      s0 /. FindMinimum[{Sum[(mx[\[Tau], s0]/avgM - 1)^2, {\[Tau], 
            lowerbound\[Tau]initial, (upperbound\[Tau]), 
            tauInterval}]}, {s0, s0initial}, PrecisionGoal -> 0.001][[
        2]];
     root1 = 
      FindRoot[Evaluate@PC0[\[Tau], s0opt1] == .1, {\[Tau], 0.1}] // 
       Chop;
     lowerbound\[Tau]1 = \[Tau] /. root1;
     (*avgM=Sum[mx[\[Tau],s0]/Length[
     Table[\[Tau],{\[Tau],(lowerbound\[Tau]init),(upperbound\[Tau]),
     0.01}]],{\[Tau],lowerbound\[Tau]init,(upperbound\[Tau]),0.01}];*)

          s0opt2 = 
      s0 /. FindMinimum[{Sum[(mx[\[Tau], s0]/avgM - 1)^2, {\[Tau], 
            lowerbound\[Tau]1, (upperbound\[Tau]), 
            tauInterval}]}, {s0, s0initial}, PrecisionGoal -> 0.001][[
        2]];
     root2 = 
      FindRoot[Evaluate@PC0[\[Tau], s0opt2] == .1, {\[Tau], 0.1}] // 
       Chop;
     lowerbound\[Tau]2 = \[Tau] /. root2;
     s0initial = s0opt2;
     lowerbound\[Tau]initial = lowerbound\[Tau]2
     ];
   ];
  {s0opt1, s0opt2, lowerbound\[Tau]1, lowerbound\[Tau]2}
  ]

End[]
Protect @@ Names["sumrules`*"];
EndPackage[]
