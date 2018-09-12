(* ::Package:: *)

BeginPackage["sumrules`"]
Unprotect @@ Names["sumrules`*"];
ClearAll @@ Names["sumrules`*"];

opts0::usage="This procedures recursively determines the lower bound of the Borel parameter \[Tau] and the optimized continuum parameter \!\(\*SubscriptBox[\(s\), \(0\)]\) such that the two values are self-consistent within specified tolerances. The function takes three real-numbered arguements: a seed value for \!\(\*SubscriptBox[\(s\), \(0\)]\), a seed value for \!\(\*SubscriptBox[\(\[Tau]\), \(min\)]\), and an integer-valued sum rule weight, where 0 indicates the lowest-weighted sum rule."
lowerbound\[Tau]::usage="Defines the lower bound of the Borel parameter \[Tau] as defined by a 10% pole contribution."
upperbound\[Tau]::usage="Defines the upper bound of the Borel parameter \[Tau] as defined by requiring OPE convergence."
Begin["`Private`"]
Needs["sumruleExpressions`"]
(* THIS SHOULD BE REMOVED IN A FUTURE VERSION TO SOMETHING MORE ROBUST *)
applyNumericValues = {m -> ((Subscript[\[Alpha], s][\[Alpha], A, \[Nu], MZ]/Subscript[\[Alpha], s][\[Alpha], A, 2, MZ])^(1/A)) 0.003861525679074026`, M ->runningMass[4.18, \[Alpha], A, \[Nu],MZ], \[Nu] -> 4.18,vev[g^3 GGG] -> (8.2) vev[\[Alpha]GG],vev[\[Psi]\[Sigma]G\[Psi]] -> I  (1/M)(massratio(0.8)(vev[mq 
\!\(\*OverscriptBox[\(q\), \(_\)]\)])) (*Coefficient and ratio translates LPN VEV to our convention*), vev[\[Alpha]GG] -> 0.075,vev[mq 
\!\(\*OverscriptBox[\(q\), \(_\)]\)] -> -1/2 (0.1304/Sqrt[2.])^2 (0.13957018)^2, \[Alpha] -> 0.1184, A -> 23/12, \[Nu] -> 4.18, Mtau ->1.77686,MZ->91.188,vev[GG] -> vev[\[Alpha]GG]/\[Alpha], massratio -> 1460.7(*1371.25*)}
(* END OF THINGS TO BE REMOVED IN A FUTURE VERSION*)


(*The following are definitions on the upper and lower bound of the Borel parameter \[Tau]. *)
upperbound\[Tau][k_]:=
Module[{percent,cond1,cond2},
   percent=1/4;
   cond1=\[Tau]/.FindRoot[Abs[Re@dim4Asym[\[Tau],k]/Re@pertAsym[\[Tau],k]]==percent//.applyNumericValues,{\[Tau],0.5}];
   cond2=\[Tau]/.FindRoot[Abs[Re@(dim6Asym[\[Tau],k])/Re@dim4Asym[\[Tau],k]]==percent//.applyNumericValues,{\[Tau],0.5}];
   Piecewise[{{cond2,cond1>cond2},{cond1,cond1<cond2}}]
]
PC0[\[Tau]_,s0_,k_]:=(c10[\[Tau],s0,k]+c40[\[Tau],s0,k]+c60[\[Tau],s0,k]+c50[\[Tau],s0,k])/((pertAsym[\[Tau],k]+dim4Asym[\[Tau],k]+dim6Asym[\[Tau],k]+dim5Asym[\[Tau],k])//Chop) //.applyNumericValues//ReleaseHold//Chop
lowerbound\[Tau][k_?NumericQ,s0_?NumericQ]:=\[Tau]/.FindRoot[{PC0[\[Tau],s0,k]==.1},{\[Tau],0.1}]

opts0[s0init_?NumericQ,lowerbound\[Tau]init_?NumericQ,k_?IntegerQ]:=
   Module[{s0initial,s0opt1,s0opt2,lowerbound\[Tau]initial,lowerbound\[Tau]1,lowerbound\[Tau]2},
      (*Define initial values for initial continuum parameter s0 and seed value for the lower bound of the Borel parameter \[Tau]*)
      (* Below I'm just spacing out the initial values of s0 and the lower bound of \[Tau] for each of the two iterations of *)
      (* the recursive process. The values are not specific, and shouldn't matter a ton *)
      s0initial=s0init;
      s0opt1=s0initial+10;
      s0opt2=s0initial-10;
      lowerbound\[Tau]initial=lowerbound\[Tau]init;
      lowerbound\[Tau]1=lowerbound\[Tau]initial+0.1;
      lowerbound\[Tau]2=lowerbound\[Tau]initial+0.05;
      (* Determining the lower bound of the Borel parameter \[Tau] as well as performing a chi-squared minimization to *)
      (* determine the optimized value of the continuum parameter s0 are dependant on one another, so here I've tried *)
      (* to do both recursively until Mathematica arrives at a self-consistent value of both the lower bound on \[Tau] and the optimized value of s0 *)
      (* The loop continues until either the continuum parameter optimization or the lower bounds on \[Tau] reaches a precise enough agreement; there's no specific reason for these thresholds *)
      While[Abs[s0opt1-s0opt2]>0.05||Abs[lowerbound\[Tau]1-lowerbound\[Tau]2]>0.005,
         Module[{avgM,s0,\[Tau],root1,root2,tauInterval,compiled\[Chi]},
            (*sample 20 points in the tau range*)
            tauInterval=Abs[upperbound\[Tau][k]-lowerbound\[Tau]initial]/20;
            (* In order to perform the chi-sqared fitting, we compare against a "smeared" average of the s0 vs. hadronic mass curves. Problems could arise if there are extraneous non-intersecting curves added to the fit, and a median might be better suited here in a future implementation *)
            avgM=Sum[mx[\[Tau],s0,k]/Length[Table[\[Tau],{\[Tau],(lowerbound\[Tau]initial),(upperbound\[Tau][k]),tauInterval}]],{\[Tau],lowerbound\[Tau]initial,(upperbound\[Tau][k]),tauInterval}];
            compiled\[Chi]=Compile[{s0,\[Tau]min},Sum[(mx[\[Tau],s0,0]/avgM-1)^2,{\[Tau],\[Tau]min,(upperbound\[Tau][0]),tauInterval}]];
            s0opt1=s0/.FindMinimum[{compiled\[Chi][s0,lowerbound\[Tau]initial],s0>50},{s0,s0init},PrecisionGoal->0.001][[2]];
            root1=FindRoot[PC0[\[Tau],s0opt1,k]==.1,{\[Tau],0.01}];
            lowerbound\[Tau]1=\[Tau]/.root1//Chop;
            Echo[{s0opt1,lowerbound\[Tau]1}];
            s0opt2=s0/.FindMinimum[{compiled\[Chi][s0,lowerbound\[Tau]1],s0>50},{s0,s0init},PrecisionGoal->0.001][[2]];
            root2=FindRoot[PC0[\[Tau],s0opt2,k]==.1,{\[Tau],0.01}]//Chop;
            lowerbound\[Tau]2=\[Tau]/.root2;
            Echo[{s0opt2,lowerbound\[Tau]2}];
            s0initial=s0opt2;
            lowerbound\[Tau]initial=lowerbound\[Tau]2//Chop
        ];
      ];
   {s0opt1,s0opt2,lowerbound\[Tau]1,lowerbound\[Tau]2}
]


End[]
Protect @@ Names["sumrules`*"];
EndPackage[]
