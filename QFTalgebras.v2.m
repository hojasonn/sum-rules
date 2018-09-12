(* ::Package:: *)

BeginPackage["QFTAlgebrasv2`"]
Unprotect @@ Names["QFTAlgebrasv2`*"];
ClearAll @@ Names["QFTAlgebrasv2`*"];

\[Gamma]::usage = "Dirac Matrix. \[Gamma][\[Mu]], \[Gamma][5], \[Gamma][id] all valid options."
Begin["`Private`"]
(* Delta Functions - Dirac Indices*)
(* Delta Functions - Color Indices*)

(*Gamma Matrices*)
\[Gamma]/:CircleTimes[\[Gamma][\[Mu]_],\[Gamma][\[Mu]_]]/;!SameQ[\[Mu],5]:=d
\[Gamma]/:\[Gamma][\[Mu]_,i_,j_]\[Gamma][\[Mu]_,j_,k_]/;!SameQ[i,k]&&!SameQ[\[Mu],5]:=d \[Delta][i,k,"dirac"]
sigma[\[Mu]_,\[Nu]_]:=(I/2)*(CircleTimes[\[Gamma][\[Mu]],\[Gamma][\[Nu]]]-CircleTimes[\[Gamma][\[Nu]],\[Gamma][\[Mu]]])
sigma[\[Mu]_,\[Nu]_,i_,j_]:=Module[{int=Unique[k]},(I/2)*(\[Gamma][\[Mu],i,int]\[Gamma][\[Nu],int,j]-\[Gamma][\[Nu],i,int]\[Gamma][\[Mu],int,j])]

(*Metrics*)
g/:g[\[Mu]_,\[Mu]_]:=d
g/:g[\[Mu]_,\[Lambda]_]g[\[Lambda]_,\[Nu]_]:=g[\[Mu],\[Nu]]
g/:g[\[Mu]_,\[Lambda]_]g[\[Nu]_,\[Lambda]_]:=g[\[Mu],\[Nu]]
g/:g[\[Lambda]_,\[Mu]_]g[\[Lambda]_,\[Nu]_]:=g[\[Mu],\[Nu]]
g/:g[\[Mu]_,\[Nu]_] vec[p_,\[Nu]_] :=vec[p,\[Mu]] 
g/:g[\[Mu]_,\[Nu]_] vec[p_,\[Mu]_]:=vec[p,\[Nu]] 
g/:g[\[Mu]_,\[Nu]_] vec[p_,\[Nu]_] vec[k_,\[Mu]_]:=vec[p,\[Nu]] vec[k,\[Nu]]
g/:g[\[Mu]_,\[Rho]_] vec[p_,\[Nu]_] vec[k_,\[Rho]_]:=vec[p,\[Nu]] vec[k,\[Mu]]
g/:g[\[Mu]_,\[Nu]_]g[\[Mu]_,\[Nu]_]:=g[\[Mu],\[Mu]]
g/:g[\[Mu]_,\[Nu]_]^2:=g[\[Mu],\[Mu]]
g/:g[\[Mu]_,\[Nu]_]CircleTimes[\[Gamma][\[Nu]_],A__]/;!SameQ[\[Nu],5]:=CircleTimes[\[Gamma][\[Mu]],A]
End[]
Protect @@ Names["QFTAlgebrasv2`*"];
EndPackage[]



