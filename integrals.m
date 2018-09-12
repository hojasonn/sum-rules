(* ::Package:: *)

BeginPackage["integrals`"]
Unprotect @@ Names["integrals`*"];
ClearAll @@ Names["integrals`*"];

integralDef::usage = "Massless integral definitions for rank-0 to rank-3 type numerator structure."
d::usage = "Number of dimensions in dim-reg scheme. Assume d = 4 + 2\[Epsilon]"
\[Epsilon]::usage = "Expansion parameter in dim-reg scheme."
\[Nu]::usage = "Renormalization scale."
Begin["`Private`"]
integralDef[p_,r_,s_]:=(I \[Nu]^(2\[Epsilon]))/(4\[Pi])^2 (-p.p/(4\[Pi] \[Nu]^2))^\[Epsilon] (2\[Pi])^d/(p.p)^(r+s-2) (Gamma[2-r+\[Epsilon]]Gamma[2-s+\[Epsilon]]Gamma[r+s-2-\[Epsilon]])/(Gamma[r]Gamma[s]Gamma[4-r-s+2\[Epsilon]])
integralDef[p_,r_,s_,\[Mu]_]:=(I \[Nu]^(2\[Epsilon]))/(4\[Pi])^2 (-p.p/(4\[Pi] \[Nu]^2))^\[Epsilon] (2\[Pi])^d/(p.p)^(r+s-2) momentum[p,\[Mu]] (Gamma[3-r+\[Epsilon]]Gamma[2-s+\[Epsilon]]Gamma[r+s-2-\[Epsilon]])/(Gamma[r]Gamma[s]Gamma[5-r-s+2\[Epsilon]])
integralDef[p_,r_,s_,\[Mu]_,\[Nu]1_]:=(I \[Nu]^(2\[Epsilon]))/(4\[Pi])^2 (-p.p/(4\[Pi] \[Nu]^2))^\[Epsilon] (2\[Pi])^d/(p.p)^(r+s-2) (g[\[Mu],\[Nu]1]p.p (Gamma[3-r+\[Epsilon]]Gamma[3-s+\[Epsilon]]Gamma[r+s-3-\[Epsilon]])/(2 Gamma[r]Gamma[s]Gamma[6-r-s+2\[Epsilon]])+momentum[p,\[Mu]]momentum[p,\[Nu]1] (Gamma[4-r+\[Epsilon]]Gamma[2-s+\[Epsilon]]Gamma[r+s-2-\[Epsilon]])/(Gamma[r]Gamma[s]Gamma[6-r-s+2\[Epsilon]]))
integralDef[p_,r_,s_,\[Mu]_,\[Nu]1_,\[Sigma]_]:=
((\[Pi]^(2+\[Epsilon]) I^(1-4-2\[Epsilon]))(p.p)^(2+\[Epsilon]-r-s))/(Gamma[r]Gamma[s]Gamma[4+2\[Epsilon]-r-s+3]) (momentum[p,\[Mu]]momentum[p,\[Nu]1]momentum[p,\[Sigma]] Gamma[2+\[Epsilon]-s]Gamma[3+2+\[Epsilon]-r]Gamma[r+s-2-\[Epsilon]]+1/2 p.p(g[\[Mu],\[Nu]1]momentum[p,\[Sigma]]+g[\[Nu]1,\[Sigma]]momentum[p,\[Mu]]+g[\[Sigma],\[Mu]]momentum[p,\[Nu]1])Gamma[2+\[Epsilon]+1-s]Gamma[2+\[Epsilon]+3-1-r]Gamma[r+s-2-\[Epsilon]-1])
End[]
Protect @@ Names["integrals`*"];
EndPackage[]






