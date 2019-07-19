(* ::Package:: *)

(*Dim-reg expressions for local vacuum expectation values (VEVs)*)
BeginPackage["vevs`"]
Unprotect @@ Names["vevs`*"];
ClearAll @@ Names["vevs`*"];

$Assumptions=d\[Element]Complexes;
$Assumptions=m\[Element]Complexes;
$Assumptions=M\[Element]Complexes;
$Assumptions=gs\[Element]Complexes;
Protect[{m,M,gs,d}];

vev::usage = "Local VEV (all indices contracted)"
GG::usage = "Dimension four gluon condensate."
DDGG::usage = "Expansion of the dim-4 gluon condensate"
GGG::usage = "Dimension six gluon condensate."
QQQQ::usage = "Dimension six quark condensate."
g::usage = "Metric tensor: g[\[Mu],\[Nu]]"

Begin["`Private`"]
GG[lorentz1_, lorentz2_, lorentz3_, lorentz4_]:=(1/(8(d-1)d))(g[lorentz1,lorentz3]g[lorentz2, lorentz4]-g[lorentz1, lorentz4]g[lorentz2,lorentz3])vev[GG]
DDGG[derivindex1_, derivindex2_,lorentz1_, lorentz2_, lorentz3_, lorentz4_]:=
(-1/((d+2)d(d-2)))vev[gggGGG]g[derivindex1,derivindex2]*(g[lorentz1,lorentz3]g[lorentz2,lorentz4]-g[lorentz1,lorentz4]g[lorentz2,lorentz3])+(-1/(2(d+2)d(d-2)))vev[gggGGG](g[derivindex2,lorentz3](g[derivindex1,lorentz1]g[lorentz2,lorentz4]-g[derivindex1,lorentz2]g[lorentz1,lorentz4])
-g[derivindex2,lorentz4](g[derivindex1,lorentz1]g[lorentz2,lorentz3]-g[derivindex1,lorentz2]g[lorentz1,lorentz3]))+(3/(2(d+2)d(d-1)(d-2)))vev[gggGGG](g[derivindex1,lorentz3](g[derivindex2,lorentz1]g[lorentz2,lorentz4]-g[derivindex2,lorentz2]g[lorentz1,lorentz4])
-g[derivindex1,lorentz4](g[derivindex2,lorentz1]g[lorentz2,lorentz3]-g[derivindex2,lorentz2]g[lorentz1,lorentz3]))
GGG[lorentz1_,lorentz2_,lorentz3_,lorentz4_,lorentz5_,lorentz6_]:=(vev[gggGGG]/(d(d-1)(d-2)))*(g[lorentz1,lorentz3]g[lorentz2,lorentz5]g[lorentz6,lorentz4]-g[lorentz2,lorentz3]g[lorentz1,lorentz5]g[lorentz6,lorentz4]-g[lorentz1,lorentz4]g[lorentz2,lorentz5]g[lorentz6,lorentz3]+g[lorentz2,lorentz4]g[lorentz1,lorentz5]g[lorentz6,lorentz3]-g[lorentz1,lorentz3]g[lorentz2,lorentz6]g[lorentz5,lorentz4]+g[lorentz2,lorentz3]g[lorentz1,lorentz6]g[lorentz4,lorentz5]+g[lorentz1,lorentz4]g[lorentz2,lorentz6]g[lorentz5,lorentz3]-g[lorentz2,lorentz4]g[lorentz1,lorentz6]g[lorentz5,lorentz3])
QQQQ[flavor1_,flavor2_,dirac1_,dirac2_,dirac3_,dirac4_,color1_,color2_,color3_,color4_]:=vev[OverBar[flavor1]flavor1]vev[OverBar[flavor2]flavor2] 1/144 (\[Delta][dirac1,dirac4]\[Delta][color1,color4]\[Delta][dirac2,dirac3]\[Delta][color2,color3]-\[Delta][dirac1,dirac3]\[Delta][color1,color3]\[Delta][dirac2,dirac4]\[Delta][color2,color4])

End[]
Protect @@ Names["vevs`*"];
EndPackage[]

