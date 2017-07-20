(* ::Package:: *)

BeginPackage["QFTAlgebras`"]
Unprotect @@ Names["QFTAlgebras`*"];
ClearAll @@ Names["QFTAlgebras`*"];

$Assumptions=d\[Element]Complexes;
$Assumptions=m\[Element]Complexes;
$Assumptions=M\[Element]Complexes;
$Assumptions=gs\[Element]Complexes;
Protect[{m,M,gs,d}];

\[Gamma]::usage = "Dirac Matrix. \[Gamma][\[Mu]], \[Gamma][5], \[Gamma][id] all valid options."
\[Gamma]T::usage = "Dirac Matrix transpose."
tr::usage = "Defines a Trace environment for QFT calculations"
\[Delta]::usage = "\[Delta][\[Mu],\[Nu]] defines the Kronecker Delta."
g::usage = "Metric tensor: g[\[Mu],\[Nu]]"
id::usage = "Denotes the dxd identity via \[Gamma][id]"
(*Levi::usage = "Levi Civita symbol in 4D"*)
lam::usage = "Expresses Gell-Man matricies with colour index. lam[a]"
momentum::usage = "momentum[p,\[Mu]] denotes a momentum four vector with Lorentz index \[Mu] and momentum p."
sigma::usage = "\!\(\*SubscriptBox[\(\[Sigma]\), \(\[Mu]\[Nu]\)]\) = \!\(\*FractionBox[\(i\), \(2\)]\)[\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Mu]\)]\),\!\(\*SubscriptBox[\(\[Gamma]\), \(\[Nu]\)]\)]"
C::usage = "Charge conjugation parity operator."
Begin["`Private`"]
(*Kronecker Delta in color-space*)
\[Delta][\[Mu]_,\[Nu]_]/;\[Mu]\[Element]Integers&&\[Nu]\[Element]Integers:=Piecewise[{\[Mu]=\[Nu],1},0]
\[Delta]/:j_[a___,\[Mu]_,b___,\[Nu]_,c___]\[Delta][\[Mu]_,\[Rho]_]/;FreeQ[j,\[Gamma]]:=j[a,\[Rho],b,\[Nu],c]
\[Delta]/:j_[a___,\[Mu]_,b___,\[Nu]_,c___]\[Delta][\[Rho]_,\[Mu]_]/;FreeQ[j,\[Gamma]]:=j[a,\[Rho],b,\[Nu],c]
\[Delta]/:j_[a___,\[Mu]_,b___,\[Nu]_,c___]\[Delta][\[Nu]_,\[Rho]_]/;FreeQ[j,\[Gamma]]:=j[a,\[Mu],b,\[Rho],c]
\[Delta]/:j_[a___,\[Mu]_,b___,\[Nu]_,c___]\[Delta][\[Rho]_,\[Nu]_]/;FreeQ[j,\[Gamma]]:=j[a,\[Mu],b,\[Rho],c]
\[Delta]/:\[Delta][a_,b_]^2:=\[Delta][a,a]
\[Delta]/:\[Delta][a_,a_]:=8
(*Gamma Matrices*)
\[Gamma]/:\[Gamma][\[Mu]_]**\[Gamma][\[Mu]_]/;!SameQ[\[Mu],5]:=d
\[Gamma]/:\[Gamma][\[Mu]_]**\[Gamma][\[Nu]_]**\[Gamma][\[Mu]_]/;!SameQ[\[Mu],5]&&!SameQ[\[Nu],5]:=(2-d)\[Gamma][\[Nu]]
\[Gamma]/:\[Gamma][\[Mu]_]**\[Gamma][\[Nu]_]**\[Gamma][\[Lambda]_]**\[Gamma][\[Mu]_]/;!SameQ[\[Mu],5]&&!SameQ[\[Nu],5]&&!SameQ[\[Lambda],5]:=d*g[\[Nu],\[Lambda]]
\[Gamma]/:\[Gamma][\[Mu]_]**\[Gamma][\[Nu]_]**\[Gamma][\[Lambda]_]**\[Gamma][\[Sigma]_]**\[Gamma][\[Mu]_]/;!SameQ[\[Mu],5]&&!SameQ[\[Nu],5]&&!SameQ[\[Lambda],5]&&!SameQ[\[Sigma],5]:=(2-d)\[Gamma][\[Nu]]**\[Gamma][\[Lambda]]**\[Gamma][\[Sigma]]
\[Gamma]/:\[Gamma][id]**A_/;FreeQ[A,_Complex]:=A
\[Gamma]/:A_**\[Gamma][id]/;FreeQ[A,_Complex]:=A
\[Gamma]/:A_**\[Gamma][id]**B_:=A**B
\[Gamma]/:\[Gamma][5]**\[Gamma][5]:=\[Gamma][id]
tr/:tr[\[Gamma][5]]:=0
tr/:tr[\[Gamma][5]**\[Gamma][\[Mu]_]**\[Gamma][\[Nu]_]]/;!SameQ[\[Mu],5]&&!SameQ[\[Nu],5]:=0
tr/:tr[\[Gamma][5]**\[Gamma][\[Mu]_]**\[Gamma][\[Nu]_]**\[Gamma][\[Lambda]_]**\[Gamma][\[Sigma]_]]/;!SameQ[\[Mu],5]&&!SameQ[\[Nu],5]&&!SameQ[\[Lambda],5]&&!SameQ[\[Sigma],5]:=4*I*Levi[\[Mu],\[Nu],\[Lambda],\[Sigma]]
tr/:tr[\[Gamma][id]]:=d
tr/:tr[NonCommutativeMultiply[a:\[Gamma][_]..,\[Gamma][5],b___]]/;!MemberQ[{a},\[Gamma][5]]:=((-1)^(Length[{a}]))tr[\[Gamma][5]**a**b]
tr/:tr[NonCommutativeMultiply[a___,\[Gamma][5],b:\[Gamma][_]..,\[Gamma][5],c___]]/;!MemberQ[{b},\[Gamma][5]]:=((-1)^(Length[{b}]))tr[a**b**c]
sigma[\[Mu]_,\[Nu]_]:=(I/2)*(\[Gamma][\[Mu]]**\[Gamma][\[Nu]]-\[Gamma][\[Nu]]**\[Gamma][\[Mu]])
(*Define a new Trace*)
tr/:tr[x_+y_]:=tr[x]+tr[y]
tr/:tr[\[Alpha]_ **A__]/;NumberQ[\[Alpha]]:=\[Alpha] tr[A]
(*Trace Theorems*)
tr/:tr[A__]/;!MemberQ[A,\[Gamma][5]]&&FreeQ[A,lam[__]] &&FreeQ[A,m]&&FreeQ[A,M]&& FreeQ[A,_Complex]&& FreeQ[A,C]&& FreeQ[A,\[Gamma]T] && OddQ[Length[A]] && MemberQ[A,Repeated[\[Gamma][_]]]:=0
tr/:tr[\[Gamma][a_]**\[Gamma][b_]]/;!IntegerQ[a]&&!IntegerQ[b]:=4*g[a,b]
tr/:tr[x__]/;EvenQ[Length[x]]&&Length[x]>2&&!MemberQ[x,\[Gamma][5]]&&MemberQ[x,Repeated[\[Gamma][_]]]&&FreeQ[x,lam[__]]&&FreeQ[x,C]&&FreeQ[x,\[Gamma]T[_]] && FreeQ[x,_Complex]&&FreeQ[x,momentum[__]]&&FreeQ[x,m]&&FreeQ[x,M]:=Sum[((-1)^n)*g[Level[x[[1]],1][[1]],Level[x[[n]],1][[1]]]*tr[Drop[Drop[x,{1}],{n-1}]],{n,2,Length[x]}]
tr/:tr[A___**a_*B_**C___]/;MemberQ[B,\[Delta][_,_]]||MemberQ[B,momentum[_,_]]||MemberQ[B,d]:=B*tr[A**a**C]
tr/:tr[\[Gamma][a_]]:=0
(*Charge Conjugation Operator*)
tr/:tr[A___**C**C**B___]:=-tr[A**B]
(*tr/:tr[A___**C**\[Gamma]T[\[Mu]_]**C**B___]:=tr[A**\[Gamma][\[Mu]]**B]*)
tr/:tr[A___**C**\[Gamma][\[Mu]_]**C**B___]:=tr[A**\[Gamma]T[\[Mu]]**B]
(*Non-commutative Multiply*)
Unprotect[NonCommutativeMultiply];
ClearAll[NonCommutativeMultiply]
NonCommutativeMultiply[a___,n_?NumericQ*c_,b___]:=n a**c**b
NonCommutativeMultiply[a___,d*c_,b___]:=d*a**c**b
NonCommutativeMultiply[a___,0,b___]:=0
NonCommutativeMultiply[]:=\[Gamma][id]
NonCommutativeMultiply/:NonCommutativeMultiply[A_]:=A
NonCommutativeMultiply/:NonCommutativeMultiply[A___,B_?NumericQ,C_?NumericQ,D___]:=B*C NonCommutativeMultiply[A,D]
NonCommutativeMultiply/:NonCommutativeMultiply[Z___,A___/B__,X___]:=NonCommutativeMultiply[Z,A,X]/B
NonCommutativeMultiply/:A___**(gs*C__**E___)**B___:=gs*A**C**E**B 
NonCommutativeMultiply/:NonCommutativeMultiply[A___,lam[a_,b_,c_],B__]:=lam[a,b,c]A**B 
NonCommutativeMultiply/:NonCommutativeMultiply[A___,momentum[a_,b_],B___]:=momentum[a,b]A**B 
NonCommutativeMultiply/:NonCommutativeMultiply[A___,\[Delta][a_,b_],B___]:=\[Delta][a,b]A**B 
NonCommutativeMultiply/:NonCommutativeMultiply[A___,g[a_,b_],B___]:=g[a,b]A**B 
NonCommutativeMultiply/:NonCommutativeMultiply[A___,lam[a_,b_,c_]*C__,B___]:=lam[a,b,c]A**C**B
NonCommutativeMultiply/:NonCommutativeMultiply[A___,\[Delta][a_,b_]*C__,B___]:=\[Delta][a,b]A**C**B
NonCommutativeMultiply/:NonCommutativeMultiply[A___,g[a_,b_]*C__,B___]:=g[a,b]A**C**B
NonCommutativeMultiply/:NonCommutativeMultiply[A___,Levi[a_,b_,c_,d_]*C__,B___]:=Levi[a,b,c,d]A**C**B
NonCommutativeMultiply/:F___**(A_+B_)**C___:=F**A**C+F**B**C
NonCommutativeMultiply/:NonCommutativeMultiply[A__,B_,C___]/;MemberQ[{B},g[_,_]]||MemberQ[{B},lam[_,_,_]]||MemberQ[{B},\[Delta][_,_]]||MemberQ[{B},Levi[_,_,_,_]]:=Times[B,NonCommutativeMultiply[A,C]]
NonCommutativeMultiply/:NonCommutativeMultiply[A___,Times[x_,C_],D___]/;SameQ[x,gs]||SameQ[x,m]||SameQ[x,M]||MemberQ[{x},momentum[_,_]]:=Times[x,NonCommutativeMultiply[A,C,D]]
NonCommutativeMultiply/:NonCommutativeMultiply[A___,momentum[p_,\[Mu]_]\[Gamma][\[Mu]_],B___]:=momentum[p,\[Mu]]*NonCommutativeMultiply[A,\[Gamma][\[Mu]],B]
SetAttributes[NonCommutativeMultiply,{Flat,OneIdentity}]
Protect[NonCommutativeMultiply];
(*Four Vectors*)
tr/:tr[A___**p_**B___]/;FreeQ[p,\[Gamma][_]]&&FreeQ[p,lam[_]]&&FreeQ[p,\[Gamma][id]]&&FreeQ[p,C]&&FreeQ[p,\[Gamma]T[_]]:=p tr[A**B]
tr/:tr[A___*p_*B___]/;FreeQ[p,\[Gamma][_]]&&FreeQ[p,lam[_]]&&FreeQ[p,\[Gamma][id]]&&FreeQ[p,C]&&FreeQ[p,\[Gamma]T[_]]:=p tr[A**B]
momentum/:momentum[a_*p_,\[Nu]_]/;a\[Element]Complexes:=a*momentum[p,\[Nu]]
momentum/:momentum[p_,\[Mu]_]momentum[q_,\[Mu]_]:=p.q
momentum/:momentum[p_,\[Mu]_]^2:=p.p
momentum/:momentum[p_+q_,\[Mu]_]:=momentum[p,\[Mu]]+momentum[q,\[Mu]]
(*Metric Identites*)
g/:g[\[Mu]_,\[Mu]_]:=d
g/:g[\[Mu]_,\[Lambda]_]g[\[Lambda]_,\[Nu]_]:=g[\[Mu],\[Nu]]
g/:g[\[Mu]_,\[Lambda]_]g[\[Nu]_,\[Lambda]_]:=g[\[Mu],\[Nu]]
g/:g[\[Lambda]_,\[Mu]_]g[\[Lambda]_,\[Nu]_]:=g[\[Mu],\[Nu]]
g/:g[\[Mu]_,\[Nu]_] momentum[p_,\[Nu]_] :=momentum[p,\[Mu]] 
g/:g[\[Mu]_,\[Nu]_] momentum[p_,\[Mu]_]:=momentum[p,\[Nu]] 
g/:g[\[Mu]_,\[Nu]_] momentum[p_,\[Nu]_] momentum[k_,\[Mu]_]:=momentum[p,\[Nu]] momentum[k,\[Nu]]
g/:g[\[Mu]_,\[Rho]_] momentum[p_,\[Nu]_] momentum[k_,\[Rho]_]:=momentum[p,\[Nu]] momentum[k,\[Mu]]
g/:g[\[Mu]_,\[Nu]_]g[\[Mu]_,\[Nu]_]:=g[\[Mu],\[Mu]]
g/:g[\[Mu]_,\[Nu]_]^2:=g[\[Mu],\[Mu]]
g/:g[\[Mu]_,\[Nu]_]NonCommutativeMultiply[\[Gamma][\[Nu]_],A__]/;!SameQ[\[Nu],5]:=NonCommutativeMultiply[\[Gamma][\[Mu]],A]
(*Gell-Mann Matrices*)
lam/:lam[a_,\[Alpha]1_,\[Alpha]4_]lam[b_,\[Alpha]4_,\[Alpha]1_]:=tr[lam[a]**lam[b]]
lam/:lam[a_,\[Alpha]1_,\[Alpha]2_]lam[b_,\[Alpha]2_,\[Alpha]3_]lam[c_,\[Alpha]3_,\[Alpha]1_]:=tr[lam[a]**lam[b]**lam[c]]
tr/:tr[A___**lam[a_,b_,c_]**B___]:=lam[a,b,c]tr[A**B]

tr/:tr[lam[a_]**lam[b_]]:=2\[Delta][a,b] 
tr/:tr[lam[a_]**lam[b_]**lam[c_]]:=2(dSym[a,b,c]+I*f[a,b,c])
tr/:tr[lam[a_]**lam[b_]**lam[c_]**lam[d_]]:=(4/3)(\[Delta][a,b]\[Delta][c,d]-\[Delta][a,c]\[Delta][b,d]+\[Delta][a,d]\[Delta][b,c])+2(dSym[a,b,r]dSym[c,d,r]-dSym[a,c,r]dSym[d,b,r]+dSym[a,d,r]dSym[b,c,r])+2I(dSym[a,b,r]f[c,a,r]-dSym[a,c,r]f[a,b,r]+dSym[a,d,r]f[b,c,r])

lam/:lam[a_,\[Alpha]1_,\[Alpha]4_]**\[Gamma][id]:=lam[a,\[Alpha]1,\[Alpha]4]
lam/:\[Gamma][id]**lam[a_,\[Alpha]1_,\[Alpha]4_]:=lam[a,\[Alpha]1,\[Alpha]4]
(*Structure Constants - INCOMPLETE*)
f[a_,b_,c_]/;a\[Element]Integers&&b\[Element]Integers&&c\[Element]Integers:=Signature[{a,b,c}]*Piecewise[{{1,Sort[{a,b,c}][[3]]==3},{Sqrt[3]/2,Sort[{a,b,c}][[3]]==8},{-1/2,Sort[{a,b,c}]=={3,6,7}||Sort[{a,b,c}]=={1,5,6}},{1/2,Sort[{a,b,c}]=={1,4,7}||Sort[{a,b,c}]=={2,4,6}||Sort[{a,b,c}]=={2,5,7}||Sort[{a,b,c}]=={3,4,5}}},0]

(*dSym[a_,b_,c_]/;a\[Element]Integers&&b\[Element]Integers&&c\[Element]Integers:=tr[(1/2)(lam[a].lam[b]+lam[b].lam[a]-(4/3)\[Delta][a,b]**\[Gamma][id]).lam[c]]/tr[lam[c]**lam[c]];*)

f/:f[a_,b_,c_]/;Signature[{a,b,c}]==-1:=Signature[{a,b,c}]f[Sort[{a,b,c}][[1]],Sort[{a,b,c}][[2]],Sort[{a,b,c}][[3]]]
(*Double Contraction*)
f/:f[a_,b_,c_]f[d_,b_,c_]:=3\[Delta][a,d]
f/:f[a_,b_,c_]f[a_,b_,d_]:=3\[Delta][c,d]
f/:f[a_,b_,c_]f[a_,c_,d_]:=3\[Delta][b,d]
(*Single Contraction*)
(*f/:f[a_,b_,r_]f[c_,d_,r_]/;!SameQ[b,d]||!SameQ[a,c]||!SameQ[b,c]||!SameQ[a,d]:=(2/3)(\[Delta][a,c]\[Delta][b,d]-\[Delta][a,d]\[Delta][b,c])+dSym[a,c,r]dSym[d,b,r]-dSym[a,d,r]dSym[b,c,r]
f/:f[r_,a_,b_]f[r_,c_,d_]/;!SameQ[b,d]||!SameQ[a,c]||!SameQ[b,c]||!SameQ[a,d]:=(2/3)(\[Delta][a,c]\[Delta][b,d]-\[Delta][a,d]\[Delta][b,c])+dSym[a,c,r]dSym[d,b,r]-dSym[a,d,r]dSym[b,c,r]*)

(*StrucSum1[a_,b_,c_,d_,r_]:=f[a,b,r]f[c,d,r]+f[a,c,r]f[d,b,r]+f[a,d,r]f[b,c,r]
StrucSum1/:StrucSum1[a_,b_,c_,d_,r_]:=0

StrucSum2[a_,b_,c_,d_,r_]:=f[a,b,r]dSym[c,d,r]+f[a,c,r]dSym[d,b,r]+f[a,d,r]dSym[b,c,r]
StrucSum2/:StrucSum2[a_,b_,c_,d_,r_]:=0

dSym/:dSym[a_,b_,c_]/;Signature[{a,b,c}]==0:=0
dSym/:dSym[a_,b_,c_]dSym[d_,b_,c_]:=(3-3/4)\[Delta][a,b]*)
dSym/:dSym[a_,b_,c_]/;Signature[{a,b,c}]==0:=0
(*dSym/:dSym[a_,b_,c_]dSym[d_,b_,c_]:=(3-3/4)\[Delta][a,b]*)

(*Levi-Civita Tensor*)
(*Defined as if all indices were up*)
Levi/:Levi[\[Mu]_,\[Nu]_,\[Rho]_,\[Sigma]_]Levi[\[Mu]_,\[Upsilon]_,\[Eta]_,\[Lambda]_]/;!SameQ[\[Nu],\[Upsilon]]:=-(d-3)(g[\[Nu],\[Upsilon]](g[\[Rho],\[Eta]]g[\[Sigma],\[Lambda]]-g[\[Rho],\[Lambda]]g[\[Sigma],\[Eta]])+g[\[Nu],\[Eta]](g[\[Rho],\[Lambda]]g[\[Sigma],\[Upsilon]]-g[\[Rho],\[Upsilon]]g[\[Sigma],\[Lambda]])+g[\[Nu],\[Lambda]](g[\[Rho],\[Upsilon]]g[\[Sigma],\[Eta]]-g[\[Rho],\[Eta]]g[\[Sigma],\[Upsilon]]))
Levi/:Levi[\[Mu]_,\[Nu]_,\[Rho]_,\[Sigma]_]Levi[\[Mu]_,\[Nu]_,\[Eta]_,\[Lambda]_]:=-(d-3)(d-2)(g[\[Rho],\[Eta]]g[\[Sigma],\[Lambda]]-g[\[Rho],\[Lambda]]g[\[Sigma],\[Eta]])
Levi/:Levi[a___,\[Mu]_,b___]Levi[c___,\[Mu]_,d___]/;!SameQ[a,c]:=(-1)^(Length[List[a]]+Length[List[c]])*Levi[\[Mu],a,b]Levi[\[Mu],c,d]
Levi/:Levi[a__]/;Signature[{a}]===0:=0
Levi/:Levi[\[Mu]_,\[Nu]_,\[Rho]_,\[Sigma]_]/;Signature[{\[Mu],\[Nu],\[Rho],\[Sigma]}]==-1:=Signature[{\[Mu],\[Nu],\[Rho],\[Sigma]}]Levi[Sort[{\[Mu],\[Nu],\[Rho],\[Sigma]}][[1]],Sort[{\[Mu],\[Nu],\[Rho],\[Sigma]}][[2]],Sort[{\[Mu],\[Nu],\[Rho],\[Sigma]}][[3]],Sort[{\[Mu],\[Nu],\[Rho],\[Sigma]}][[4]]]


Levi/:g[a_,b_]Levi[\[Alpha]_,\[Beta]_,\[Delta]_,\[Sigma]_]/;MemberQ[{\[Alpha],\[Beta],\[Delta],\[Sigma]},a]:=Piecewise[{{Levi[\[Alpha],\[Beta],\[Delta],\[Sigma]]/.\[Alpha]->b, a===\[Alpha]},{Levi[\[Alpha],\[Beta],\[Delta],\[Sigma]]/.\[Beta]->b,a===\[Beta]},{Levi[\[Alpha],\[Beta],\[Delta],\[Sigma]]/.\[Delta]->b,a===\[Delta]},{Levi[\[Alpha],\[Beta],\[Delta],\[Sigma]]/.\[Sigma]->b,a===\[Sigma]}},0]
Levi/:g[a_,b_]Levi[\[Alpha]_,\[Beta]_,\[Delta]_,\[Sigma]_]/;MemberQ[{\[Alpha],\[Beta],\[Delta],\[Sigma]},b]:=Piecewise[{{Levi[\[Alpha],\[Beta],\[Delta],\[Sigma]]/.\[Alpha]->a,b===\[Alpha]},{Levi[\[Alpha],\[Beta],\[Delta],\[Sigma]]/.\[Beta]->a,b===\[Beta]},{Levi[\[Alpha],\[Beta],\[Delta],\[Sigma]]/.\[Delta]->a,b===\[Delta]},{Levi[\[Alpha],\[Beta],\[Delta],\[Sigma]]/.\[Sigma]->a,b===\[Sigma]}},0]
(*Custom Commands and Rules*)
momentumConvert:={momentum[p_,\[Mu]_]->p[\[Mu]], p.p->p^2, q.q->q^2,k.k->k^2}
momentumUnConvert:={p_[\[Mu]_]->momentum[p,\[Mu]], p^2->p.p, q^2->q.q,k^2->k.k}
toEpsilon:= {d->4+2\[Epsilon]}
(*result from Davydychev and Boos*)
(*conventions changed oct 17 2013 such that C(-M^2)^\[Alpha] defined with C(-1)^\[Alpha](M^2)^\[Alpha]*)
evalInt:={(*TJI[d,Dot[p_,p_]|Power[p_,2],List[List[\[Nu]_,M_],List[1,0],List[1,0]]]\[RuleDelayed] ((-1)^(2+\[Epsilon]+1)/(I^(4+2\[Epsilon])))(M^2)^(2+2\[Epsilon]-\[Nu])(-1)^(2+2\[Epsilon]-\[Nu])((Gamma[2-1+\[Epsilon]]Gamma[1+\[Nu]-3-2\[Epsilon]]Gamma[1+\[Epsilon]]Gamma[-\[Epsilon]])/(Gamma[2+\[Epsilon]]Gamma[\[Nu]]))*(HoldForm[HypergeometricPFQ[{-\[Epsilon],1+\[Nu]-3-2\[Epsilon]},{2+\[Epsilon]},(p/M)^2]])*)TJI[d,Dot[p_,p_]|Power[p_,2],List[List[\[Nu]_,M_/;!NumberQ[M]],List[1,0],List[1,0]]]:> (-1)(M^2)^(2+2\[Epsilon]-\[Nu])(-1)^(-\[Nu])((Gamma[1+\[Epsilon]]Gamma[\[Nu]-2-2\[Epsilon]]Gamma[1+\[Epsilon]]Gamma[-\[Epsilon]])/(Gamma[2+\[Epsilon]]Gamma[\[Nu]]))*(HoldForm[HypergeometricPFQ[{-\[Epsilon],\[Nu]-2-2\[Epsilon]},{2+\[Epsilon]},(p/M)^2]]),TBI[d,Dot[q_,q_]|Power[q_,2],List[List[\[Beta]_,M_/;!NumberQ[M]],List[\[Alpha]_,0]]]:>I*(-1)^(-\[Alpha]-\[Beta])*(M^2)^(d/2-\[Alpha]-\[Beta])*((Gamma[d/2-\[Beta]]Gamma[\[Alpha]+\[Beta]-d/2])/(Gamma[d/2]Gamma[\[Alpha]]))*HypergeometricPFQ[{\[Beta],\[Alpha]+\[Beta]-d/2},{d/2},q^2/M^2],
TAI[d,0,{{1,M_}}]:>-I*(M^2)^(d/2-1)*Gamma[1-d/2],
TJI[d,Dot[p_,p_]|Power[p_,2],List[List[\[Nu]_,0],List[1,0],List[1,0]]]:> - ((-p^2)^(d-4)/(p^2)^(\[Nu]-2)) (Gamma[\[Nu]+2-d]Gamma[d/2-1]Gamma[d/2-1]Gamma[d/2-\[Nu]])/(Gamma[\[Nu]]Gamma[1]Gamma[1]Gamma[3/2 d-\[Nu]-2]),
TBI[d,Dot[q_,q_]|Power[q_,2],List[List[\[Beta]_,0],List[\[Alpha]_,0]]]:>I (-q^2)^(d/2-2)/(q^2)^(\[Beta]+\[Alpha]-2) (Gamma[d/2-\[Beta]]Gamma[d/2-\[Alpha]]Gamma[\[Alpha]+\[Beta]-d/2])/(Gamma[\[Beta]]Gamma[\[Alpha]]Gamma[d-\[Alpha]-\[Beta]])
(*TBI[d,Dot[q_,q_]|Power[q_,2],List[List[\[Beta]_,M_],List[\[Alpha]_,M_]]]:>I^(1-d)*(-1)^(d/2-\[Alpha]-\[Beta])*(M^2)^(d/2-\[Alpha]-\[Beta])*((Gamma[\[Alpha]+\[Beta]-d/2])/(Gamma[\[Alpha]+\[Beta]]))*HypergeometricPFQ[{\[Alpha],\[Beta],\[Alpha]+\[Beta]-d/2},{(\[Alpha]+\[Beta])/2,(\[Alpha]+\[Beta]+1)/2},q^2/(4M^2)]*)}
evalDerivative:={D[TJI[d,Power[p_,2],List[List[\[Nu]_,M_],List[\[Mu]_,m_],List[1,0]]],{m_,1}]:> 2\[Mu]*m(TJI[d,Power[p,2],List[List[\[Nu],M],List[\[Mu]+1,m],List[1,0]]]) ,D[TJI[d,Power[p_,2],List[List[\[Nu]_,M_],List[\[Mu]_,m_],List[1,0]]],{m_,2}]:> \[Mu](TJI[d,Power[p,2],List[List[\[Nu],M],List[1+\[Mu],m],List[1,0]]]+m(TJI[d,Power[p,2],List[List[\[Nu],M],List[\[Mu]+2,m],List[1,0]]])) }
(*expandPFQ:={Hold[HypergeometricPFQ[{A__},{B__},z_]]\[Rule](1+((Product[Gamma[{B}[[i]]],{i,1,Length[{B}]}])/(Product[Gamma[{A}[[i]]],{i,1,Length[{A}]}]))
Sum[
((Product[Gamma[{A}[[i]]],{i,1,Length[{A}]}])/(Product[Gamma[{B}[[i]]],{i,1,Length[{B}]}]))*(z^n/n!),{n,1,Infinity}])}*)
(*expandPFQ:={HypergeometricPFQ[{A__},{B__},z_]:>(1+((Product[Gamma[{B}[[i]]],{i,1,Length[{B}]}])/(Product[Gamma[{A}[[i]]],{i,1,Length[{A}]}]))
Sum[
((Product[Gamma[{A}[[i]]+n],{i,1,Length[{A}]}])/(Product[Gamma[{B}[[i]]+n],{i,1,Length[{B}]}]))*(Hold[z^n/n!]),{n,1,Infinity}])}
expandPFQ2:={HypergeometricPFQ[{A__},{B__},z_]:>(1+((Product[Gamma[{B}[[i]]],{i,1,Length[{B}]}])/(Product[Gamma[{A}[[i]]],{i,1,Length[{A}]}]))*((Product[Gamma[{A}[[i]]+1],{i,1,Length[{A}]}])/(Product[Gamma[{B}[[i]]+1],{i,1,Length[{B}]}]))*(z^1/1!)+
((Product[Gamma[{B}[[i]]],{i,1,Length[{B}]}])/(Product[Gamma[{A}[[i]]],{i,1,Length[{A}]}]))Sum[
((Product[Gamma[{A}[[i]]+n],{i,1,Length[{A}]}])/(Product[Gamma[{B}[[i]]+n],{i,1,Length[{B}]}]))*Hold[(z^n/n!)],{n,2,Infinity}])}
expand2F1:={Hold[Hypergeometric2F1[A_,B_,C_,z_]]:>(1+((Gamma[C]Gamma[A+1]Gamma[B+1])/(Gamma[A]Gamma[B]Gamma[C+1]))(z)+
((Gamma[C])/(Gamma[A]Gamma[B]))Sum[
((Gamma[A+n]Gamma[B+n])/Gamma[C+n])*(z^n/n!),{n,2,Infinity}])}*)

(*Formatting*)
Unprotect[KroneckerDelta];
Format[\[Delta][a_,b_],StandardForm]:=Superscript["\[Delta]",{a,b}]
Protect[KroneckerDelta];
Format[f[a_,b_,c_],StandardForm]:=Superscript["f",{a,b,c}]
Format[momentum[p_,\[Mu]_],StandardForm]:=Subscript[Style[p,Bold],\[Mu]]
Format[lam[a_,\[Alpha]1_,\[Alpha]2_],StandardForm]:=Subscript[Superscript["\[Lambda]",a],{\[Alpha]1,\[Alpha]2}]
Format[Levi[a_,b_,c_,d_],StandardForm]:=Superscript["\[Epsilon]",{a,b,c,d}]
Format[g[a_,b_],StandardForm]:=Subscript["g",{a,b}]
Format[\[Gamma][b_],StandardForm]:=Superscript[Style["\[Gamma]",Bold],{b}]

End[]
Protect @@ Names["QFTAlgebras`*"];
EndPackage[]







