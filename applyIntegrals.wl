(* ::Package:: *)

BeginPackage["applyIntegrals`"]
Needs[""]

appyIntegrals::usage = "Massless integration procedure. Syntax applyIntegrals[expression_,rootarg_,arg1_,arg2_]"
momentum::usage = "Momentum commands as defined in QFTpackage"
d::usage = "Number of dimensions in dim-reg scheme. Assume d = 4 + 2\[Epsilon]"
\[Epsilon]::usage = "Expansion parameter in dim-reg scheme."
g::usage = "metric"
integralDef::usage = "Massless integral definitions for rank-0 to rank-3 type numerator structure."
\[Nu]::usage = "Renormalization scale."

Begin["`Private`"]

(*all checked-2018/04/11*)
(*defined with relative minus sign between p and q*)
integralDef[p_,r_,s_]:=(I \[Nu]^(2\[Epsilon]))/(4\[Pi])^2 (-p.p/(4\[Pi] \[Nu]^2))^\[Epsilon] (2\[Pi])^d/(p.p)^(r+s-2) (Gamma[2-r+\[Epsilon]]Gamma[2-s+\[Epsilon]]Gamma[r+s-2-\[Epsilon]])/(Gamma[r]Gamma[s]Gamma[4-r-s+2\[Epsilon]])
integralDef[p_,r_,s_,\[Mu]_]:=(I \[Nu]^(2\[Epsilon]))/(4\[Pi])^2 (-p.p/(4\[Pi] \[Nu]^2))^\[Epsilon] (2\[Pi])^d/(p.p)^(r+s-2) momentum[p,\[Mu]] (Gamma[3-r+\[Epsilon]]Gamma[2-s+\[Epsilon]]Gamma[r+s-2-\[Epsilon]])/(Gamma[r]Gamma[s]Gamma[5-r-s+2\[Epsilon]])
integralDef[p_,r_,s_,\[Mu]_,\[Nu]1_]:=(I \[Nu]^(2\[Epsilon]))/(4\[Pi])^2 (-p.p/(4\[Pi] \[Nu]^2))^\[Epsilon] (2\[Pi])^d/(p.p)^(r+s-2) (g[\[Mu],\[Nu]1]p.p (Gamma[3-r+\[Epsilon]]Gamma[3-s+\[Epsilon]]Gamma[r+s-3-\[Epsilon]])/(2 Gamma[r]Gamma[s]Gamma[6-r-s+2\[Epsilon]])+momentum[p,\[Mu]]momentum[p,\[Nu]1] (Gamma[4-r+\[Epsilon]]Gamma[2-s+\[Epsilon]]Gamma[r+s-2-\[Epsilon]])/(Gamma[r]Gamma[s]Gamma[6-r-s+2\[Epsilon]]))
integralDef[p_,r_,s_,\[Mu]_,\[Nu]1_,\[Sigma]_]:=
((\[Pi]^(2+\[Epsilon]) I^(1-4-2\[Epsilon]))(p.p)^(2+\[Epsilon]-r-s))/(Gamma[r]Gamma[s]Gamma[4+2\[Epsilon]-r-s+3]) \[Nu]^(2\[Epsilon])/(\[Nu]^2)^\[Epsilon] (momentum[p,\[Mu]]momentum[p,\[Nu]1]momentum[p,\[Sigma]] Gamma[2+\[Epsilon]-s]Gamma[3+2+\[Epsilon]-r]Gamma[r+s-2-\[Epsilon]]+1/2 p.p(g[\[Mu],\[Nu]1]momentum[p,\[Sigma]]+g[\[Nu]1,\[Sigma]]momentum[p,\[Mu]]+g[\[Sigma],\[Mu]]momentum[p,\[Nu]1])Gamma[2+\[Epsilon]+1-s]Gamma[2+\[Epsilon]+3-1-r]Gamma[r+s-2-\[Epsilon]-1])
integralDef[p_,r_,s_,\[Alpha]1_,\[Alpha]2_,\[Alpha]3_,\[Alpha]4_]:=
((\[Pi]^(2+\[Epsilon]) I^(1-4-2\[Epsilon]))(p.p)^(2+\[Epsilon]-r-s))/(Gamma[r]Gamma[s]Gamma[4+2\[Epsilon]-r-s+4]) \[Nu]^(2\[Epsilon])/(\[Nu]^2)^\[Epsilon] (momentum[p,\[Alpha]1]momentum[p,\[Alpha]2]momentum[p,\[Alpha]3]momentum[p,\[Alpha]4] Gamma[2+\[Epsilon]-s]Gamma[4+2+\[Epsilon]-r]Gamma[r+s-2-\[Epsilon]]+1/2 p.p(g[\[Alpha]1,\[Alpha]2]momentum[p,\[Alpha]3]momentum[p,\[Alpha]4]+g[\[Alpha]2,\[Alpha]3]momentum[p,\[Alpha]1]momentum[p,\[Alpha]4]+g[\[Alpha]3,\[Alpha]4]momentum[p,\[Alpha]2]momentum[p,\[Alpha]1]+g[\[Alpha]1,\[Alpha]4]momentum[p,\[Alpha]2]momentum[p,\[Alpha]3])Gamma[2+\[Epsilon]+1-s]Gamma[2+\[Epsilon]+3-r]Gamma[r+s-2-\[Epsilon]-1]+1/4 (p.p)^2 (g[\[Alpha]1,\[Alpha]2]g[\[Alpha]3,\[Alpha]4]+g[\[Alpha]2,\[Alpha]3]g[\[Alpha]1,\[Alpha]4])Gamma[2+\[Epsilon]-s+2]Gamma[2+\[Epsilon]+2-r]Gamma[r+s-2-\[Epsilon]-2])
integralDef[p_,r_,s_,\[Alpha]1_,\[Alpha]2_,\[Alpha]3_,\[Alpha]4_,\[Alpha]5_]:=(*integralDef[p,r,s,\[Alpha]1,\[Alpha]2,\[Alpha]3,\[Alpha]4,\[Alpha]5]*)(-(\[Pi]^(2+\[Epsilon]) I^(1-4-2\[Epsilon])) (p.p)^(2+\[Epsilon]-r-s))/(Gamma[r]Gamma[s]Gamma[4+2\[Epsilon]-r-s+5]) \[Nu]^(2\[Epsilon])/(\[Nu]^2)^\[Epsilon] (momentum[p,\[Alpha]1]momentum[p,\[Alpha]2]momentum[p,\[Alpha]3]momentum[p,\[Alpha]4]momentum[p,\[Alpha]5] Gamma[2+\[Epsilon]-s]Gamma[5+2+\[Epsilon]-r]Gamma[r+s-2-\[Epsilon]]+1/2 p.p(g[\[Alpha]1,\[Alpha]2]momentum[p,\[Alpha]3]momentum[p,\[Alpha]4]momentum[p,\[Alpha]5]+g[\[Alpha]2,\[Alpha]3]momentum[p,\[Alpha]1]momentum[p,\[Alpha]4]momentum[p,\[Alpha]5]+g[\[Alpha]3,\[Alpha]4]momentum[p,\[Alpha]5]momentum[p,\[Alpha]2]momentum[p,\[Alpha]1]+g[\[Alpha]4,\[Alpha]5]momentum[p,\[Alpha]1]momentum[p,\[Alpha]2]momentum[p,\[Alpha]3]+g[\[Alpha]5,\[Alpha]1]momentum[p,\[Alpha]4]momentum[p,\[Alpha]2]momentum[p,\[Alpha]3])Gamma[2+\[Epsilon]+1-s]Gamma[2+\[Epsilon]+4-r]Gamma[r+s-2-\[Epsilon]-1]+1/4 (p.p)^2 ((g[\[Alpha]1,\[Alpha]2]g[\[Alpha]3,\[Alpha]4]+g[\[Alpha]1,\[Alpha]3]g[\[Alpha]2,\[Alpha]4]+g[\[Alpha]1,\[Alpha]4]g[\[Alpha]2,\[Alpha]3])momentum[p,\[Alpha]5]+(g[\[Alpha]2,\[Alpha]3]g[\[Alpha]4,\[Alpha]5]+g[\[Alpha]2,\[Alpha]4]g[\[Alpha]3,\[Alpha]5]+g[\[Alpha]2,\[Alpha]5]g[\[Alpha]3,\[Alpha]4])momentum[p,\[Alpha]1]+(g[\[Alpha]1,\[Alpha]3]g[\[Alpha]4,\[Alpha]5]+g[\[Alpha]1,\[Alpha]4]g[\[Alpha]3,\[Alpha]5]+g[\[Alpha]1,\[Alpha]5]g[\[Alpha]3,\[Alpha]4])momentum[p,\[Alpha]2]+(g[\[Alpha]2,\[Alpha]1]g[\[Alpha]4,\[Alpha]5]+g[\[Alpha]2,\[Alpha]4]g[\[Alpha]1,\[Alpha]5]+g[\[Alpha]2,\[Alpha]5]g[\[Alpha]1,\[Alpha]4])momentum[p,\[Alpha]3]+(g[\[Alpha]2,\[Alpha]3]g[\[Alpha]1,\[Alpha]5]+g[\[Alpha]2,\[Alpha]1]g[\[Alpha]3,\[Alpha]5]+g[\[Alpha]2,\[Alpha]5]g[\[Alpha]3,\[Alpha]1])momentum[p,\[Alpha]4])Gamma[2+\[Epsilon]-s+2]Gamma[2+\[Epsilon]+2-r]Gamma[r+s-2-\[Epsilon]-2])

applyIntegrals[expression_,rootarg_,arg1_,arg2_]:=Module[{denom1,numer1,denom2,numer2,exponentListDenom,exponentListNum,coeffRank,momentumList,indexList,uniqueIndexList,(*index1,index2,index3,*)rootArgVecCount,arg2VecCount,rootArgMomentum,arg2Momentum},
(*find any numerator and denomenator structure with the integration variable*)
denom1=Cases[Denominator[expression],Power[Dot[arg1,arg1],_]|Power[Dot[arg1,_],_]|Power[Dot[_,arg1],_]|Power[Dot[Plus[_. arg1,_],Plus[_. arg1,_]],_]|Dot[arg1,_]|Dot[_,arg1]|Dot[Plus[_. arg1,_],Plus[_. arg1,_]](*|((-Dot[Plus[_. arg1,_],Plus[_. arg1,_]])^_)*),{0,1}];
(*numer1=Cases[Numerator[expression]/.(a_+b_).(a_+b_)\[RuleDelayed]a.a+b.b+2a.b,Power[Dot[arg1,arg1],_]|Dot[arg1,arg1]|Power[Dot[arg1,_],_]|Dot[arg1,_]|Power[Dot[_,arg1],_]|Dot[_,arg1]|momentum[arg1,_],{0,1}];*)
numer1=Cases[Numerator[expression]/.(a_+b_).(a_+b_):>a.a+b.b+2a.b,Power[Dot[arg1,arg1],_]|Power[Dot[arg1,_],_]|Power[Dot[_,arg1],_],{0,1}];
numer2=Join[numer1,Cases[Numerator[expression/(Times@@numer1)],Dot[arg1,arg1]|Dot[arg1,_]|Dot[_,arg1]|momentum[arg1,_],{0,1}]];
exponentListDenom=Exponent[Times@@denom1,{Dot[arg1,arg1],(-rootarg+arg1).(-rootarg+arg1),(rootarg-arg1).(rootarg-arg1)}];
exponentListNum=Exponent[Times@@numer2,{Dot[arg1,arg1],Dot[arg1,rootarg],Dot[rootarg,arg1],Dot[arg2,arg1],Dot[arg1,arg2]}];
(*for four-vectors*)
momentumList=Cases[Numerator[expression],momentum[arg1,_],{0,1}];
indexList=Table[momentumList[[i,2]],{i,1,Length[momentumList]}];
(*{denom1,numer1, exponentListDenom,exponentListNum};*)
coeffRank=2exponentListNum[[1]]+exponentListNum[[2]]+exponentListNum[[3]]+exponentListNum[[4]]+exponentListNum[[5]]+Length[momentumList];
(*by integrating, we rip apart the dot products and leave four-vectors that match with the integralDef indices.*)
rootArgVecCount=(exponentListNum[[2]]+exponentListNum[[3]]);
arg2VecCount=(exponentListNum[[4]]+exponentListNum[[5]]);
uniqueIndexList={index1,index2,index3,index4,index5, index6};
rootArgMomentum=If[rootArgVecCount>0,Times@@Table[momentum[rootarg,uniqueIndexList[[i]]],{i,1,rootArgVecCount}],1];
arg2Momentum=If[arg2VecCount>0,Times@@Table[momentum[arg2,uniqueIndexList[[i]]],{i,(1+rootArgVecCount),(arg2VecCount+rootArgVecCount)}],1];
expression/.A___ Times@@numer2/Times@@denom1/;FreeQ[{A},arg1]:>
A rootArgMomentum arg2Momentum 
Switch[coeffRank,
0,integralDef[rootarg,exponentListDenom[[1]],exponentListDenom[[2]]+exponentListDenom[[3]]],
1,Switch[Length[momentumList],0,integralDef[rootarg,exponentListDenom[[1]],exponentListDenom[[2]]+exponentListDenom[[3]],index1],1,integralDef[rootarg,exponentListDenom[[1]],exponentListDenom[[2]]+exponentListDenom[[3]],indexList[[1]]]],
2,Switch[Length[momentumList],
0,(g[uniqueIndexList[[rootArgVecCount+arg2VecCount+1]],uniqueIndexList[[rootArgVecCount+arg2VecCount+2]]])^exponentListNum[[1]] integralDef[rootarg,exponentListDenom[[1]],exponentListDenom[[2]]+exponentListDenom[[3]],index1,index2],
1,integralDef[rootarg,exponentListDenom[[1]],exponentListDenom[[2]]+exponentListDenom[[3]],index1,indexList[[1]]],2,integralDef[rootarg,exponentListDenom[[1]],exponentListDenom[[2]]+exponentListDenom[[3]],indexList[[1]],indexList[[2]]]],
3,Switch[Length[momentumList],
0,(g[uniqueIndexList[[rootArgVecCount+arg2VecCount+1]],uniqueIndexList[[rootArgVecCount+arg2VecCount+2]]])^exponentListNum[[1]] integralDef[rootarg,exponentListDenom[[1]],exponentListDenom[[2]]+exponentListDenom[[3]],index1,index2,index3],1,(g[uniqueIndexList[[rootArgVecCount+arg2VecCount+1]],uniqueIndexList[[rootArgVecCount+arg2VecCount+2]]])^exponentListNum[[1]] integralDef[rootarg,exponentListDenom[[1]],exponentListDenom[[2]]+exponentListDenom[[3]],index1,index2,indexList[[1]]],2,integralDef[rootarg,exponentListDenom[[1]],exponentListDenom[[2]]+exponentListDenom[[3]],index1,indexList[[1]],indexList[[2]]],3,integralDef[rootarg,exponentListDenom[[1]],exponentListDenom[[2]]+exponentListDenom[[3]],indexList[[1]],indexList[[2]],indexList[[3]]]],
4,Switch[Length[momentumList],
0,(g[uniqueIndexList[[rootArgVecCount+arg2VecCount+1]],uniqueIndexList[[rootArgVecCount+arg2VecCount+2]]])^Boole[exponentListNum[[1]]>0] (g[uniqueIndexList[[rootArgVecCount+arg2VecCount+3]],uniqueIndexList[[rootArgVecCount+arg2VecCount+4]]])^Boole[exponentListNum[[1]]>1] integralDef[rootarg,exponentListDenom[[1]],exponentListDenom[[2]]+exponentListDenom[[3]],index1,index2,index3,index4],1,(g[uniqueIndexList[[rootArgVecCount+arg2VecCount+1]],uniqueIndexList[[rootArgVecCount+arg2VecCount+2]]])^Boole[exponentListNum[[1]]>0] (g[uniqueIndexList[[rootArgVecCount+arg2VecCount+3]],uniqueIndexList[[rootArgVecCount+arg2VecCount+4]]])^Boole[exponentListNum[[1]]>1] integralDef[rootarg,exponentListDenom[[1]],exponentListDenom[[2]]+exponentListDenom[[3]],index1,index2,index3,indexList[[1]]],2,(g[uniqueIndexList[[rootArgVecCount+arg2VecCount+1]],uniqueIndexList[[rootArgVecCount+arg2VecCount+2]]])^Boole[exponentListNum[[1]]>0] (g[uniqueIndexList[[rootArgVecCount+arg2VecCount+3]],uniqueIndexList[[rootArgVecCount+arg2VecCount+4]]])^Boole[exponentListNum[[1]]>1] integralDef[rootarg,exponentListDenom[[1]],exponentListDenom[[2]]+exponentListDenom[[3]],index1,index2,indexList[[1]],indexList[[2]]],
3,integralDef[rootarg,exponentListDenom[[1]],exponentListDenom[[2]]+exponentListDenom[[3]],index1,indexList[[1]],indexList[[2]],indexList[[3]]],
4,integralDef[rootarg,exponentListDenom[[1]],exponentListDenom[[2]]+exponentListDenom[[3]],indexList[[1]],indexList[[2]],indexList[[3]],indexList[[4]]]],
5,Switch[Length[momentumList],
0,(g[uniqueIndexList[[rootArgVecCount+arg2VecCount+1]],uniqueIndexList[[rootArgVecCount+arg2VecCount+2]]])^Boole[exponentListNum[[1]]>0] (g[uniqueIndexList[[rootArgVecCount+arg2VecCount+3]],uniqueIndexList[[rootArgVecCount+arg2VecCount+4]]])^Boole[exponentListNum[[1]]>1] integralDef[rootarg,exponentListDenom[[1]],exponentListDenom[[2]]+exponentListDenom[[3]],index1,index2,index3,index4,index5],
1,(g[uniqueIndexList[[rootArgVecCount+arg2VecCount+1]],uniqueIndexList[[rootArgVecCount+arg2VecCount+2]]])^Boole[exponentListNum[[1]]>0] (g[uniqueIndexList[[rootArgVecCount+arg2VecCount+3]],uniqueIndexList[[rootArgVecCount+arg2VecCount+4]]])^Boole[exponentListNum[[1]]>1] integralDef[rootarg,exponentListDenom[[1]],exponentListDenom[[2]]+exponentListDenom[[3]],index1,index2,index3,index4,indexList[[1]]],
2,(g[uniqueIndexList[[rootArgVecCount+arg2VecCount+1]],uniqueIndexList[[rootArgVecCount+arg2VecCount+2]]])^Boole[exponentListNum[[1]]>0] (g[uniqueIndexList[[rootArgVecCount+arg2VecCount+3]],uniqueIndexList[[rootArgVecCount+arg2VecCount+4]]])^Boole[exponentListNum[[1]]>1] integralDef[rootarg,exponentListDenom[[1]],exponentListDenom[[2]]+exponentListDenom[[3]],index1,index2,index3,indexList[[1]],indexList[[2]]],
3,integralDef[rootarg,exponentListDenom[[1]],exponentListDenom[[2]]+exponentListDenom[[3]],index1,index2,indexList[[1]],indexList[[2]],indexList[[3]]],
4,integralDef[rootarg,exponentListDenom[[1]],exponentListDenom[[2]]+exponentListDenom[[3]],index1,indexList[[1]],indexList[[2]],indexList[[3]],indexList[[4]]],
5,integralDef[rootarg,exponentListDenom[[1]],exponentListDenom[[2]]+exponentListDenom[[3]],indexList[[1]],indexList[[2]],indexList[[3]],indexList[[4]],indexList[[5]]]]
]
]
applyIntegrals/: applyIntegrals[a_Plus,x_,y_,z_]:=applyIntegrals[#,x,y,z]&/@a

End[]

EndPackage[]


(* ::Input:: *)
(**)
