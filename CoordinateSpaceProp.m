(* ::Package:: *)

BeginPackage["CoordinateSpaceProp`"]
Needs["QFTAlgebras`","/home/hojasonn/Documents/ResearchProjects/RecodingProject/sum-rules/QFTalgebras.m"]

applyIntegrals::usage = "Massless integration procedure. Syntax applyIntegrals[expression_,rootarg_,arg1_,arg2_]"
momentum::usage = "Momentum commands as defined in QFTpackage"
d::usage = "Number of dimensions in dim-reg scheme. Assume d = 4 + 2\[Epsilon]"
\[Epsilon]::usage = "Expansion parameter in dim-reg scheme."
g::usage = "metric"
integralDef::usage = "Massless integral definitions for rank-0 to rank-3 type numerator structure."
\[Nu]::usage = "Renormalization scale."

Begin["`Private`"]
gluonCoordProp[x_,\[Nu]_,\[Rho]_,\[Mu]_,\[Sigma]_]:=(I 2^(d-2) \[Pi]^(d/2) (1-d/2)Gamma[d/2-1])/(2\[Pi])^d (g[\[Rho],\[Sigma]]((-d/2)(-x.x)^(-1-d/2) (4 momentum[x,\[Nu]]momentum[x,\[Mu]])-2g[\[Mu],\[Nu]](-x.x)^(-d/2))-g[\[Rho],\[Mu]]((-d/2)(-x.x)^(-1-d/2) (4 momentum[x,\[Nu]]momentum[x,\[Sigma]])-2g[\[Sigma],\[Nu]](-x.x)^(-d/2))-g[\[Nu],\[Sigma]]((-d/2)(-x.x)^(-1-d/2) (4 momentum[x,\[Rho]]momentum[x,\[Mu]])-2g[\[Mu],\[Rho]](-x.x)^(-d/2))+g[\[Mu],\[Nu]]((-d/2)(-x.x)^(-1-d/2) (4 momentum[x,\[Rho]]momentum[x,\[Sigma]])-2g[\[Sigma],\[Rho]](-x.x)^(-d/2)))//Simplify
fermionCoordProp[x_,index_]:=-((2^(d-1) \[Pi]^(d/2))/(2\[Pi])^d)(1-d/2)(-x.x)^(-d/2) Gamma[d/2-1]momentum[x,index]
contract["BB",coord1_,coord2_,color1_,color2_,lorentz1_,lorentz2_]:=I \[Delta][color1,color2]((I g[lorentz1,lorentz2])/(2\[Pi])^d 2^(d-2) \[Pi]^(d/2) (-(coord1-coord2).(coord1-coord2))^(1-d/2) Gamma[d/2-1]);
contract["GB",coord1_,coord2_,color1_,color2_,lorentz1_,lorentz2_,lorentz3_]:=-2I \[Delta][color1,color2]((I 2^(d-2) \[Pi]^(d/2))/(2\[Pi])^d)(1-d/2)Gamma[d/2-1](-(coord1-coord2).(coord1-coord2))^(-d/2) (g[lorentz3,lorentz1](momentum[coord1,lorentz2]-momentum[coord2,lorentz2])-g[lorentz1,lorentz2](momentum[coord1,lorentz3]-momentum[coord2,lorentz3]))
contract["BG",coord1_,coord2_,color1_,color2_,lorentz1_,lorentz2_,lorentz3_]:=2I \[Delta][color1,color2]((I 2^(d-2) \[Pi]^(d/2))/(2\[Pi])^d)(1-d/2)Gamma[d/2-1](-(coord1-coord2).(coord1-coord2))^(-d/2) (g[lorentz3,lorentz1](momentum[coord1,lorentz2]-momentum[coord2,lorentz2])-g[lorentz1,lorentz2](momentum[coord1,lorentz3]-momentum[coord2,lorentz3]))
contract["GG",coord1_,coord2_,color1_,color2_,lorentz1_,lorentz2_,lorentz3_,lorentz4_]:=I \[Delta][color1,color2] (-I 2^(d-2) \[Pi]^(d/2) (1-d/2))/(2\[Pi])^d (((-4d)/2)(-(coord1-coord2).(coord1-coord2))^(-1-d/2) (momentum[coord1-coord2,lorentz1]momentum[coord1-coord2,lorentz3]g[lorentz2,lorentz4]-momentum[coord1-coord2,lorentz1]momentum[coord1-coord2,lorentz4]g[lorentz2,lorentz3]-momentum[coord1-coord2,lorentz2]momentum[coord1-coord2,lorentz3]g[lorentz1,lorentz4]+momentum[coord1-coord2,lorentz2]momentum[coord1-coord2,lorentz4]g[lorentz1,lorentz3])+(-2)(-(coord1-coord2).(coord1-coord2))^(-(d/2)) (g[lorentz1,lorentz3]g[lorentz2,lorentz4]-g[lorentz1,lorentz4]g[lorentz2,lorentz3]-g[lorentz2,lorentz3]g[lorentz1,lorentz4]+g[lorentz2,lorentz4]g[lorentz1,lorentz3]))Gamma[d/2-1]
End[]

EndPackage[]


(* ::Input:: *)
(**)
