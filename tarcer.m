(* ::Package:: *)

BeginPackage["tarcer`"]
Unprotect @@ Names["tarcer`*"];
ClearAll @@ Names["tarcer`*"];
(*Definitions copied from "2016-10-28 Master Sheet V3.7.nb"*)
d::usage = "Globa dim-reg parameter."
Tarce::usage = "Takes two-loop massive/massless integrals and formats them into 
an appropriate form for TARCER to digest. Expression is defined according to conventions 
summarized in original TARCER paper (https://arxiv.org/pdf/hep-ph/9801383.pdf). 
The integrals are over momenta \!\(\*SubscriptBox[\(k\), \(1\)]\) and \!\(\*SubscriptBox[\(k\), \(2\)]\) with p representing an external momentum; 
momenta \!\(\*SubscriptBox[\(k\), \(3\)]\)-\!\(\*SubscriptBox[\(k\), \(5\)]\) are abbreviations (\!\(\*SubscriptBox[\(k\), \(3\)]\)=\!\(\*SubscriptBox[\(k\), \(1\)]\)-p, \!\(\*SubscriptBox[\(k\), \(4\)]\)=\!\(\*SubscriptBox[\(k\), \(2\)]\)-p, \!\(\*SubscriptBox[\(k\), \(5\)]\)=\!\(\*SubscriptBox[\(k\), \(1\)]\)-\!\(\*SubscriptBox[\(k\), \(2\)]\)). 
Tarce takes inputs Tarce[x_,p_,k1_,k2_,m1_,m2_, m3_, m4_,m5_], where x is the expression, p, \!\(\*SubscriptBox[\(k\), \(1\)]\)-\!\(\*SubscriptBox[\(k\), \(5\)]\) have been described, and \!\(\*SubscriptBox[\(m\), \(1\)]\)-\!\(\*SubscriptBox[\(m\), \(5\)]\) are corresponding masses."
numeratorReplace::usage = "Kludgy tool to simplify stuck Dot[] expressions in numerator. 
Expressions with sums of momentum tend not to be distributed correctly in previous versions of code (i.e. (p-q).(p-q)). 
This expands and simplifies those expressions (but only in the numerator)."


Begin["`Private`"]

Tarce/:Tarce[x_+y_,m__]:=Tarce[x,m]+Tarce[y,m]
Tarce/:Tarce[B___*c_*A___,p_,q_,k_,m__]/;FreeQ[c,q]&&FreeQ[c,k]:=c*Tarce[B*A,p,q,k,m]
Tarce/:Tarce[B___*c_*A__,p_,q_,k_,m__]/;FreeQ[c,q]&&FreeQ[c,k]:=c*Tarce[B*A,p,q,k,m]
Tarce/:Tarce[B___*p_. p_*A__,p_,q_,k_,m__]:=(p.p)*Tarce[B*A,p,q,k,m]
Tarce/:Tarce[A_/c_,p_,q_,k_,m__]/;FreeQ[c,q]&&FreeQ[c,k]:=(1/c)*Tarce[A,p,q,k,m]

Tarce/:Tarce[x_,p_,k1_,k2_,m1_,m2_, m3_, m4_,m5_]:=
((\[Pi])^(d))(-1)^(Exponent[x,m1^2-k1.k1]+Exponent[x,m2^2-k2.k2]+Exponent[x,m3^2-(p-k1).(p-k1)]+Exponent[x,m3^2-(k1-p).(k1-p)]+Exponent[x,m4^2-(-k2+p).(-k2+p)]+Exponent[x,m4^2-(k2-p).(k2-p)]+Exponent[x,m5^2-(-k2+k1).(-k2+k1)]+Exponent[x,m5^2-(k2-k1).(k2-k1)])TFI[d,p.p,
{
If[Exponent[x,k1.k1]>0,Exponent[x,k1.k1],0],
If[Exponent[x,k2.k2]>0,Exponent[x,k2.k2],0],
If[Exponent[x,p.k1]>0,Exponent[x,p.k1],If[Exponent[x,k1.p]>0,Exponent[x,k1.p],0]],
If[Exponent[x,p.k2]>0,Exponent[x,p.k2],If[Exponent[x,k2.p]>0,Exponent[x,k2.p],0]],
If[Exponent[x,k2.k1]>0,Exponent[x,k2.k1],If[Exponent[x,k1.k2]>0,Exponent[x,k1.k2],0]]
},
{
{If[Exponent[x,k1.k1-m1^2]<0,-Exponent[x,k1.k1-m1^2],If[Exponent[x,m1^2-k1.k1]<0,-Exponent[x,m1^2-k1.k1],0]],m1},
{If[Exponent[x,k2.k2-m2^2]<0,-Exponent[x,k2.k2-m2^2],If[Exponent[x,m2^2-k2.k2]<0,-Exponent[x,m2^2-k2.k2],0]],m2},
{If[Exponent[x,(p-k1).(p-k1)-m3^2]<0,-Exponent[x,(p-k1).(p-k1)-m3^2],
If[Exponent[x,(k1-p).(k1-p)-m3^2]<0,-Exponent[x,(k1-p).(k1-p)-m3^2],
If[Exponent[x,m3^2-(p-k1).(p-k1)]<0,-Exponent[x,m3^2-(p-k1).(p-k1)],
If[Exponent[x,m3^2-(k1-p).(k1-p)]<0,-Exponent[x,m3^2-(k1-p).(k1-p)],0]]]],m3},
{If[Exponent[x,(-k2+p).(-k2+p)-m4^2]<0,-Exponent[x,(-k2+p).(-k2+p)-m4^2],
If[Exponent[x,(k2-p).(k2-p)-m4^2]<0,-Exponent[x,(k2-p).(k2-p)-m4^2],
If[Exponent[x,m4^2-(-k2+p).(-k2+p)]<0,-Exponent[x,m4^2-(-k2+p).(-k2+p)],
If[Exponent[x,m4^2-(k2-p).(k2-p)]<0,-Exponent[x,m4^2-(k2-p).(k2-p)],0]]]],m4},
{If[Exponent[x,(-k2+k1).(-k2+k1)-m5^2]<0,-Exponent[x,(-k2+k1).(-k2+k1)-m5^2],
If[Exponent[x,(k2-k1).(k2-k1)-m5^2]<0,-Exponent[x,(k2-k1).(k2-k1)-m5^2],
If[Exponent[x,m5^2-(-k2+k1).(-k2+k1)]<0,-Exponent[x,m5^2-(-k2+k1).(-k2+k1)],
If[Exponent[x,m5^2-(k2-k1).(k2-k1)]<0,-Exponent[x,m5^2-(k2-k1).(k2-k1)],0]]]],m5}
}
]


numeratorReplace/:numeratorReplace[A_+B_,C__]:=numeratorReplace[A,C]+numeratorReplace[B,C]
(*numeratorReplace[A_,p_,q_,k_]:=If[MemberQ[Numerator[A],(p-q).(p-q)|(q-p).(q-p),Infinity],(Numerator[A]/.(p-q).(p-q)|(q-p).(q-p)\[Rule]p.p+q.q-2q.p)/Denominator[A],A]*)
numeratorReplace[A_,p_,q_,k_]:=Which[MemberQ[Numerator[A],(p-q).(p-q)|(q-p).(q-p),Infinity],(Numerator[A]/.(p-q).(p-q)|(q-p).(q-p)->p.p+q.q-2q.p)/Denominator[A],MemberQ[Numerator[A],(p+q).(p+q)|(q+p).(q+p),Infinity],(Numerator[A]/.(p+q).(p+q)|(q+p).(q+p)->p.p+q.q+2q.p)/Denominator[A],True,A]


End[]
Protect @@ Names["tarcer`*"];
EndPackage[]
