(* ::Package:: *)

BeginPackage["QCDParameters`"]
\[Alpha]GG::usage="Dimension 4 gluon condensate. Defined as ... (citation)"
mqq::usage="Dimension 4 light quark condensate. Defined as ... (citation)"
mss::usage="Dimension 4 strange quark condensate. Defined as ... (citation)"
qq::usage="Dimension 3 light quark condensate. Defined as ... (citation)"
ss::usage="Dimension 3 strange quark condensate. Defined as ... (citation)"
\[Psi]\[Sigma]G\[Psi]::usage="Dimension 5 'mixed' condensate. Defined as ... (citation)"
gggGGG::usage="Dimension 6 gluon condensate. Defined as ... (citation)"
qqqq::usage="Dimension 6 quark condensate. Defined as ... (citation)"

\[Alpha][\[Mu]]::usage="Running coupling with renormalization scale \[Mu]. Defined as ... (citation)"
mq[\[Mu]]::usage="Running light quark mass with renormalization scale \[Mu]. Defined as ... (citation)"
ms[\[Mu]]::usage="Running strange quark mass with renormalization scale \[Mu]. Defined as ... (citation)"
mc[\[Mu]]::usage="Running charm quark mass with renormalization scale \[Mu]. Defined as ... (citation)"
mb[\[Mu]]::usage="Running bottom quark mass with renormalization scale \[Mu]. Defined as ... (citation)"

Begin["`Private`"]
(*Constant Parameters*)

(*Bare Light Quark Masses @2GeV- PDG 2016*)
lightBare = 0.0035(*GeV*)(*PDG average - (u+d)/2*)
\[delta]lightBare = {+0.0007,-0.0003}
strangeBare = 0.096 (*GeV*)
\[delta]strangeBare = {+0.008,-0.004}

(*Bare Heavy Quark Masses @ Bare Mass - PDG 2016*)
charmBare = 1.27 (*GeV*)
\[delta]charmBare = {0.03,-0.03}
bottomBare = 4.18 (*GeV*)
\[delta]bottomBare = {0.04,-0.03}

(*Quark Mass Ratios*)
strangelightRatio = 0.0273 (*PDGLive 2016*)
\[delta]strangelightRatio  = {0.0007,-0.0007} (*PDGLive 2016*)
charmlightRatio = 322.6 (*Calculated in J.Ho, D. Harnett, T.G. Steele 2017*)
\[delta]charmlightRatio = {13.6,-13.6} (*Calculated in J.Ho, D. Harnett, T.G. Steele 2017*)
charmstrangeRatio = 11.72 (*PDGLive 2016*)
\[delta]charmstrangeRatio = {0.25,-0.25} (*PDGLive 2016*)
bottomlightRatio = 1460.7 (*Calculated in J.Ho, D. Harnett, T.G. Steele 2017*)
\[delta]bottomlightRatio = {64,-64} (*Calculated in J.Ho, D. Harnett, T.G. Steele 2017*)
bottomstrangeRatio = 52.55 (*HPQCD Collab, PRD91 054508*)
\[delta]bottomstrangeRatio = {1.30,-1.30} (*HPQCD Collab, PRD91 054508*)
bottomcharmRatio = 4.53 (*HPQCD Collab, PRD91 054508*)
\[delta]bottomcharmRatio = {0.05,-0.05} (*HPQCD Collab, PRD91 054508*)

(*Condensate Parameters - Sources as stated*)
\[Alpha]GG = 0.07 (*GeV^4; arXiv:1511.05903v3*)
\[delta]\[Alpha]GG = {0.020,-0.020}

Subscript[c,6] = 8.2 (*GeV^2*)
\[delta]Subscript[c,6] = {1,-1}
gggGGG = Subscript[c,6] \[Alpha]GG (*arXiv:1511.05903v3*)

Subscript[m,Pi] = 0.13957018 (*GeV*)
\[delta]Subscript[m,Pi] = {0.00000035,-0.00000035}
Subscript[m,k] = 0.493677 (*GeV*)
\[delta]Subscript[m,k] = {0.000013,-0.000013}
Subscript[f,Pi] = 0.0922 (*arXiv:1509.02220 [hep-ph]*)
\[delta]Subscript[f,Pi] = {0.0035,-0.0035} (*arXiv:1509.02220 [hep-ph]*)
Subscript[f,k]= 110.0(*arXiv:1509.02220 [hep-ph]*)
\[delta]Subscript[f,k] = 110.0(*arXiv:1509.02220 [hep-ph]*)
mqq = (-(1/2)Subscript[f,Pi]^2 Subscript[m,Pi]^2)
mss = (-(1/2)Subscript[f,k]^2 Subscript[m,k]^2)

M0squared = 0.8 (*GeV^2; arXiv:1511.05903v3*)
\[delta]M0squared = {0.1,-0.1}
q\[Sigma]Gq = M0squared(charmlightRatio)mqq
s\[Sigma]Gs =  M0squared(charmstrangeRatio)mss

\[Kappa] = 0.74 (*arXiv:1511.05903v3*)
\[Alpha]qqqq = 3.5*10^-4 (*Latorre, Pascual, Narison, Z.Phys.C34 (1987) 347*)
\[Alpha]ssss = \[Kappa]^2 \[Alpha]qqqq (*Latorre, Pascual, Narison, Z.Phys.C34 (1987) 347*)

(*Running Coupling*)
M\[Tau] = 1.77699 (*GeV*)
\[Alpha]\[Tau] = 0.326
\[delta]\[Alpha]\[Tau] = {0.014,-0.014}
\[Alpha]Z = 0.1185
\[delta]\[Alpha]Z = {0.0006,-0.0006}
\[Alpha][\[Mu]_] = Module[
                          {nf = 4, b0},
                            b0 = 1/(12 \[Pi]) (33 - 2 nf);
                            \[Alpha]\[Tau]/(1 + b0 \[Alpha]\[Tau] Log[\[Mu]^2/M\[Tau]^2])
                        ];
(*Running Mass*)
mq[\[Mu]_]:= lightBare (\[Alpha][\[Mu]/\[Alpha][2(*GeV*)]])^(4/9)
ms[\[Mu]_]:= strangeBare (\[Alpha][\[Mu]/\[Alpha][2(*GeV*)]])^(4/9)
mc[\[Mu]_]:= charmBare (\[Alpha][\[Mu]/\[Alpha][charmBare]])^(12/25)
mb[\[Mu]_]:= bottomBare (\[Alpha][\[Mu]/\[Alpha][charmBare]])^(12/25)
End[]
EndPackage[]
