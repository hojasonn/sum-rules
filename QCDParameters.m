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
strangelightRatio = 0.0273
\[delta]strangelightRatio  = {0.0007,-0.0007}
charmstrangeRatio = 11.72
\[delta]charmstrangeRatio = {0.25,-0.25}
bottomcharmRatio = 4.53
\[delta]bottomcharmRatio = {0.05,-0.05}

(*Condensate Parameters - Sources as stated*)


(*Running Coupling*)
\[Alpha][\[Mu]_] = Module[
                          {M\[Tau] = 1.77699, \[Alpha]\[Tau] = 0.330, nf = 4, b0},
                            b0 = 1/(12 \[Pi]) (33 - 2 nf);
                            \[Alpha]\[Tau]/(1 + b0 \[Alpha]\[Tau] Log[\[Mu]^2/M\[Tau]^2])
                        ];
(*Running Mass*)
mq[\[Mu]_]:= lightBare(\[Alpha][\[Mu]/\[Alpha][2(*GeV*)]])^(4/9)
ms[\[Mu]_]:= strangeBare(\[Alpha][\[Mu]/\[Alpha][2(*GeV*)]])^(4/9)
mc[\[Mu]_]:= charmBare(\[Alpha][\[Mu]/\[Alpha][charmBare]])^(12/25)
mb[\[Mu]_]:= bottomBare(\[Alpha][\[Mu]/\[Alpha][charmBare]])^(12/25)
End[]
EndPackage[]
