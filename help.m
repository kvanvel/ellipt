(* ::Package:: *)

ClearAll["Global`*"]


ClearAll["Global`*"]
SetDirectory["~/Workspace/ellipt/src/step-1"];

massU = Import["MassU.mtx"];
massQ = Import["MassQ.mtx"];
inverseMassU = Import["InverseMassU.mtx"];
inverseMassQ = Import["InverseMassQ.mtx"];
traceUdir = Import["TraceUDir.mtx"];
traceQNeu = Import["TraceQNeu.mtx"]; 
totalUFromQ=Import["TotalUFromQ.mtx"];
totalQFromU=Import["TotalQFromU.mtx"];
(*fluxBoundaryUFromQ=Import["FluxBoundaryUFromQ.mtx"];
fluxBoundaryQFromU=Import["FluxBoundaryQFromU.mtx"];
massElec=Import["MassElec.mtx"];*)





massU//MatrixPlot
stiff=totalQFromU;
stiff//MatrixPlot


mmQ=massQ[[1;;8,1;;8]];
mmU=massU[[1;;4,1;;4]]


mm//MatrixPlot


ss=stiff[[1;;8,1;;4]]


ss//MatrixPlot


sss=MatrixPower[Normal[mmQ],-1/2].ss.MatrixPower[Normal[mmU],-1/2]//Chop


sss//MatrixPlot


sss.Transpose[sss]


PseudoInverse[sss]-Transpose[sss].PseudoInverse[sss.Transpose[sss]] 


NullSpace[Transpose[sss]]


PseudoInverse[sss]





\[CapitalXi]//Normal


?traceUdir


stiff = Import["StiffQFromU.mtx"];

S=totalQFromU;
S=ConjugateTranspose[S];
minusSstar= totalUFromQ;
R=ConjugateTranspose[totalUFromQ];
\[CapitalGamma]=traceUdir;
\[CapitalXi]=traceQNeu;
(*
Needs["GraphUtilities`"]
{r,c} = MinimumBandwidthOrdering[\[CapitalGamma]];

permutationMatrix[p_List]:=SparseArray@IdentityMatrix[Length[p]][[p]]

zz=SparseArray[Transpose[SparseArray[Chop[NullSpace[\[CapitalGamma][[r,c]]]]].permutationMatrix[c]]];

\[Omega]=1/10;
ZZ=SparseArray[ArrayFlatten[{{-I \[Omega] massU,S,ConjugateTranspose[\[CapitalGamma]],0},{-ConjugateTranspose[S],(-I \[Omega]+1) massQ,0,ConjugateTranspose[\[CapitalXi]]},{\[CapitalGamma],0,0,0},{0,\[CapitalXi],0,0}}]];
rhs=1.0RandomInteger[{-10,10},Length[ZZ]]+I 1.0RandomInteger[{-10,10},Length[ZZ]];

soln=LinearSolve[ZZ,rhs];

uSoln=soln[[1;;Length[massU]]];
uR=Re/@ uSoln;
uI=Im/@uSoln;
qSoln=soln[[Length[massU]+1;;Length[massU]+Length[massQ]]];
qR=Re/@qSoln;
qI=Im/@qSoln;
\[Lambda]Soln=soln[[Length[massU]+1+Length[massQ];;Length[massU]+Length[massQ]+Length[\[CapitalGamma]]]];

\[Lambda]R=Re/@\[Lambda]Soln;
\[Lambda]I=Im/@\[Lambda]Soln;
\[Mu]Soln=
soln[[
Length[massU]+1+Length[massQ]+Length[\[CapitalGamma]];;Length[massU]+Length[massQ]+Length[\[CapitalGamma]]+Length[\[CapitalXi]]
]];
\[Mu]R=Re/@\[Mu]Soln;
\[Mu]I=Im/@\[Mu]Soln;
rhs1=rhs[[1;;Length[massU]]];
rhs2=rhs[[Length[massU]+1;;Length[massU]+Length[massQ]]];
rhs3=rhs[[Length[massU]+Length[massQ]+1;;Length[massU]+Length[massQ]+Length[\[CapitalGamma]]]];
rhs4=rhs[[Length[massU]+Length[massQ]+Length[\[CapitalGamma]]+1;;-1]];

rhs1R=Re/@rhs1;
rhs1I=Im/@rhs1;
rhs2R=Re/@rhs2;
rhs2I=Im/@rhs2;
rhs3R=Re/@rhs3;
rhs3I=Im/@rhs3;
rhs4R=Re/@rhs4;
rhs4I=Im/@rhs4;*)


?\[CapitalGamma]


 


V=Transpose[NullSpace[\[CapitalGamma]]];
W=Transpose[NullSpace[\[CapitalXi]]];

\[ScriptL]=LinearSolve[ConjugateTranspose[V].massU.V,ConjugateTranspose[V].massU.uSoln];
\[Alpha]=LinearSolve[\[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]],\[CapitalGamma].uSoln];
\[ScriptM]=LinearSolve[ConjugateTranspose[W].massQ.W,ConjugateTranspose[W].massQ.qSoln];
\[Beta]=LinearSolve[\[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]],\[CapitalXi].qSoln];
\[ScriptL]R=Re/@ \[ScriptL];
\[ScriptL]I=Im/@ \[ScriptL];
\[Alpha]R= Re/@ \[Alpha];
\[Alpha]I= Im/@ \[Alpha];
\[ScriptM]R=Re/@ \[ScriptM];
\[ScriptM]I=Im/@\[ScriptM];
\[Beta]R= Re/@ \[Beta];
\[Beta]I=Im/@ \[Beta];


Z=ArrayFlatten[
  {
   {-I \[Omega] massU,S,ConjugateTranspose[\[CapitalGamma]],0},
   {-ConjugateTranspose[S],-I \[Omega] massQ+massQ,0,ConjugateTranspose[\[CapitalXi]]},
   {\[CapitalGamma],0,0,0},{0,\[CapitalXi],0,0}
   }
];
ZZ=ArrayFlatten[{
 {-I \[Omega] massU.V, -I \[Omega] massU.inverseMassU.ConjugateTranspose[\[CapitalGamma]], S.W, S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], ConjugateTranspose[\[CapitalGamma]], 0},
 {-ConjugateTranspose[S].V, -ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], -I \[Omega] massQ.W+massQ.W, -I \[Omega] massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]]+massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, ConjugateTranspose[\[CapitalXi]]},
 {\[CapitalGamma].V, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, 0},
 {0, 0, \[CapitalXi].W, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0}
}];


ZZZ=ArrayFlatten[{
 {-I \[Omega] ConjugateTranspose[V].massU.V, -I \[Omega] ConjugateTranspose[V].massU.inverseMassU.ConjugateTranspose[\[CapitalGamma]], ConjugateTranspose[V].S.W, ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], ConjugateTranspose[V].ConjugateTranspose[\[CapitalGamma]], 0},
 {-I \[Omega] \[CapitalGamma].inverseMassU.massU.V, -I \[Omega] \[CapitalGamma].inverseMassU.massU.inverseMassU.ConjugateTranspose[\[CapitalGamma]], \[CapitalGamma].inverseMassU.S.W, \[CapitalGamma].inverseMassU.S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0},
 {-\[CapitalXi].inverseMassQ.ConjugateTranspose[S].V, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], -I \[Omega] \[CapitalXi].inverseMassQ.massQ.W+\[CapitalXi].inverseMassQ.massQ.W, -I \[Omega] \[CapitalXi].inverseMassQ.massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]]+\[CapitalXi].inverseMassQ.massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]]},
 {-ConjugateTranspose[W].ConjugateTranspose[S].V, -ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], -I \[Omega] ConjugateTranspose[W].massQ.W+ConjugateTranspose[W].massQ.W, -I \[Omega] ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]]+ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, ConjugateTranspose[W].ConjugateTranspose[\[CapitalXi]]},
 {\[CapitalGamma].V, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, 0},
 {0, 0, \[CapitalXi].W, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0}
}];


ZZZZ=ArrayFlatten[{
 {-I \[Omega] ConjugateTranspose[V].massU.V, 0, ConjugateTranspose[V].S.W, ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0},
 {0, -I \[Omega] \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], \[CapitalGamma].inverseMassU.S.W, \[CapitalGamma].inverseMassU.S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0},
 {-ConjugateTranspose[W].ConjugateTranspose[S].V, -ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], -I \[Omega] ConjugateTranspose[W].massQ.W+ConjugateTranspose[W].massQ.W, ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0},
 {-\[CapitalXi].inverseMassQ.ConjugateTranspose[S].V, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], \[CapitalXi].inverseMassQ.massQ.W, -I \[Omega] \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]]+\[CapitalXi].inverseMassQ.massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]]},
 {0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, 0},
 {0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0}
}];



BBBB=Join[ConjugateTranspose[V].rhs1,\[CapitalGamma].inverseMassU.rhs1,ConjugateTranspose[W].rhs2,\[CapitalXi].inverseMassQ.rhs2,rhs3,rhs4];


ZZZZZ=ArrayFlatten[{
 {-\[Omega] ConjugateTranspose[V].massU.V, 0, 0, 0, 0, ConjugateTranspose[V].S.W, 0, ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, 0, 0},
 {0, \[Omega] ConjugateTranspose[V].massU.V, 0, 0, ConjugateTranspose[V].S.W, 0, ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, 0, 0, 0},
 {0, 0, -\[Omega] \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, \[CapitalGamma].inverseMassU.S.W, 0, \[CapitalGamma].inverseMassU.S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0},
 {0, 0, 0, \[Omega] \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], \[CapitalGamma].inverseMassU.S.W, 0, \[CapitalGamma].inverseMassU.S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0},
 {-ConjugateTranspose[W].ConjugateTranspose[S].V, 0, -ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, ConjugateTranspose[W].massQ.W, \[Omega] ConjugateTranspose[W].massQ.W, ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, 0, 0, 0},
 {0, -ConjugateTranspose[W].ConjugateTranspose[S].V, 0, -ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], -\[Omega] ConjugateTranspose[W].massQ.W, ConjugateTranspose[W].massQ.W, 0, ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, 0, 0},
 {-\[CapitalXi].inverseMassQ.ConjugateTranspose[S].V, 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, \[CapitalXi].inverseMassQ.massQ.W, 0, \[CapitalXi].inverseMassQ.massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], \[Omega] \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0},
 {0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].V, 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, \[CapitalXi].inverseMassQ.massQ.W, -\[Omega] \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], \[CapitalXi].inverseMassQ.massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]]},
 {0, 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, 0, 0, 0},
 {0, 0, 0, 0, 0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, 0, 0}
}];

BBBBB=
Join[
ConjugateTranspose[V].rhs1I,
ConjugateTranspose[V].rhs1R,
\[CapitalGamma].inverseMassU.rhs1I,
\[CapitalGamma].inverseMassU.rhs1R,
ConjugateTranspose[W].rhs2R,
ConjugateTranspose[W].rhs2I,
\[CapitalXi].inverseMassQ.rhs2R,
\[CapitalXi].inverseMassQ.rhs2I,
rhs3R,
rhs3I,
rhs4R,
rhs4I];

ZZZZZ.Join[\[ScriptL]R,\[ScriptL]I,\[Alpha]R,\[Alpha]I,\[ScriptM]R,\[ScriptM]I,\[Beta]R,\[Beta]I,\[Lambda]R,\[Lambda]I,\[Mu]R,\[Mu]I]-BBBBB//Chop//Norm



sys=SparseArray@ArrayFlatten[{
 {-\[Omega] ConjugateTranspose[V].massU.V, 0, 0, 0, 0, ConjugateTranspose[V].S.W, 0, ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, 0, 0, Transpose[{ConjugateTranspose[V].rhs1I}]},
 {0, \[Omega] ConjugateTranspose[V].massU.V, 0, 0, ConjugateTranspose[V].S.W, 0, ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, 0, 0, 0, Transpose[{ConjugateTranspose[V].rhs1R}]},
 {0, 0, -\[Omega] \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, \[CapitalGamma].inverseMassU.S.W, 0, \[CapitalGamma].inverseMassU.S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, Transpose[{\[CapitalGamma].inverseMassU.rhs1I}]},
 {0, 0, 0, \[Omega] \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], \[CapitalGamma].inverseMassU.S.W, 0, \[CapitalGamma].inverseMassU.S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, Transpose[{\[CapitalGamma].inverseMassU.rhs1R}]},
 {-ConjugateTranspose[W].ConjugateTranspose[S].V, 0, -ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, ConjugateTranspose[W].massQ.W, \[Omega] ConjugateTranspose[W].massQ.W, ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, 0, 0, 0, Transpose[{ConjugateTranspose[W].rhs2R}]},
 {0, -ConjugateTranspose[W].ConjugateTranspose[S].V, 0, -ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], -\[Omega] ConjugateTranspose[W].massQ.W, ConjugateTranspose[W].massQ.W, 0, ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, 0, 0, Transpose[{ConjugateTranspose[W].rhs2I}]},
 {-\[CapitalXi].inverseMassQ.ConjugateTranspose[S].V, 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, \[CapitalXi].inverseMassQ.massQ.W, 0, \[CapitalXi].inverseMassQ.massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], \[Omega] \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, Transpose[{\[CapitalXi].inverseMassQ.rhs2R}]},
 {0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].V, 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, \[CapitalXi].inverseMassQ.massQ.W, -\[Omega] \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], \[CapitalXi].inverseMassQ.massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], Transpose[{\[CapitalXi].inverseMassQ.rhs2I}]},
 {0, 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, 0, 0, 0, 0, 0, 0, Transpose[{rhs3R}]},
 {0, 0, 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, 0, 0, 0, 0, 0, Transpose[{rhs3I}]},
 {0, 0, 0, 0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, 0, 0, 0, Transpose[{rhs4R}]},
 {0, 0, 0, 0, 0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, 0, 0, Transpose[{rhs4I}]}
}];
Map[Last,RowReduce[sys]]-Join[\[ScriptL]R,\[ScriptL]I,\[Alpha]R,\[Alpha]I,\[ScriptM]R,\[ScriptM]I,\[Beta]R,\[Beta]I,\[Lambda]R,\[Lambda]I,\[Mu]R,\[Mu]I]//Chop//Norm


\[Omega]


sys=SparseArray@ArrayFlatten[{
 { ConjugateTranspose[V].massU.V, 0, 0, 0, 0, -(1/\[Omega])ConjugateTranspose[V].S.W, 0, -1/\[Omega] ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, 0, 0, -1/\[Omega] Transpose[{ConjugateTranspose[V].rhs1I}]},
 {0,  ConjugateTranspose[V].massU.V, 0, 0, 1/\[Omega] ConjugateTranspose[V].S.W, 0, 1/\[Omega] ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, 0, 0, 0, 1/\[Omega] Transpose[{ConjugateTranspose[V].rhs1R}]},
 {0, 0,  \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, -1/\[Omega] \[CapitalGamma].inverseMassU.S.W, 0, -1/\[Omega] \[CapitalGamma].inverseMassU.S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, -1/\[Omega] \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, -1/\[Omega] Transpose[{\[CapitalGamma].inverseMassU.rhs1I}]},
 {0, 0, 0,  \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 1/\[Omega] \[CapitalGamma].inverseMassU.S.W, 0, 1/\[Omega] \[CapitalGamma].inverseMassU.S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 1/\[Omega] \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, 1/\[Omega] Transpose[{\[CapitalGamma].inverseMassU.rhs1R}]},
 {-ConjugateTranspose[W].ConjugateTranspose[S].V, 0, -ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, ConjugateTranspose[W].massQ.W, \[Omega] ConjugateTranspose[W].massQ.W, ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, 0, 0, 0, Transpose[{ConjugateTranspose[W].rhs2R}]},
 {0, -ConjugateTranspose[W].ConjugateTranspose[S].V, 0, -ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], -\[Omega] ConjugateTranspose[W].massQ.W, ConjugateTranspose[W].massQ.W, 0, ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, 0, 0, Transpose[{ConjugateTranspose[W].rhs2I}]},
 {-\[CapitalXi].inverseMassQ.ConjugateTranspose[S].V, 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, \[CapitalXi].inverseMassQ.massQ.W, 0, \[CapitalXi].inverseMassQ.massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], \[Omega] \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, Transpose[{\[CapitalXi].inverseMassQ.rhs2R}]},
 {0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].V, 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, \[CapitalXi].inverseMassQ.massQ.W, -\[Omega] \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], \[CapitalXi].inverseMassQ.massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], Transpose[{\[CapitalXi].inverseMassQ.rhs2I}]},
 {0, 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, 0, 0, 0, 0, 0, 0, Transpose[{rhs3R}]},
 {0, 0, 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, 0, 0, 0, 0, 0, Transpose[{rhs3I}]},
 {0, 0, 0, 0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, 0, 0, 0, Transpose[{rhs4R}]},
 {0, 0, 0, 0, 0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, 0, 0, Transpose[{rhs4I}]}
}];
Map[Last,RowReduce[sys]]-Join[\[ScriptL]R,\[ScriptL]I,\[Alpha]R,\[Alpha]I,\[ScriptM]R,\[ScriptM]I,\[Beta]R,\[Beta]I,\[Lambda]R,\[Lambda]I,\[Mu]R,\[Mu]I]//Chop//Norm


sys=SparseArray@ArrayFlatten[{
 { ConjugateTranspose[V].massU.V, 0, 0, -(1/\[Omega])ConjugateTranspose[V].S.W, 0, 0, 0, 0, 0, 0, 0, -1/\[Omega] ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], -1/\[Omega] Transpose[{ConjugateTranspose[V].rhs1I}]},
 {0,  ConjugateTranspose[V].massU.V, 1/\[Omega] ConjugateTranspose[V].S.W, 0, 0, 0, 0, 0, 0, 0, 1/\[Omega] ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 1/\[Omega] Transpose[{ConjugateTranspose[V].rhs1R}]},
 {0, 0, 0, -1/\[Omega] \[CapitalGamma].inverseMassU.S.W, 0, -1/\[Omega] \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0,  \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, -1/\[Omega] \[CapitalGamma].inverseMassU.S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], -1/\[Omega] Transpose[{\[CapitalGamma].inverseMassU.rhs1I}]},
 {0, 0, 1/\[Omega] \[CapitalGamma].inverseMassU.S.W, 0, 1/\[Omega] \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, 0,  \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 1/\[Omega] \[CapitalGamma].inverseMassU.S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 1/\[Omega] Transpose[{\[CapitalGamma].inverseMassU.rhs1R}]},
 {-ConjugateTranspose[W].ConjugateTranspose[S].V, 0, ConjugateTranspose[W].massQ.W, \[Omega] ConjugateTranspose[W].massQ.W, 0, 0, 0, 0, -ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, Transpose[{ConjugateTranspose[W].rhs2R}]},
 {0, -ConjugateTranspose[W].ConjugateTranspose[S].V, -\[Omega] ConjugateTranspose[W].massQ.W, ConjugateTranspose[W].massQ.W, 0, 0, 0, 0, 0, -ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], Transpose[{ConjugateTranspose[W].rhs2I}]},
 {-\[CapitalXi].inverseMassQ.ConjugateTranspose[S].V, 0, \[CapitalXi].inverseMassQ.massQ.W, 0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, \[CapitalXi].inverseMassQ.massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], \[Omega] \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], Transpose[{\[CapitalXi].inverseMassQ.rhs2R}]},
 {0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].V, 0, \[CapitalXi].inverseMassQ.massQ.W, 0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], -\[Omega] \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], \[CapitalXi].inverseMassQ.massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], Transpose[{\[CapitalXi].inverseMassQ.rhs2I}]},
 {0, 0, 0, 0, 0, 0, 0, 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, Transpose[{rhs3R}]},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, Transpose[{rhs3I}]},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, Transpose[{rhs4R}]},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], Transpose[{rhs4I}]}
}];
Map[Last,RowReduce[sys]]-Join[\[ScriptL]R,\[ScriptL]I,\[ScriptM]R,\[ScriptM]I,\[Lambda]R,\[Lambda]I,\[Mu]R,\[Mu]I,\[Alpha]R,\[Alpha]I,\[Beta]R,\[Beta]I]//Chop//Norm


sys=SparseArray@ArrayFlatten[{
 { ConjugateTranspose[V].massU.V, 0, 0, -(1/\[Omega])ConjugateTranspose[V].S.W, 0, 0, 0, 0, 0, 0, 0, -1/\[Omega] ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], -1/\[Omega] Transpose[{ConjugateTranspose[V].rhs1I}]},
 {0,  ConjugateTranspose[V].massU.V, 1/\[Omega] ConjugateTranspose[V].S.W, 0, 0, 0, 0, 0, 0, 0, 1/\[Omega] ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 1/\[Omega] Transpose[{ConjugateTranspose[V].rhs1R}]},
 {0, 0, 0, -1/\[Omega] \[CapitalGamma].inverseMassU.S.W, 0, -1/\[Omega] \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0,  \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, -1/\[Omega] \[CapitalGamma].inverseMassU.S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], -1/\[Omega] Transpose[{\[CapitalGamma].inverseMassU.rhs1I}]},
 {0, 0, 1/\[Omega] \[CapitalGamma].inverseMassU.S.W, 0, 1/\[Omega] \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, 0,  \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 1/\[Omega] \[CapitalGamma].inverseMassU.S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 1/\[Omega] Transpose[{\[CapitalGamma].inverseMassU.rhs1R}]},
 {-ConjugateTranspose[W].ConjugateTranspose[S].V, 0, ConjugateTranspose[W].massQ.W, \[Omega] ConjugateTranspose[W].massQ.W, 0, 0, 0, 0, -ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, Transpose[{ConjugateTranspose[W].rhs2R}]},
 {0, -ConjugateTranspose[W].ConjugateTranspose[S].V, -\[Omega] ConjugateTranspose[W].massQ.W, ConjugateTranspose[W].massQ.W, 0, 0, 0, 0, 0, -ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], Transpose[{ConjugateTranspose[W].rhs2I}]},
 {-\[CapitalXi].inverseMassQ.ConjugateTranspose[S].V, 0, \[CapitalXi].inverseMassQ.massQ.W, 0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, \[CapitalXi].inverseMassQ.massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], \[Omega] \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], Transpose[{\[CapitalXi].inverseMassQ.rhs2R}]},
 {0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].V, 0, \[CapitalXi].inverseMassQ.massQ.W, 0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], -\[Omega] \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], \[CapitalXi].inverseMassQ.massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], Transpose[{\[CapitalXi].inverseMassQ.rhs2I}]},
 {0, 0, 0, 0, 0, 0, 0, 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, Transpose[{rhs3R}]},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, Transpose[{rhs3I}]},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, Transpose[{rhs4R}]},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], Transpose[{rhs4I}]}
}];
Map[Last,RowReduce[sys]]-Join[\[ScriptL]R,\[ScriptL]I,\[ScriptM]R,\[ScriptM]I,\[Lambda]R,\[Lambda]I,\[Mu]R,\[Mu]I,\[Alpha]R,\[Alpha]I,\[Beta]R,\[Beta]I]//Chop//Norm

sys=SparseArray@ArrayFlatten[{
 {0, -1/\[Omega] \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, 0, 0, -1/\[Omega] \[CapitalGamma].inverseMassU.S.W,  \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, -1/\[Omega] \[CapitalGamma].inverseMassU.S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], -1/\[Omega] Transpose[{\[CapitalGamma].inverseMassU.rhs1I}]},
 {1/\[Omega] \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, 0, 0, 1/\[Omega] \[CapitalGamma].inverseMassU.S.W, 0, 0,  \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 1/\[Omega] \[CapitalGamma].inverseMassU.S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 1/\[Omega] Transpose[{\[CapitalGamma].inverseMassU.rhs1R}]},
 {0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].V, 0, \[CapitalXi].inverseMassQ.massQ.W, 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, \[CapitalXi].inverseMassQ.massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], \[Omega] \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], Transpose[{\[CapitalXi].inverseMassQ.rhs2R}]},
 {0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].V, 0, \[CapitalXi].inverseMassQ.massQ.W, 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], -\[Omega] \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], \[CapitalXi].inverseMassQ.massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], Transpose[{\[CapitalXi].inverseMassQ.rhs2I}]},
 {0, 0, 0, 0,  ConjugateTranspose[V].massU.V, 0, 0, -(1/\[Omega])ConjugateTranspose[V].S.W, 0, 0, 0, -1/\[Omega] ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], -1/\[Omega] Transpose[{ConjugateTranspose[V].rhs1I}]},
 {0, 0, 0, 0, 0,  ConjugateTranspose[V].massU.V, 1/\[Omega] ConjugateTranspose[V].S.W, 0, 0, 0, 1/\[Omega] ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 1/\[Omega] Transpose[{ConjugateTranspose[V].rhs1R}]},
 {0, 0, 0, 0, 0, 0, ConjugateTranspose[W].massQ.W, -(1/\[Omega])ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W+\[Omega] ConjugateTranspose[W].massQ.W, -ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], -1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], -1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].Transpose[{ConjugateTranspose[V].rhs1I}]+Transpose[{ConjugateTranspose[W].rhs2R}]},
 {0, 0, 0, 0, 0, -ConjugateTranspose[W].ConjugateTranspose[S].V, -\[Omega] ConjugateTranspose[W].massQ.W, ConjugateTranspose[W].massQ.W, 0, -ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], Transpose[{ConjugateTranspose[W].rhs2I}]},
 {0, 0, 0, 0, 0, 0, 0, 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, Transpose[{rhs3R}]},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, Transpose[{rhs3I}]},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, Transpose[{rhs4R}]},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], Transpose[{rhs4I}]}
}];
Map[Last,RowReduce[sys]]-Join[\[Lambda]R,\[Lambda]I,\[Mu]R,\[Mu]I,\[ScriptL]R,\[ScriptL]I,\[ScriptM]R,\[ScriptM]I,\[Alpha]R,\[Alpha]I,\[Beta]R,\[Beta]I]//Chop//Norm


sys=SparseArray@ArrayFlatten[{
 {0, -1/\[Omega] \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, 0, 0, -1/\[Omega] \[CapitalGamma].inverseMassU.S.W,  \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, -1/\[Omega] \[CapitalGamma].inverseMassU.S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], -1/\[Omega] Transpose[{\[CapitalGamma].inverseMassU.rhs1I}]},
 {1/\[Omega] \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, 0, 0, 1/\[Omega] \[CapitalGamma].inverseMassU.S.W, 0, 0,  \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 1/\[Omega] \[CapitalGamma].inverseMassU.S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 1/\[Omega] Transpose[{\[CapitalGamma].inverseMassU.rhs1R}]},
 {0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].V, 0, \[CapitalXi].inverseMassQ.massQ.W, 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, \[CapitalXi].inverseMassQ.massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], \[Omega] \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], Transpose[{\[CapitalXi].inverseMassQ.rhs2R}]},
 {0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].V, 0, \[CapitalXi].inverseMassQ.massQ.W, 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], -\[Omega] \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], \[CapitalXi].inverseMassQ.massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], Transpose[{\[CapitalXi].inverseMassQ.rhs2I}]},
 {0, 0, 0, 0,  ConjugateTranspose[V].massU.V, 0, 0, -(1/\[Omega])ConjugateTranspose[V].S.W, 0, 0, 0, -1/\[Omega] ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], -1/\[Omega] Transpose[{ConjugateTranspose[V].rhs1I}]},
 {0, 0, 0, 0, 0,  ConjugateTranspose[V].massU.V, 1/\[Omega] ConjugateTranspose[V].S.W, 0, 0, 0, 1/\[Omega] ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 1/\[Omega] Transpose[{ConjugateTranspose[V].rhs1R}]},
 {0, 0, 0, 0, 0, 0, ConjugateTranspose[W].massQ.W, -(1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega] ConjugateTranspose[W].massQ.W), -ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], -1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], -1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].Transpose[{ConjugateTranspose[V].rhs1I}]+Transpose[{ConjugateTranspose[W].rhs2R}]},
 {0, 0, 0, 0, 0, 0, (1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega] ConjugateTranspose[W].massQ.W), ConjugateTranspose[W].massQ.W, 0, -ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].Transpose[{ConjugateTranspose[V].rhs1R}]+Transpose[{ConjugateTranspose[W].rhs2I}]},
 {0, 0, 0, 0, 0, 0, 0, 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, Transpose[{rhs3R}]},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, Transpose[{rhs3I}]},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, Transpose[{rhs4R}]},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], Transpose[{rhs4I}]}
}];
Map[Last,RowReduce[sys]]-Join[\[Lambda]R,\[Lambda]I,\[Mu]R,\[Mu]I,\[ScriptL]R,\[ScriptL]I,\[ScriptM]R,\[ScriptM]I,\[Alpha]R,\[Alpha]I,\[Beta]R,\[Beta]I]//Chop//Norm


sys//MatrixPlot;
sys=SparseArray@ArrayFlatten[{
 {0, -1/\[Omega] \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, 0, 0, -1/\[Omega] \[CapitalGamma].inverseMassU.S.W,  \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, -1/\[Omega] \[CapitalGamma].inverseMassU.S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], -1/\[Omega] Transpose[{\[CapitalGamma].inverseMassU.rhs1I}]},
 {1/\[Omega] \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, 0, 0, 1/\[Omega] \[CapitalGamma].inverseMassU.S.W, 0, 0,  \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 1/\[Omega] \[CapitalGamma].inverseMassU.S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 1/\[Omega] Transpose[{\[CapitalGamma].inverseMassU.rhs1R}]},
 {0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].V, 0, \[CapitalXi].inverseMassQ.massQ.W, 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, \[CapitalXi].inverseMassQ.massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], \[Omega] \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], Transpose[{\[CapitalXi].inverseMassQ.rhs2R}]},
 {0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].V, 0, \[CapitalXi].inverseMassQ.massQ.W, 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], -\[Omega] \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], \[CapitalXi].inverseMassQ.massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], Transpose[{\[CapitalXi].inverseMassQ.rhs2I}]},
 {0, 0, 0, 0,  ConjugateTranspose[V].massU.V, 0, 0, -(1/\[Omega])ConjugateTranspose[V].S.W, 0, 0, 0, -1/\[Omega] ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], -1/\[Omega] Transpose[{ConjugateTranspose[V].rhs1I}]},
 {0, 0, 0, 0, 0,  ConjugateTranspose[V].massU.V, 1/\[Omega] ConjugateTranspose[V].S.W, 0, 0, 0, 1/\[Omega] ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 1/\[Omega] Transpose[{ConjugateTranspose[V].rhs1R}]},
 {0, 0, 0, 0, 0, 0, ConjugateTranspose[W].massQ.W, -(1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega] ConjugateTranspose[W].massQ.W), -ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], -1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], -1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].Transpose[{ConjugateTranspose[V].rhs1I}]+Transpose[{ConjugateTranspose[W].rhs2R}]},
 {0, 0, 0, 0, 0, 0, 0, -(1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega] ConjugateTranspose[W].massQ.W).Inverse[ConjugateTranspose[W].massQ.W].(1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega] ConjugateTranspose[W].massQ.W)-ConjugateTranspose[W].massQ.W, -(1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega] ConjugateTranspose[W].massQ.W).Inverse[ConjugateTranspose[W].massQ.W].ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], (1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega] ConjugateTranspose[W].massQ.W).Inverse[ConjugateTranspose[W].massQ.W].ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]]-1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], -1/\[Omega] (1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega] ConjugateTranspose[W].massQ.W).Inverse[ConjugateTranspose[W].massQ.W].ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]]-ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], (1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega] ConjugateTranspose[W].massQ.W).Inverse[ConjugateTranspose[W].massQ.W].(-1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].Transpose[{ConjugateTranspose[V].rhs1I}]+Transpose[{ConjugateTranspose[W].rhs2R}])-(1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].Transpose[{ConjugateTranspose[V].rhs1R}]+Transpose[{ConjugateTranspose[W].rhs2I}])},
 {0, 0, 0, 0, 0, 0, 0, 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, Transpose[{rhs3R}]},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, Transpose[{rhs3I}]},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, Transpose[{rhs4R}]},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], Transpose[{rhs4I}]}
}];
Map[Last,RowReduce[sys]]-Join[\[Lambda]R,\[Lambda]I,\[Mu]R,\[Mu]I,\[ScriptL]R,\[ScriptL]I,\[ScriptM]R,\[ScriptM]I,\[Alpha]R,\[Alpha]I,\[Beta]R,\[Beta]I]//Chop//Norm

sys//MatrixPlot;

sys=SparseArray@ArrayFlatten[{
 {\[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, 0, 0, \[CapitalGamma].inverseMassU.S.W, 0, 0, \[Omega] \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], \[CapitalGamma].inverseMassU.S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, Transpose[{\[CapitalGamma].inverseMassU.rhs1R}]},
 {0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, 0, 0, \[CapitalGamma].inverseMassU.S.W,  -\[Omega] \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, \[CapitalGamma].inverseMassU.S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], Transpose[{\[CapitalGamma].inverseMassU.rhs1I}]},
 {0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].V, 0, \[CapitalXi].inverseMassQ.massQ.W, 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, \[CapitalXi].inverseMassQ.massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], \[Omega] \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], Transpose[{\[CapitalXi].inverseMassQ.rhs2R}]},
 {0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].V, 0, \[CapitalXi].inverseMassQ.massQ.W, 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], -\[Omega] \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], \[CapitalXi].inverseMassQ.massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], Transpose[{\[CapitalXi].inverseMassQ.rhs2I}]},
 {0, 0, 0, 0,  ConjugateTranspose[V].massU.V, 0, 0, -(1/\[Omega])ConjugateTranspose[V].S.W, 0, 0, 0, -1/\[Omega] ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], -1/\[Omega] Transpose[{ConjugateTranspose[V].rhs1I}]},
 {0, 0, 0, 0, 0,  ConjugateTranspose[V].massU.V, 1/\[Omega] ConjugateTranspose[V].S.W, 0, 0, 0, 1/\[Omega] ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 1/\[Omega] Transpose[{ConjugateTranspose[V].rhs1R}]},
 {0, 0, 0, 0, 0, 0, ConjugateTranspose[W].massQ.W, -(1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega] ConjugateTranspose[W].massQ.W), -ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], -1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], -1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].Transpose[{ConjugateTranspose[V].rhs1I}]+Transpose[{ConjugateTranspose[W].rhs2R}]},
 {0, 0, 0, 0, 0, 0, 0, -\[Omega]^2(1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega] ConjugateTranspose[W].massQ.W).Inverse[ConjugateTranspose[W].massQ.W].(1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega] ConjugateTranspose[W].massQ.W)-\[Omega]^2 ConjugateTranspose[W].massQ.W, -\[Omega]^2(1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega] ConjugateTranspose[W].massQ.W).Inverse[ConjugateTranspose[W].massQ.W].ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], \[Omega]^2 ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], \[Omega]^2 (1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega] ConjugateTranspose[W].massQ.W).Inverse[ConjugateTranspose[W].massQ.W].ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]]-\[Omega]^2 1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], \[Omega]^2 (-1)/\[Omega] (1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega] ConjugateTranspose[W].massQ.W).Inverse[ConjugateTranspose[W].massQ.W].ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]]-\[Omega]^2  ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], \[Omega]^2 (1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega] ConjugateTranspose[W].massQ.W).Inverse[ConjugateTranspose[W].massQ.W].(-1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].Transpose[{ConjugateTranspose[V].rhs1I}]+Transpose[{ConjugateTranspose[W].rhs2R}])-\[Omega]^2 (1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].Transpose[{ConjugateTranspose[V].rhs1R}]+Transpose[{ConjugateTranspose[W].rhs2I}])},
 {0, 0, 0, 0, 0, 0, 0, 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, Transpose[{rhs3R}]},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, Transpose[{rhs3I}]},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, Transpose[{rhs4R}]},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], Transpose[{rhs4I}]}
}];
Map[Last,RowReduce[sys]]-Join[\[Lambda]R,\[Lambda]I,\[Mu]R,\[Mu]I,\[ScriptL]R,\[ScriptL]I,\[ScriptM]R,\[ScriptM]I,\[Alpha]R,\[Alpha]I,\[Beta]R,\[Beta]I]//Chop//Norm

sys//MatrixPlot;



sys=SparseArray@ArrayFlatten[{
 {\[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, 0, 0, \[CapitalGamma].inverseMassU.S.W, 0, 0, \[Omega] \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], \[CapitalGamma].inverseMassU.S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, Transpose[{\[CapitalGamma].inverseMassU.rhs1R}]},
 {0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, 0, 0, \[CapitalGamma].inverseMassU.S.W,  -\[Omega] \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, \[CapitalGamma].inverseMassU.S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], Transpose[{\[CapitalGamma].inverseMassU.rhs1I}]},
 {0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].V, 0, \[CapitalXi].inverseMassQ.massQ.W, 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, \[CapitalXi].inverseMassQ.massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], \[Omega] \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], Transpose[{\[CapitalXi].inverseMassQ.rhs2R}]},
 {0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].V, 0, \[CapitalXi].inverseMassQ.massQ.W, 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], -\[Omega] \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], \[CapitalXi].inverseMassQ.massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], Transpose[{\[CapitalXi].inverseMassQ.rhs2I}]},
 {0, 0, 0, 0,  ConjugateTranspose[V].massU.V, 0, 0, -(1/\[Omega])ConjugateTranspose[V].S.W, 0, 0, 0, -1/\[Omega] ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], -1/\[Omega] Transpose[{ConjugateTranspose[V].rhs1I}]},
 {0, 0, 0, 0, 0,  ConjugateTranspose[V].massU.V, 1/\[Omega] ConjugateTranspose[V].S.W, 0, 0, 0, 1/\[Omega] ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 1/\[Omega] Transpose[{ConjugateTranspose[V].rhs1R}]},
 {0, 0, 0, 0, 0, 0, ConjugateTranspose[W].massQ.W, -(1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega] ConjugateTranspose[W].massQ.W), -ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], -1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], -1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].Transpose[{ConjugateTranspose[V].rhs1I}]+Transpose[{ConjugateTranspose[W].rhs2R}]},
 {0, 0, 0, 0, 0, 0, 0, -(ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega]^2 ConjugateTranspose[W].massQ.W).Inverse[ConjugateTranspose[W].massQ.W].(ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega]^2 ConjugateTranspose[W].massQ.W)-\[Omega]^2 ConjugateTranspose[W].massQ.W, -\[Omega](ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega]^2 ConjugateTranspose[W].massQ.W).Inverse[ConjugateTranspose[W].massQ.W].ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], \[Omega]^2 ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], \[Omega] (ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega]^2 ConjugateTranspose[W].massQ.W).Inverse[ConjugateTranspose[W].massQ.W].ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]]-\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], -(ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega]^2 ConjugateTranspose[W].massQ.W).Inverse[ConjugateTranspose[W].massQ.W].ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]]-\[Omega]^2  ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], (ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega]^2 ConjugateTranspose[W].massQ.W).Inverse[ConjugateTranspose[W].massQ.W].(-ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].Transpose[{ConjugateTranspose[V].rhs1I}]+\[Omega] Transpose[{ConjugateTranspose[W].rhs2R}])-(\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].Transpose[{ConjugateTranspose[V].rhs1R}]+\[Omega]^2 Transpose[{ConjugateTranspose[W].rhs2I}])},
 {0, 0, 0, 0, 0, 0, 0, 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, Transpose[{rhs3R}]},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, Transpose[{rhs3I}]},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, Transpose[{rhs4R}]},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], Transpose[{rhs4I}]}
}];
Map[Last,RowReduce[sys]]-Join[\[Lambda]R,\[Lambda]I,\[Mu]R,\[Mu]I,\[ScriptL]R,\[ScriptL]I,\[ScriptM]R,\[ScriptM]I,\[Alpha]R,\[Alpha]I,\[Beta]R,\[Beta]I]//Chop//Norm

sys//MatrixPlot;




sys=SparseArray@ArrayFlatten[{
 {\[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, 0, 0, \[CapitalGamma].inverseMassU.S.W, 0, 0, \[Omega] \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], \[CapitalGamma].inverseMassU.S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, Transpose[{\[CapitalGamma].inverseMassU.rhs1R}]},
 {0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, 0, 0, \[CapitalGamma].inverseMassU.S.W,  -\[Omega] \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, \[CapitalGamma].inverseMassU.S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], Transpose[{\[CapitalGamma].inverseMassU.rhs1I}]},
 {0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].V, 0, \[CapitalXi].inverseMassQ.massQ.W, 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, \[CapitalXi].inverseMassQ.massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], \[Omega] \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], Transpose[{\[CapitalXi].inverseMassQ.rhs2R}]},
 {0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].V, 0, \[CapitalXi].inverseMassQ.massQ.W, 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], -\[Omega] \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], \[CapitalXi].inverseMassQ.massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], Transpose[{\[CapitalXi].inverseMassQ.rhs2I}]},
 {0, 0, 0, 0,  ConjugateTranspose[V].massU.V, 0, 0, -(1/\[Omega])ConjugateTranspose[V].S.W, 0, 0, 0, -1/\[Omega] ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], -1/\[Omega] Transpose[{ConjugateTranspose[V].rhs1I}]},
 {0, 0, 0, 0, 0,  ConjugateTranspose[V].massU.V, 1/\[Omega] ConjugateTranspose[V].S.W, 0, 0, 0, 1/\[Omega] ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 1/\[Omega] Transpose[{ConjugateTranspose[V].rhs1R}]},
 {0, 0, 0, 0, 0, 0, ConjugateTranspose[W].massQ.W, -(1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega] ConjugateTranspose[W].massQ.W), -ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], -1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], -1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].Transpose[{ConjugateTranspose[V].rhs1I}]+Transpose[{ConjugateTranspose[W].rhs2R}]},
 {0, 0, 0, 0, 0, 0, 0, (ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega]^2 ConjugateTranspose[W].massQ.W).Inverse[ConjugateTranspose[W].massQ.W].(ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega]^2 ConjugateTranspose[W].massQ.W)+\[Omega]^2 ConjugateTranspose[W].massQ.W, \[Omega] (ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega]^2 ConjugateTranspose[W].massQ.W).Inverse[ConjugateTranspose[W].massQ.W].ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], -\[Omega]^2ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], -\[Omega](ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega]^2 ConjugateTranspose[W].massQ.W).Inverse[ConjugateTranspose[W].massQ.W].ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]]+\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], (ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega]^2 ConjugateTranspose[W].massQ.W).Inverse[ConjugateTranspose[W].massQ.W].ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]]+\[Omega]^2  ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], -(ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega]^2 ConjugateTranspose[W].massQ.W).Inverse[ConjugateTranspose[W].massQ.W].(-ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].Transpose[{ConjugateTranspose[V].rhs1I}]+\[Omega] Transpose[{ConjugateTranspose[W].rhs2R}])+(\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].Transpose[{ConjugateTranspose[V].rhs1R}]+\[Omega]^2 Transpose[{ConjugateTranspose[W].rhs2I}])},
 {0, 0, 0, 0, 0, 0, 0, 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, Transpose[{rhs3R}]},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, Transpose[{rhs3I}]},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, Transpose[{rhs4R}]},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], Transpose[{rhs4I}]}
}];
sysRAW={
 {\[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, 0, 0, \[CapitalGamma].inverseMassU.S.W, 0, 0, \[Omega] \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], \[CapitalGamma].inverseMassU.S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, Transpose[{\[CapitalGamma].inverseMassU.rhs1R}]},
 {0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, 0, 0, \[CapitalGamma].inverseMassU.S.W,  -\[Omega] \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, \[CapitalGamma].inverseMassU.S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], Transpose[{\[CapitalGamma].inverseMassU.rhs1I}]},
 {0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].V, 0, \[CapitalXi].inverseMassQ.massQ.W, 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, \[CapitalXi].inverseMassQ.massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], \[Omega] \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], Transpose[{\[CapitalXi].inverseMassQ.rhs2R}]},
 {0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].V, 0, \[CapitalXi].inverseMassQ.massQ.W, 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], -\[Omega] \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], \[CapitalXi].inverseMassQ.massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], Transpose[{\[CapitalXi].inverseMassQ.rhs2I}]},
 {0, 0, 0, 0,  ConjugateTranspose[V].massU.V, 0, 0, -(1/\[Omega])ConjugateTranspose[V].S.W, 0, 0, 0, -1/\[Omega] ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], -1/\[Omega] Transpose[{ConjugateTranspose[V].rhs1I}]},
 {0, 0, 0, 0, 0,  ConjugateTranspose[V].massU.V, 1/\[Omega] ConjugateTranspose[V].S.W, 0, 0, 0, 1/\[Omega] ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 1/\[Omega] Transpose[{ConjugateTranspose[V].rhs1R}]},
 {0, 0, 0, 0, 0, 0, ConjugateTranspose[W].massQ.W, -(1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega] ConjugateTranspose[W].massQ.W), -ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], -1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], -1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].Transpose[{ConjugateTranspose[V].rhs1I}]+Transpose[{ConjugateTranspose[W].rhs2R}]},
 {0, 0, 0, 0, 0, 0, 0, (ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega]^2 ConjugateTranspose[W].massQ.W).Inverse[ConjugateTranspose[W].massQ.W].(ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega]^2 ConjugateTranspose[W].massQ.W)+\[Omega]^2 ConjugateTranspose[W].massQ.W, \[Omega] (ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega]^2 ConjugateTranspose[W].massQ.W).Inverse[ConjugateTranspose[W].massQ.W].ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], -\[Omega]^2ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], -\[Omega](ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega]^2 ConjugateTranspose[W].massQ.W).Inverse[ConjugateTranspose[W].massQ.W].ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]]+\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], (ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega]^2 ConjugateTranspose[W].massQ.W).Inverse[ConjugateTranspose[W].massQ.W].ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]]+\[Omega]^2  ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], -(ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega]^2 ConjugateTranspose[W].massQ.W).Inverse[ConjugateTranspose[W].massQ.W].(-ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].Transpose[{ConjugateTranspose[V].rhs1I}]+\[Omega] Transpose[{ConjugateTranspose[W].rhs2R}])+(\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].Transpose[{ConjugateTranspose[V].rhs1R}]+\[Omega]^2 Transpose[{ConjugateTranspose[W].rhs2I}])},
 {0, 0, 0, 0, 0, 0, 0, 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, Transpose[{rhs3R}]},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, Transpose[{rhs3I}]},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, Transpose[{rhs4R}]},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], Transpose[{rhs4I}]}
};

Map[Last,RowReduce[sys]]-Join[\[Lambda]R,\[Lambda]I,\[Mu]R,\[Mu]I,\[ScriptL]R,\[ScriptL]I,\[ScriptM]R,\[ScriptM]I,\[Alpha]R,\[Alpha]I,\[Beta]R,\[Beta]I]//Chop//Norm

sys//MatrixPlot;







ConjugateTranspose[W].massQ.W//Chop//MatrixPlot


CholeskyDecomposition[ConjugateTranspose[W].massQ.W]//Chop//MatrixPlot


?LinearSolve


(*Operators Solvers and Preconditioners*)
Vsmuv[x_] := Dot[ConjugateTranspose[V],Dot[massU,Dot[V.x]]]
VsmuVPre[x_] := Dot[ConjugateTranspose[V],Dot[inverseMassU,Dot[V.x]]]
VsmuvInv[x_] := LinearSolve[Vsmuv,x,Method->{"Krylov",Method->"ConjugateGradient",Preconditioner-> VsmuVPre,"Tolerance"-> 10^-16}]
VsmuvInv /: Dot[VsmuvInv,x_List] := VsmuvInv[x]
Op1[x_] := ConjugateTranspose[W].ConjugateTranspose[S].V.VsmuvInv.(ConjugateTranspose[V].S.W.x)-(\[Omega]^2 ConjugateTranspose[W].massQ.W).x
Op1/: Dot[Op1,x_List] := Op1[x];

Wsmqw[x_] := Dot[ConjugateTranspose[W],Dot[massQ,Dot[W,x]]]
WsmqwPre[x_] := Dot[ConjugateTranspose[W],Dot[inverseMassQ,Dot[W,x]]]
WsmqwInv[x_] := LinearSolve[Wsmqw,x,Method->{"Krylov",Method->"ConjugateGradient",Preconditioner->WsmqwPre}]
WsmqwInv /: Dot[WsmqwInv,x_List] := WsmqwInv[x]

XiInvmqXis[x_] := \[CapitalXi].(inverseMassQ.(ConjugateTranspose[\[CapitalXi]].x))
XiInvmqXisPre[x_] := \[CapitalXi].(massQ.(ConjugateTranspose[\[CapitalXi]].x))
XiInvmqXisInv[x_] := LinearSolve[XiInvmqXis,x,Method->{"Krylov",Method->"ConjugateGradient",Preconditioner->XiInvmqXisPre}]
XiInvmqXisInv /: Dot[XiInvmqXisInv,x_List] := XiInvmqXisInv[x]

GaInvmuGas[x_] := \[CapitalGamma].(inverseMassU.(ConjugateTranspose[\[CapitalGamma]].x));
GaInvmuGasPre[x_] := \[CapitalGamma].(massU.(ConjugateTranspose[\[CapitalGamma]].x));
GaInvmuGasInv[x_] := LinearSolve[GaInvmuGas,x,Method->{"Krylov",Method->"ConjugateGradient",Preconditioner->GaInvmuGasPre}]
GaInvmuGasInv /: Dot[GaInvmuGasInv,x_List] := GaInvmuGasInv[x]

big[x_] := \[Omega]^2 (ConjugateTranspose[W].massQ.W).x+Op1.(WsmqwInv.(Op1.x))
bigPre[x_] := x;
bigInv[x_] := LinearSolve[big,x,Method->{"Krylov",Method->"ConjugateGradient",Tolerance->10^-12}]
bigInv /: Dot[bigInv,x_List] := bigInv[x];


(ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega]^2 ConjugateTranspose[W].massQ.W).Inverse[ConjugateTranspose[W].massQ.W].(ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega]^2 ConjugateTranspose[W].massQ.W)+\[Omega]^2 ConjugateTranspose[W].massQ.W





?LinearSolve


(*Construct New RHS's*)
On[Assert]
RHS1=\[CapitalGamma].inverseMassU.rhs1R;
RHS2=\[CapitalGamma].inverseMassU.rhs1I;
RHS3=\[CapitalXi].inverseMassQ.rhs2R;
RHS4=\[CapitalXi].inverseMassQ.rhs2I;
RHS5=-1/\[Omega] ConjugateTranspose[V].rhs1I;
RHS6=1/\[Omega] ConjugateTranspose[V].rhs1R;
RHS7=-1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].rhs1I+ConjugateTranspose[W].rhs2R;
RHS7Test=-1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.VsmuvInv[ConjugateTranspose[V].rhs1I]+ConjugateTranspose[W].rhs2R;
Assert[Norm[RHS7Test-RHS7]<10^-10]
RHS8=
(
-(ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega]^2 ConjugateTranspose[W].massQ.W).Inverse[ConjugateTranspose[W].massQ.W].(-ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].rhs1I+\[Omega] ConjugateTranspose[W].rhs2R)
+(\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].rhs1R+\[Omega]^2 ConjugateTranspose[W].rhs2I)
);
RHS8Test=
(-(ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega]^2 ConjugateTranspose[W].massQ.W).Inverse[ConjugateTranspose[W].massQ.W].(-ConjugateTranspose[W].ConjugateTranspose[S].V.VsmuvInv.(ConjugateTranspose[V].rhs1I)+\[Omega] ConjugateTranspose[W].rhs2R)
+(\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.VsmuvInv.(ConjugateTranspose[V].rhs1R)+\[Omega]^2 ConjugateTranspose[W].rhs2I));


RHS8Test=
(-Op1.(Inverse[ConjugateTranspose[W].massQ.W].(-ConjugateTranspose[W].ConjugateTranspose[S].V.VsmuvInv.(ConjugateTranspose[V].rhs1I)+\[Omega] ConjugateTranspose[W].rhs2R))
+(\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.VsmuvInv.(ConjugateTranspose[V].rhs1R)+\[Omega]^2 ConjugateTranspose[W].rhs2I));

RHS8Test=
(-Op1.(WsmqwInv.(-ConjugateTranspose[W].ConjugateTranspose[S].V.VsmuvInv.(ConjugateTranspose[V].rhs1I)+\[Omega] ConjugateTranspose[W].rhs2R))
+(\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.VsmuvInv.(ConjugateTranspose[V].rhs1R)+\[Omega]^2 ConjugateTranspose[W].rhs2I));

Assert[Norm[RHS8Test-RHS8]<10^-10]
RHS9=rhs3R;
RHS10=rhs3I;
RHS11=rhs4R;
RHS12=rhs4I; 


(*Find \[Beta]I *)
LinearSolve[\[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]],RHS12]-\[Beta]I//Chop//Norm
XiInvmqXisInv.RHS12-\[Beta]I//Chop//Norm
(*Change RHS *)
RHS1=RHS1;
RHS2=RHS2-\[CapitalGamma].inverseMassU.S.inverseMassQ.ConjugateTranspose[\[CapitalXi]].\[Beta]I;

RHS3=RHS3-\[Omega] \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]].\[Beta]I;
RHS4-=\[CapitalXi].inverseMassQ.massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]].\[Beta]I;
RHS5-=-1/\[Omega] ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]].\[Beta]I;

RHS6=RHS6;


(*RHS7-=0(-1/\[Omega]W^\[ConjugateTranspose].S^\[ConjugateTranspose].V.Inverse[V^\[ConjugateTranspose].massU.V].V^\[ConjugateTranspose].S.inverseMassQ.\[CapitalXi]^\[ConjugateTranspose]).\[Beta]I;*)
RHS7-=(-1/\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.VsmuvInv.(ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]]).\[Beta]I);

(*RHS8-=0((W^\[ConjugateTranspose].S^\[ConjugateTranspose].V.Inverse[V^\[ConjugateTranspose].massU.V].V^\[ConjugateTranspose].S.W-\[Omega]^2 W^\[ConjugateTranspose].massQ.W).Inverse[W^\[ConjugateTranspose].massQ.W].W^\[ConjugateTranspose].S^\[ConjugateTranspose].V.Inverse[V^\[ConjugateTranspose].massU.V].V^\[ConjugateTranspose].S.inverseMassQ.\[CapitalXi]^\[ConjugateTranspose]+\[Omega]^2 W^\[ConjugateTranspose].massQ.inverseMassQ.\[CapitalXi]^\[ConjugateTranspose]).\[Beta]I;*)

RHS8-=(Op1.(WsmqwInv.(ConjugateTranspose[W].ConjugateTranspose[S].V.(VsmuvInv.(ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]].\[Beta]I)))))+\[Omega]^2  ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]].\[Beta]I;

RHS9=RHS9;
RHS10=RHS10;
RHS11=RHS11;

(*Find \[Beta]R *)
LinearSolve[\[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]],RHS11]-\[Beta]R//Chop//Norm
XiInvmqXisInv.RHS11-\[Beta]R//Chop//Norm

(*Change RHS *)
RHS1-=(\[CapitalGamma].inverseMassU.S.inverseMassQ.ConjugateTranspose[\[CapitalXi]]).\[Beta]R;
RHS2=RHS2;
RHS3-=(\[CapitalXi].inverseMassQ.massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]]).\[Beta]R;
RHS4-=(-\[Omega] \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]]).\[Beta]R;
RHS5=RHS5;
RHS6-=1/\[Omega] ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]].\[Beta]R;
RHS7-=-ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]].\[Beta]R;


(*new=-\[Omega](Op1.(WsmqwInv.(W^\[ConjugateTranspose].massQ.inverseMassQ.\[CapitalXi]^\[ConjugateTranspose].\[Beta]R)))+\[Omega] W^\[ConjugateTranspose].S^\[ConjugateTranspose].V.VsmuvInv.(V^\[ConjugateTranspose].S.inverseMassQ.\[CapitalXi]^\[ConjugateTranspose].\[Beta]R);*)

(*RHS8-=(-\[Omega](W^\[ConjugateTranspose].S^\[ConjugateTranspose].V.Inverse[V^\[ConjugateTranspose].massU.V].V^\[ConjugateTranspose].S.W-\[Omega]^2 W^\[ConjugateTranspose].massQ.W).Inverse[W^\[ConjugateTranspose].massQ.W].W^\[ConjugateTranspose].massQ.inverseMassQ.\[CapitalXi]^\[ConjugateTranspose]+\[Omega] W^\[ConjugateTranspose].S^\[ConjugateTranspose].V.Inverse[V^\[ConjugateTranspose].massU.V].V^\[ConjugateTranspose].S.inverseMassQ.\[CapitalXi]^\[ConjugateTranspose]).\[Beta]R;*)

RHS8-=(-\[Omega](Op1.(WsmqwInv.(ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]].\[Beta]R)))+\[Omega] ConjugateTranspose[W].ConjugateTranspose[S].V.VsmuvInv.(ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]].\[Beta]R));
RHS9=RHS9;
RHS10=RHS10;

(*Find \[Alpha]I *)
LinearSolve[\[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]],RHS10]-\[Alpha]I//Chop//Norm
GaInvmuGasInv.RHS10-\[Alpha]I//Chop//Norm
(*Change RHS *)

RHS1-=(\[Omega] \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]]).\[Alpha]I;
RHS2=RHS2;
RHS3=RHS3;
RHS4-=(-\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]]).\[Alpha]I;
RHS5=RHS5;
RHS6=RHS6;
RHS7=RHS7;
RHS8-=(-\[Omega]^2ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]]).\[Alpha]I;
RHS9=RHS9;

(*Fine \[Alpha]R *)

LinearSolve[\[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]],RHS9]-\[Alpha]R//Chop//Norm
GaInvmuGasInv.RHS9-\[Alpha]R//Chop//Norm
(*Change RHS *)

RHS1=RHS1;
RHS2-=(-\[Omega] \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]]).\[Alpha]R;
RHS3-=(-\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]]).\[Alpha]R;
RHS4=RHS4;
RHS5=RHS5;
RHS6=RHS6;
RHS7-=(-ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]]).\[Alpha]R;
(*RHS8-=(\[Omega](W^\[ConjugateTranspose].S^\[ConjugateTranspose].V.Inverse[V^\[ConjugateTranspose].massU.V].V^\[ConjugateTranspose].S.W-\[Omega]^2 W^\[ConjugateTranspose].massQ.W).Inverse[W^\[ConjugateTranspose].massQ.W].W^\[ConjugateTranspose].S^\[ConjugateTranspose].inverseMassU.\[CapitalGamma]^\[ConjugateTranspose]).\[Alpha]R;*)
RHS8-= \[Omega] Op1.WsmqwInv.(ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]].\[Alpha]R);

(*Find \[ScriptM]I *)
\[ScriptCapitalC]=SparseArray[CholeskyDecomposition[ConjugateTranspose[W].massQ.W]];



LinearSolve[(ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega]^2 ConjugateTranspose[W].massQ.W).Inverse[ConjugateTranspose[W].massQ.W].(ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W-\[Omega]^2 ConjugateTranspose[W].massQ.W)+\[Omega]^2 ConjugateTranspose[W].massQ.W,RHS8]-\[ScriptM]I//Chop//Norm
ans=(bigInv.RHS8);
ans2=LinearSolve[big,RHS8,Method->{"Krylov",Method->"ConjugateGradient",Preconditioner->bigPre,StartingVector->ans }];
ans3=LinearSolve[big,RHS8,Method->{"Krylov",Method->"ConjugateGradient",Preconditioner->bigPre,StartingVector->ans2 }];

(*Change RHS *)
RHS1=RHS1;
RHS2-=(\[CapitalGamma].inverseMassU.S.W).\[ScriptM]I;
RHS3=RHS3;
RHS4-=(\[CapitalXi].inverseMassQ.massQ.W).\[ScriptM]I;
RHS5-=(-(1/\[Omega])ConjugateTranspose[V].S.W).\[ScriptM]I;
RHS6=RHS6;

(*RHS7-=(-(1/\[Omega]W^\[ConjugateTranspose].S^\[ConjugateTranspose].V.Inverse[V^\[ConjugateTranspose].massU.V].V^\[ConjugateTranspose].S.W-\[Omega] W^\[ConjugateTranspose].massQ.W)).\[ScriptM]I;*)

RHS7-=-1/\[Omega] Op1.\[ScriptM]I;


\[ScriptCapitalR]=4;


{\[CapitalLambda],P}=Eigensystem[{ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W,ConjugateTranspose[W].massQ.W}];
\[CapitalLambda]=DiagonalMatrix[\[CapitalLambda]];
P=Transpose[P];


(ConjugateTranspose[W].massQ.W).P.\[CapitalLambda].Inverse[P]-ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W//Chop


Transpose[P].(ConjugateTranspose[W].massQ.W).P//Chop
Transpose[P].ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W.P//Chop


/


(*Find \[ScriptM]R *)
LinearSolve[ConjugateTranspose[W].massQ.W,RHS7]-\[ScriptM]R//Chop//Norm

(*Change RHS *)

RHS1-= (\[CapitalGamma].inverseMassU.S.W).\[ScriptM]R;
RHS2= RHS2;
RHS3-=(\[CapitalXi].inverseMassQ.massQ.W).\[ScriptM]R;
RHS4=RHS4;
RHS5=RHS5;
RHS6-=(1/\[Omega] ConjugateTranspose[V].S.W).\[ScriptM]R;

(* Find \[ScriptL]I *)

LinearSolve[ConjugateTranspose[V].massU.V,RHS6]-\[ScriptL]I//Chop//Norm

(*Change RHS *)

RHS1=RHS1;
RHS2=RHS2;
RHS3=RHS3;
RHS4-=(-\[CapitalXi].inverseMassQ.ConjugateTranspose[S].V).\[ScriptL]I;
RHS5=RHS5;

(*Find \[ScriptL]R *)

LinearSolve[ConjugateTranspose[V].massU.V,RHS5]-\[ScriptL]R//Chop//Norm
(*Change RHS *)
RHS1=RHS1;
RHS2=RHS2;
RHS3-=(-\[CapitalXi].inverseMassQ.ConjugateTranspose[S].V).\[ScriptL]R;
RHS4=RHS4;

(*Find \[Mu]I *)

LinearSolve[\[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]],RHS4]-\[Mu]I//Chop//Norm
(*Change RHS *)
RHS1=RHS1;
RHS2=RHS2;
RHS3=RHS3;

(*Find \[Mu]R *)

LinearSolve[\[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]],RHS3]-\[Mu]R//Chop//Norm

(*Change RHS *)

RHS1=RHS1;
RHS2=RHS2;

(*Fine \[Lambda]I *)

LinearSolve[\[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]],RHS2]-\[Lambda]I//Chop//Norm
(*Change RHS *)
RHS1=RHS1;
LinearSolve[\[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]],RHS1]-\[Lambda]R//Chop//Norm


\[DoubleStruckCapitalA]=(ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W);
\[DoubleStruckCapitalB]=ConjugateTranspose[W].massQ.W;
inv\[DoubleStruckCapitalB]=Inverse[ConjugateTranspose[W].massQ.W];
\[DoubleStruckCapitalC]=SparseArray[CholeskyDecomposition[\[DoubleStruckCapitalB]]];
inv\[DoubleStruckCapitalC]=SparseArray@Inverse[\[DoubleStruckCapitalC]];
\[ScriptCapitalI]=SparseArray[IdentityMatrix[Length[\[DoubleStruckCapitalA]]]];
{\[CapitalLambda],Q}=Eigensystem[Transpose[inv\[DoubleStruckCapitalC]].\[DoubleStruckCapitalA].inv\[DoubleStruckCapitalC]];
\[CapitalLambda]=SparseArray[DiagonalMatrix[\[CapitalLambda]]];





Inverse[\[ScriptCapitalI]-UU.\[CapitalSigma]\[CapitalSigma].ConjugateTranspose[VV]/\[Sigma]].Inverse[\[ScriptCapitalI]-UU.\[CapitalSigma]\[CapitalSigma].ConjugateTranspose[VV]/Conjugate[\[Sigma]]]-Conjugate[y].y//Chop
Conjugate[y].y-x//Norm
Inverse[\[ScriptCapitalI]-UU.\[CapitalSigma]\[CapitalSigma].ConjugateTranspose[VV]/\[Sigma]].Inverse[\[ScriptCapitalI]-UU.\[CapitalSigma]\[CapitalSigma].ConjugateTranspose[VV]/Conjugate[\[Sigma]]]-x//Norm


Inverse[\[ScriptCapitalI]-UU.\[CapitalSigma]\[CapitalSigma].ConjugateTranspose[VV]/\[Sigma]].Inverse[\[ScriptCapitalI]-UU.\[CapitalSigma]\[CapitalSigma].ConjugateTranspose[VV]/Conjugate[\[Sigma]]]-UU.\[ScriptCapitalT].Transpose[VV]//Norm


{\[CapitalLambda],P}=Eigensystem[{ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W,ConjugateTranspose[W].massQ.W},1];
{{Ua,Ub},{\[CapitalSigma]a,\[CapitalSigma]b},v}=SingularValueDecomposition[{ConjugateTranspose[W].ConjugateTranspose[S].V.Inverse[ConjugateTranspose[V].massU.V].ConjugateTranspose[V].S.W,ConjugateTranspose[W].massQ.W}];
\[CapitalLambda]=DiagonalMatrix[\[CapitalLambda]];
P=Transpose[P];


ConjugateTranspose[W].massQ.W//Length


ConjugateTranspose[Ua].Ua//Chop//MatrixPlot
ConjugateTranspose[Ua].v//Chop//MatrixPlot
ConjugateTranspose[Ub].ConjugateTranspose[W].massQ.W.v//Chop//MatrixPlot


ConjugateTranspose[Ub].ConjugateTranspose[W].massQ.W.v-ConjugateTranspose[v].v//Chop//MatrixPlot


ConjugateTranspose[v].v//Chop//MatrixPlot


ConjugateTranspose[W].massQ.W.P.\[CapitalLambda].ConjugateTranspose[P]//MatrixPlot
(*W^\[ConjugateTranspose].massQ.W.P.\[CapitalLambda].Inverse[P]-(W^\[ConjugateTranspose].S^\[ConjugateTranspose].V.Inverse[V^\[ConjugateTranspose].massU.V].V^\[ConjugateTranspose].S.W)*)


Transpose[P].ConjugateTranspose[W].massQ.W.P.\[CapitalLambda]//Chop


2.23066^2


\[Omega]=1/10;
\[Sigma]:=\[Omega]^2+I \[Omega];
\[ScriptCapitalR]=5;

\[DoubleStruckCapitalI]=IdentityMatrix[\[ScriptCapitalR]];
z=(\[DoubleStruckCapitalA]-\[Omega]^2 \[DoubleStruckCapitalB]).inv\[DoubleStruckCapitalB].(\[DoubleStruckCapitalA]-\[Omega]^2 \[DoubleStruckCapitalB])+\[Omega]^2 \[DoubleStruckCapitalB];

zz=Transpose[\[DoubleStruckCapitalC]].((Transpose[inv\[DoubleStruckCapitalC]].\[DoubleStruckCapitalA].inv\[DoubleStruckCapitalC]-\[Omega]^2 \[ScriptCapitalI]).(Transpose[inv\[DoubleStruckCapitalC]].\[DoubleStruckCapitalA].inv\[DoubleStruckCapitalC]-\[Omega]^2 \[ScriptCapitalI])+\[Omega]^2 \[ScriptCapitalI]).\[DoubleStruckCapitalC];
\[ScriptCapitalA]=Transpose[inv\[DoubleStruckCapitalC]].\[DoubleStruckCapitalA].inv\[DoubleStruckCapitalC];
zz=Transpose[\[DoubleStruckCapitalC]].((\[ScriptCapitalA]-\[Omega]^2 \[ScriptCapitalI]).(\[ScriptCapitalA]-\[Omega]^2 \[ScriptCapitalI])+\[Omega]^2 \[ScriptCapitalI]).\[DoubleStruckCapitalC];

zz=Transpose[\[DoubleStruckCapitalC]].(\[ScriptCapitalA]-\[Sigma] \[ScriptCapitalI]).(\[ScriptCapitalA]-Conjugate[\[Sigma]]\[ScriptCapitalI]).\[DoubleStruckCapitalC];

zz=\[Sigma] Conjugate[\[Sigma]] Transpose[\[DoubleStruckCapitalC]].(\[ScriptCapitalA]/\[Sigma]-\[ScriptCapitalI]).(\[ScriptCapitalA]/Conjugate[\[Sigma]]-\[ScriptCapitalI]).\[DoubleStruckCapitalC];
zz=\[Sigma] Conjugate[\[Sigma]] Transpose[\[DoubleStruckCapitalC]].(\[ScriptCapitalI]-\[ScriptCapitalA]/\[Sigma]).(\[ScriptCapitalI]-\[ScriptCapitalA]/Conjugate[\[Sigma]]).\[DoubleStruckCapitalC];
{UU,\[CapitalSigma]\[CapitalSigma],VV}=SingularValueDecomposition[\[ScriptCapitalA],\[ScriptCapitalR]];
(\[ScriptCapitalI]-\[ScriptCapitalA]/\[Sigma]).(\[ScriptCapitalI]-\[ScriptCapitalA]/Conjugate[\[Sigma]])//Norm

Inverse[\[ScriptCapitalI]-UU.\[CapitalSigma]\[CapitalSigma].ConjugateTranspose[VV]/\[Sigma]].(\[ScriptCapitalI]-\[ScriptCapitalA]/\[Sigma]).(\[ScriptCapitalI]-\[ScriptCapitalA]/Conjugate[\[Sigma]]).Inverse[\[ScriptCapitalI]-UU.\[CapitalSigma]\[CapitalSigma].ConjugateTranspose[VV]/Conjugate[\[Sigma]]]//Norm

Inverse[\[ScriptCapitalI]-UU.\[CapitalSigma]\[CapitalSigma].ConjugateTranspose[VV]/\[Sigma]].Inverse[\[ScriptCapitalI]-UU.\[CapitalSigma]\[CapitalSigma].ConjugateTranspose[VV]/Conjugate[\[Sigma]]].(\[ScriptCapitalI]-\[ScriptCapitalA]/\[Sigma]).(\[ScriptCapitalI]-\[ScriptCapitalA]/Conjugate[\[Sigma]])//Norm

z-zz//Norm

A=\[ScriptCapitalI];
U=-UU/Conjugate[\[Sigma]];
CC=\[CapitalSigma]\[CapitalSigma];
V=ConjugateTranspose[VV];

y=Inverse[\[ScriptCapitalI]-UU.\[CapitalSigma]\[CapitalSigma].ConjugateTranspose[VV]/Conjugate[\[Sigma]]];
yy=Inverse[A+U.CC.V];

yy=Inverse[A]-Inverse[A].U.Inverse[Inverse[CC]+V.Inverse[A].U].V.Inverse[A];
yy=\[ScriptCapitalI]-\[ScriptCapitalI].U.Inverse[Inverse[CC]+V.\[ScriptCapitalI].U].V.\[ScriptCapitalI];

yy=\[ScriptCapitalI]-U.Inverse[Inverse[CC]+V.U].V;
yy=\[ScriptCapitalI]+(UU/Conjugate[\[Sigma]]).Inverse[Inverse[CC]-V.(UU/Conjugate[\[Sigma]])].V;

yy=\[ScriptCapitalI]+(UU/Conjugate[\[Sigma]]).Inverse[Inverse[CC]-ConjugateTranspose[VV].(UU/Conjugate[\[Sigma]])].ConjugateTranspose[VV];

yy=\[ScriptCapitalI]+(UU/Conjugate[\[Sigma]]).Inverse[Inverse[CC]-1/Conjugate[\[Sigma]] \[DoubleStruckCapitalI]].ConjugateTranspose[VV];

\[ScriptCapitalH]=DiagonalMatrix[Map[1/#&,Map[1/#-1/Conjugate[\[Sigma]]&,Diagonal[\[CapitalSigma]\[CapitalSigma]]]]]//Chop;

yy=\[ScriptCapitalI]+(UU/Conjugate[\[Sigma]]).\[ScriptCapitalH].ConjugateTranspose[VV];


y-yy//Norm

x=Conjugate[y].y;
xx=(\[ScriptCapitalI]+(UU/\[Sigma]).Inverse[Inverse[CC]-1/\[Sigma] \[DoubleStruckCapitalI]].ConjugateTranspose[VV]).(\[ScriptCapitalI]+(UU/Conjugate[\[Sigma]]).Inverse[Inverse[CC]-1/Conjugate[\[Sigma]] \[DoubleStruckCapitalI]].ConjugateTranspose[VV]);

xx=Plus[
(\[ScriptCapitalI]+0(UU/\[Sigma]).Inverse[Inverse[CC]-1/\[Sigma] \[DoubleStruckCapitalI]].ConjugateTranspose[VV]).(\[ScriptCapitalI]+0(UU/Conjugate[\[Sigma]]).Inverse[Inverse[CC]-1/Conjugate[\[Sigma]] \[DoubleStruckCapitalI]].ConjugateTranspose[VV]),
(\[ScriptCapitalI]+0(UU/\[Sigma]).Inverse[Inverse[CC]-1/\[Sigma] \[DoubleStruckCapitalI]].ConjugateTranspose[VV]).(0\[ScriptCapitalI]+(UU/Conjugate[\[Sigma]]).Inverse[Inverse[CC]-1/Conjugate[\[Sigma]] \[DoubleStruckCapitalI]].ConjugateTranspose[VV]),
(0\[ScriptCapitalI]+(UU/\[Sigma]).Inverse[Inverse[CC]-1/\[Sigma] \[DoubleStruckCapitalI]].ConjugateTranspose[VV]).(\[ScriptCapitalI]+0(UU/Conjugate[\[Sigma]]).Inverse[Inverse[CC]-1/Conjugate[\[Sigma]] \[DoubleStruckCapitalI]].ConjugateTranspose[VV]),
(0\[ScriptCapitalI]+(UU/\[Sigma]).Inverse[Inverse[CC]-1/\[Sigma] \[DoubleStruckCapitalI]].ConjugateTranspose[VV]).(0\[ScriptCapitalI]+(UU/Conjugate[\[Sigma]]).Inverse[Inverse[CC]-1/Conjugate[\[Sigma]] \[DoubleStruckCapitalI]].ConjugateTranspose[VV])];

xx=Plus[
(\[ScriptCapitalI]).(\[ScriptCapitalI]),
(\[ScriptCapitalI]).((UU/Conjugate[\[Sigma]]).Inverse[Inverse[CC]-1/Conjugate[\[Sigma]] \[DoubleStruckCapitalI]].ConjugateTranspose[VV]),
((UU/\[Sigma]).Inverse[Inverse[CC]-1/\[Sigma] \[DoubleStruckCapitalI]].ConjugateTranspose[VV]).(\[ScriptCapitalI]),
((UU/\[Sigma]).Inverse[Inverse[CC]-1/\[Sigma] \[DoubleStruckCapitalI]].ConjugateTranspose[VV]).((UU/Conjugate[\[Sigma]]).Inverse[Inverse[CC]-1/Conjugate[\[Sigma]] \[DoubleStruckCapitalI]].ConjugateTranspose[VV])];

xx=Plus[
\[ScriptCapitalI],
((UU/Conjugate[\[Sigma]]).Inverse[Inverse[CC]-1/Conjugate[\[Sigma]] \[DoubleStruckCapitalI]].ConjugateTranspose[VV]),
((UU/\[Sigma]).Inverse[Inverse[CC]-1/\[Sigma] \[DoubleStruckCapitalI]].ConjugateTranspose[VV]),
((UU/\[Sigma]).Inverse[Inverse[CC]-1/\[Sigma] \[DoubleStruckCapitalI]].ConjugateTranspose[VV]).((UU/Conjugate[\[Sigma]]).Inverse[Inverse[CC]-1/Conjugate[\[Sigma]] \[DoubleStruckCapitalI]].ConjugateTranspose[VV])];

x-xx//Norm

xx=Plus[
\[ScriptCapitalI],
((UU/Conjugate[\[Sigma]]).Inverse[Inverse[CC]-1/Conjugate[\[Sigma]] \[DoubleStruckCapitalI]].ConjugateTranspose[VV]),
((UU/\[Sigma]).Inverse[Inverse[CC]-1/\[Sigma] \[DoubleStruckCapitalI]].ConjugateTranspose[VV]),
1/(\[Sigma] Conjugate[\[Sigma]]) UU.Inverse[Inverse[CC]-1/\[Sigma] \[DoubleStruckCapitalI]].ConjugateTranspose[VV].UU.Inverse[Inverse[CC]-1/Conjugate[\[Sigma]] \[DoubleStruckCapitalI]].ConjugateTranspose[VV]];

xx=Plus[
\[ScriptCapitalI],
((UU/Conjugate[\[Sigma]]).Inverse[Inverse[CC]-1/Conjugate[\[Sigma]] \[DoubleStruckCapitalI]].ConjugateTranspose[VV]),
((UU/\[Sigma]).Inverse[Inverse[CC]-1/\[Sigma] \[DoubleStruckCapitalI]].ConjugateTranspose[VV]),
1/(\[Sigma] Conjugate[\[Sigma]]) UU.Inverse[Inverse[CC]-1/\[Sigma] \[DoubleStruckCapitalI]].ConjugateTranspose[VV].UU.Inverse[Inverse[CC]-1/Conjugate[\[Sigma]] \[DoubleStruckCapitalI]].ConjugateTranspose[VV]];

xx=Plus[
\[ScriptCapitalI],
UU.(1/Conjugate[\[Sigma]] Inverse[Inverse[CC]-1/Conjugate[\[Sigma]] \[DoubleStruckCapitalI]]).ConjugateTranspose[VV],
UU.(1/\[Sigma] Inverse[Inverse[CC]-1/\[Sigma] \[DoubleStruckCapitalI]]).ConjugateTranspose[VV],
1/(\[Sigma] Conjugate[\[Sigma]]) UU.Inverse[Inverse[CC]-1/\[Sigma] \[DoubleStruckCapitalI]].ConjugateTranspose[VV].UU.Inverse[Inverse[CC]-1/Conjugate[\[Sigma]] \[DoubleStruckCapitalI]].ConjugateTranspose[VV]];

xx=Plus[
\[ScriptCapitalI],
UU.(1/Conjugate[\[Sigma]] Inverse[Inverse[CC]-1/Conjugate[\[Sigma]] \[DoubleStruckCapitalI]]+1/\[Sigma] Inverse[Inverse[CC]-1/\[Sigma] \[DoubleStruckCapitalI]]).ConjugateTranspose[VV],
0 UU.(1/\[Sigma] Inverse[Inverse[CC]-1/\[Sigma] \[DoubleStruckCapitalI]]).ConjugateTranspose[VV],
1/(\[Sigma] Conjugate[\[Sigma]]) UU.Inverse[Inverse[CC]-1/\[Sigma] \[DoubleStruckCapitalI]].ConjugateTranspose[VV].UU.Inverse[Inverse[CC]-1/Conjugate[\[Sigma]] \[DoubleStruckCapitalI]].ConjugateTranspose[VV]];

xx=Plus[
\[ScriptCapitalI],
UU.(Inverse[Conjugate[\[Sigma]]Inverse[CC]-\[DoubleStruckCapitalI]]+Inverse[\[Sigma] Inverse[CC]-\[DoubleStruckCapitalI]]).ConjugateTranspose[VV],
1/(\[Sigma] Conjugate[\[Sigma]]) UU.Inverse[Inverse[CC]-1/\[Sigma] \[DoubleStruckCapitalI]].ConjugateTranspose[VV].UU.Inverse[Inverse[CC]-1/Conjugate[\[Sigma]] \[DoubleStruckCapitalI]].ConjugateTranspose[VV]];

xx=Plus[
\[ScriptCapitalI],
UU.(Inverse[Conjugate[\[Sigma]]Inverse[CC]-\[DoubleStruckCapitalI]]+Inverse[\[Sigma] Inverse[CC]-\[DoubleStruckCapitalI]]).ConjugateTranspose[VV],
UU.Inverse[\[Sigma] Inverse[CC]-\[DoubleStruckCapitalI]].ConjugateTranspose[VV].UU.Inverse[Conjugate[\[Sigma]]Inverse[CC]-\[DoubleStruckCapitalI]].ConjugateTranspose[VV]];

xx=Plus[
\[ScriptCapitalI],
UU.(Inverse[Conjugate[\[Sigma]]Inverse[CC]-\[DoubleStruckCapitalI]]+Inverse[\[Sigma] Inverse[CC]-\[DoubleStruckCapitalI]]).ConjugateTranspose[VV],
UU.Inverse[\[Sigma] Inverse[CC]-\[DoubleStruckCapitalI]].Inverse[Conjugate[\[Sigma]]Inverse[CC]-\[DoubleStruckCapitalI]].ConjugateTranspose[VV]];

xx=\[ScriptCapitalI]+
Plus[
UU.(Inverse[Conjugate[\[Sigma]]Inverse[CC]-\[DoubleStruckCapitalI]]+Inverse[\[Sigma] Inverse[CC]-\[DoubleStruckCapitalI]]).ConjugateTranspose[VV],
UU.Inverse[\[Sigma] Inverse[CC]-\[DoubleStruckCapitalI]].Inverse[Conjugate[\[Sigma]]Inverse[CC]-\[DoubleStruckCapitalI]].ConjugateTranspose[VV]];


x-xx//Norm;

xx=\[ScriptCapitalI]+
UU.Plus[
(Inverse[Conjugate[\[Sigma]]Inverse[CC]-\[DoubleStruckCapitalI]]+Inverse[\[Sigma] Inverse[CC]-\[DoubleStruckCapitalI]]).ConjugateTranspose[VV],
Inverse[\[Sigma] Inverse[CC]-\[DoubleStruckCapitalI]].Inverse[Conjugate[\[Sigma]]Inverse[CC]-\[DoubleStruckCapitalI]].ConjugateTranspose[VV]];


x-xx//Norm;

xx=\[ScriptCapitalI]+
UU.Plus[
(Inverse[Conjugate[\[Sigma]]Inverse[CC]-\[DoubleStruckCapitalI]]+Inverse[\[Sigma] Inverse[CC]-\[DoubleStruckCapitalI]]),
Inverse[\[Sigma] Inverse[CC]-\[DoubleStruckCapitalI]].Inverse[Conjugate[\[Sigma]]Inverse[CC]-\[DoubleStruckCapitalI]]].ConjugateTranspose[VV];


x-xx//Norm;


xx=\[ScriptCapitalI]+
UU.Plus[
(Inverse[Conjugate[\[Sigma]]Inverse[CC]-\[DoubleStruckCapitalI]]+Inverse[\[Sigma] Inverse[CC]-\[DoubleStruckCapitalI]]),
Inverse[\[Sigma] Inverse[CC]-\[DoubleStruckCapitalI]].Inverse[Conjugate[\[Sigma]]Inverse[CC]-\[DoubleStruckCapitalI]]].ConjugateTranspose[VV];

T=Plus[
(Inverse[Conjugate[\[Sigma]]Inverse[CC]-\[DoubleStruckCapitalI]]+Inverse[\[Sigma] Inverse[CC]-\[DoubleStruckCapitalI]]),
Inverse[\[Sigma] Inverse[CC]-\[DoubleStruckCapitalI]].Inverse[Conjugate[\[Sigma]]Inverse[CC]-\[DoubleStruckCapitalI]]];

xx=\[ScriptCapitalI]+UU.T.ConjugateTranspose[VV];

x-xx//Norm;
TT=DiagonalMatrix[((-#(#-2\[Omega]^2))/((#-\[Omega]^2)^2+\[Omega]^2)& /@ Diagonal[CC])];
T-TT//Norm;


xx=\[ScriptCapitalI]+UU.TT.ConjugateTranspose[VV];

x-xx//Norm;

SingularValueList[xx.(\[ScriptCapitalI]-\[ScriptCapitalA]/\[Sigma]).(\[ScriptCapitalI]-\[ScriptCapitalA]/Conjugate[\[Sigma]])];

SingularValueList[(\[ScriptCapitalI]-\[ScriptCapitalA]/\[Sigma]).(\[ScriptCapitalI]-\[ScriptCapitalA]/Conjugate[\[Sigma]])];


Clear[\[Alpha],\[Beta]]


(\[Alpha]-\[Omega]^2 \[Beta])(\[Alpha]-\[Omega]^2 \[Beta])+\[Omega]^2 \[Beta]//Expand


(* ::InheritFromParent:: *)
(*Solve[\[Alpha]^2+\[Beta] \[Omega]^2-2 \[Alpha] \[Beta] \[Omega]^2+\[Beta]^2 \[Omega]^4==0,\[Alpha]]*)


Clear[\[Omega]]


Z1=ArrayFlatten[{
 {-\[Omega] ConjugateTranspose[V].massU.V, 0, 0, 0, 0, ConjugateTranspose[V].S.W, 0, 0ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, 0, 0},
 {0, \[Omega] ConjugateTranspose[V].massU.V, 0, 0, ConjugateTranspose[V].S.W, 0, 0ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, 0, 0, 0},
 {0, 0, -\[Omega] 0\[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, \[CapitalGamma].inverseMassU.S.W, 0, 0\[CapitalGamma].inverseMassU.S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0},
 {0, 0, 0, \[Omega] 0 \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], \[CapitalGamma].inverseMassU.S.W, 0, 0\[CapitalGamma].inverseMassU.S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0},
 {-ConjugateTranspose[W].ConjugateTranspose[S].V, 0, -0ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, ConjugateTranspose[W].massQ.W, \[Omega] ConjugateTranspose[W].massQ.W, 0ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, 0, 0, 0},
 {0, -ConjugateTranspose[W].ConjugateTranspose[S].V, 0, -0ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], -\[Omega] ConjugateTranspose[W].massQ.W, ConjugateTranspose[W].massQ.W, 0, 0ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, 0, 0},
 {-\[CapitalXi].inverseMassQ.ConjugateTranspose[S].V, 0, -0\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, \[CapitalXi].inverseMassQ.massQ.W, 0, 0\[CapitalXi].inverseMassQ.massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0\[Omega] \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0},
 {0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].V, 0, -0\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, \[CapitalXi].inverseMassQ.massQ.W, -\[Omega] 0\[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0\[CapitalXi].inverseMassQ.massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]]},
 {0, 0, 0\[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 0\[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 0, 0, 0, 0\[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, 0, 0, 0},
 {0, 0, 0, 0, 0, 0, 0, 0\[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, 0, 0}
}];
Z2=ArrayFlatten[{
 {-0\[Omega] ConjugateTranspose[V].massU.V, 0, 0, 0, 0, 0ConjugateTranspose[V].S.W, 0, ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, 0, 0},
 {0, 0\[Omega] ConjugateTranspose[V].massU.V, 0, 0, 0ConjugateTranspose[V].S.W, 0, ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, 0, 0, 0},
 {0, 0, -\[Omega] \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0\[CapitalGamma].inverseMassU.S.W, 0, \[CapitalGamma].inverseMassU.S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0\[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0},
 {0, 0, 0, \[Omega] \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0\[CapitalGamma].inverseMassU.S.W, 0, \[CapitalGamma].inverseMassU.S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0\[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0},
 {-0ConjugateTranspose[W].ConjugateTranspose[S].V, 0, -ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0ConjugateTranspose[W].massQ.W, \[Omega] 0ConjugateTranspose[W].massQ.W, ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, 0, 0, 0},
 {0, -0ConjugateTranspose[W].ConjugateTranspose[S].V, 0, -ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], -0\[Omega] ConjugateTranspose[W].massQ.W, 0ConjugateTranspose[W].massQ.W, 0, ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, 0, 0},
 {-0\[CapitalXi].inverseMassQ.ConjugateTranspose[S].V, 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0\[CapitalXi].inverseMassQ.massQ.W, 0, \[CapitalXi].inverseMassQ.massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], \[Omega] \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, 0\[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0},
 {0, -0\[CapitalXi].inverseMassQ.ConjugateTranspose[S].V, 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0\[CapitalXi].inverseMassQ.massQ.W, -\[Omega] \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], \[CapitalXi].inverseMassQ.massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, 0, 0\[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]]},
 {0, 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, 0, 0, 0},
 {0, 0, 0, 0, 0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0, 0, 0, 0}
}];


Z1.Join[\[ScriptL]R,\[ScriptL]I,\[Alpha]R,\[Alpha]I,\[ScriptM]R,\[ScriptM]I,\[Beta]R,\[Beta]I,\[Lambda]R,\[Lambda]I,\[Mu]R,\[Mu]I]+Z2.Join[\[ScriptL]R,\[ScriptL]I,\[Alpha]R,\[Alpha]I,\[ScriptM]R,\[ScriptM]I,\[Beta]R,\[Beta]I,\[Lambda]R,\[Lambda]I,\[Mu]R,\[Mu]I]-BBBBB//Chop//Norm

BB=Join[ConjugateTranspose[V].rhs1I,ConjugateTranspose[V].rhs1R,\[CapitalGamma].inverseMassU.rhs1I,\[CapitalGamma].inverseMassU.rhs1R,ConjugateTranspose[W].rhs2R,ConjugateTranspose[W].rhs2I,\[CapitalXi].inverseMassQ.rhs2R,\[CapitalXi].inverseMassQ.rhs2I];

Z1=ArrayFlatten[
{
 {-\[Omega] ConjugateTranspose[V].massU.V, 0, 0, ConjugateTranspose[V].S.W, 0, 0, 0, 0},
 {0, \[Omega] ConjugateTranspose[V].massU.V, ConjugateTranspose[V].S.W, 0, 0, 0, 0, 0},
 {0, 0, 0, \[CapitalGamma].inverseMassU.S.W, 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0},
 {0, 0, \[CapitalGamma].inverseMassU.S.W, 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0},
 {-ConjugateTranspose[W].ConjugateTranspose[S].V, 0, ConjugateTranspose[W].massQ.W, \[Omega] ConjugateTranspose[W].massQ.W, 0, 0, 0, 0},
 {0, -ConjugateTranspose[W].ConjugateTranspose[S].V, -\[Omega] ConjugateTranspose[W].massQ.W, ConjugateTranspose[W].massQ.W, 0, 0, 0, 0},
 {-\[CapitalXi].inverseMassQ.ConjugateTranspose[S].V, 0, \[CapitalXi].inverseMassQ.massQ.W, 0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0},
 {0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].V, 0, \[CapitalXi].inverseMassQ.massQ.W, 0, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]]}
}];
Z2=ArrayFlatten[{
 {0, 0, 0, ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]]},
 {0, 0, ConjugateTranspose[V].S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0},
 {-\[Omega] \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, \[CapitalGamma].inverseMassU.S.inverseMassQ.ConjugateTranspose[\[CapitalXi]]},
 {0, \[Omega] \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], \[CapitalGamma].inverseMassU.S.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0},
 {-ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0},
 {0, -ConjugateTranspose[W].ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, ConjugateTranspose[W].massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]]},
 {-\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, \[CapitalXi].inverseMassQ.massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]], \[Omega] \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]]},
 {0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], -\[Omega] \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], \[CapitalXi].inverseMassQ.massQ.inverseMassQ.ConjugateTranspose[\[CapitalXi]]}
}];
BBB=BB-Z2.Join[\[Alpha]R,\[Alpha]I,\[Beta]R,\[Beta]I];

Z1.Join[\[ScriptL]R,\[ScriptL]I,\[ScriptM]R,\[ScriptM]I,\[Lambda]R,\[Lambda]I,\[Mu]R,\[Mu]I]+Z2.Join[\[Alpha]R,\[Alpha]I,\[Beta]R,\[Beta]I]-BB//Chop//Norm
Z1.Join[\[ScriptL]R,\[ScriptL]I,\[ScriptM]R,\[ScriptM]I,\[Lambda]R,\[Lambda]I,\[Mu]R,\[Mu]I]-BBB//Chop//Norm






 


BBBBB


Y.Join[\[Beta]I,\[Beta]R,\[Alpha]I,\[Alpha]R]


\[Beta]I-LinearSolve[\[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]],rhs4I]//Chop
\[Beta]R-LinearSolve[\[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]],rhs4R]//Chop
\[Alpha]I-LinearSolve[\[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]],rhs3I]//Chop
\[Alpha]R-LinearSolve[\[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]],rhs3R]//Chop



ZZZZZ.Join[\[ScriptL]R,\[ScriptL]I,\[Alpha]R,\[Alpha]I,\[ScriptM]R,\[ScriptM]I,\[Beta]R,\[Beta]I,\[Lambda]R,\[Lambda]I,\[Mu]R,\[Mu]I]-BBBBB//Chop


\[ScriptL]R//Length
\[ScriptL]I//Length
\[Alpha]R//Length
\[Alpha]I//Length
\[ScriptM]R//Length


BBBBB//Length


ZZZZZ//M

BBBBB=


ConjugateTranspose[V].ConjugateTranspose[\[CapitalGamma]]


\[CapitalGamma].V


-I \[Omega] \[CapitalGamma].inverseMassU.massU.V


rhs1=


ZZ.Join[\[ScriptL],\[Alpha],\[ScriptM],\[Beta],\[Lambda]Soln,\[Mu]Soln]-rhs//Chop

ZZZZ.Join[\[ScriptL],\[Alpha],\[ScriptM],\[Beta],\[Lambda]Soln,\[Mu]Soln]-BBBB//Chop


ZZ[[1]]//Dimensions


Join[\[ScriptL],\[Alpha],\[ScriptM],\[Beta],\[Lambda]Soln,\[Mu]Soln]


\[CapitalXi].W//Chop


ZZZ


ZZZ[[1]].Join[V.\[ScriptL]+inverseMassU.ConjugateTranspose[\[CapitalGamma]].\[Alpha],W.\[ScriptM]+inverseMassQ.ConjugateTranspose[\[CapitalXi]].\[Beta],\[Lambda]Soln,\[Mu]Soln]-rhs//Chop


V.\[ScriptL]+inverseMassU.ConjugateTranspose[\[CapitalGamma]].\[Alpha]-uSoln


\[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]]


rhs4


rhs3


soln//Last
\[Mu]Soln//Last


\[Omega]


\[Omega] massU.uR-S.qI-ConjugateTranspose[\[CapitalGamma]].\[Lambda]I//Chop//Norm
\[Omega] massU.uI+S.qR+ConjugateTranspose[\[CapitalGamma]].\[Lambda]R-rhs1//Chop//Norm
  ConjugateTranspose[S].uI +\[Omega] massQ.qR - massQ.qI - ConjugateTranspose[\[CapitalXi]].\[Mu]I//Chop//Norm
- ConjugateTranspose[S].uR+ massQ.qR+ \[Omega] massQ.qI+  ConjugateTranspose[\[CapitalXi]].\[Mu]R- rhs2//Chop//Norm

-\[CapitalGamma].uI//Chop//Norm
\[CapitalGamma].uR-rhs3//Chop//Norm

\[CapitalXi].qR-rhs4//Chop//Norm
-\[CapitalXi].qI//Chop//Norm


AA


$Assumptions=
{
massU\[Element]Matrices[{3,3},Reals,Symmetric[{1,2}]],
S\[Element]Matrices[{3,12},Reals],
\[CapitalGamma]\[Element]Matrices[{2,3},Reals],
\[CapitalXi]\[Element]Matrices[{2,12},Reals],
massQ\[Element]Matrices[{12,12},Reals,Symmetric[{1,2}] ],
\[Omega]\[Element] Reals
}



AA=TensorReduce[{
 {massU, 0, 1/\[Omega] Transpose[\[CapitalGamma]], 0, 1/\[Omega] S, 0, 0, 0},
 {0, massU, 0, -1/\[Omega] Transpose[\[CapitalGamma]], 0, -1/\[Omega] S, 0, 0},
 {-\[CapitalGamma], 0, 0, 0, 0, 0, 0, 0},
 {0, \[CapitalGamma], 0, 0, 0, 0, 0, 0},
 {0, -Transpose[S], 0, 0, massQ,   \[Omega] massQ, Transpose[\[CapitalXi]], 0},
 {Transpose[S], 0, 0, 0,  \[Omega] massQ,  -massQ, 0, Transpose[\[CapitalXi]]},
 {0, 0, 0, 0, \[CapitalXi], 0, 0, 0},
 {0, 0, 0, 0, 0, \[CapitalXi], 0, 0}
}]
MatrixPlot[AA]

(*AA=TensorReduce[MatrixPower[massQ,0]	  \[Omega] MatrixPower[massQ,0]	-S^\[Transpose]	0	0	0	0	\[CapitalXi]^\[Transpose]
 \[Omega] MatrixPower[massQ,0]	 -MatrixPower[massQ,0]	0	S^\[Transpose]	0	0	\[CapitalXi]^\[Transpose]	0
S	0	0	\[Omega] MatrixPower[massU,0]	\[CapitalGamma]^\[Transpose]	0	0	0
0	-S	\[Omega] MatrixPower[massU,0]	0	0	-\[CapitalGamma]^\[Transpose]	0	0
0	0	\[CapitalGamma]	0	0	0	0	0
0	0	0	-\[CapitalGamma]	0	0	0	0
0	\[CapitalXi]	0	0	0	0	0	0
\[CapitalXi]	0	0	0	0	0	0	0

];*)
(*MatrixPlot[AA]
AA=TensorReduce[massQ	  \[Omega] massQ	\[CapitalXi]^\[Transpose]	0	-S^\[Transpose]	0	0	0
 \[Omega] massQ	 -massQ	0	\[CapitalXi]^\[Transpose]	0	S^\[Transpose]	0	0
\[CapitalXi]	0	0	0	0	0	0	0
0	\[CapitalXi]	0	0	0	0	0	0
S	0	0	0	0	\[Omega] massU	\[CapitalGamma]^\[Transpose]	0
0	-S	0	0	\[Omega] massU	0	0	-\[CapitalGamma]^\[Transpose]
0	0	0	0	\[CapitalGamma]	0	0	0
0	0	0	0	0	-\[CapitalGamma]	0	0

];*)
(*MatrixPlot[AA]

AA=TensorReduce[massQ	\[CapitalXi]^\[Transpose]	  \[Omega] massQ	-S^\[Transpose]	0	0	0	0
\[CapitalXi]	0	0	0	0	0	0	0
 -\[Omega] massQ	0	 massQ	0	-\[CapitalXi]^\[Transpose]	-S^\[Transpose]	0	0
S	0	0	0	0	\[Omega] massU	\[CapitalGamma]^\[Transpose]	0
0	0	\[CapitalXi]	0	0	0	0	0
0	0	-S	\[Omega] massU	0	0	0	-\[CapitalGamma]^\[Transpose]
0	0	0	\[CapitalGamma]	0	0	0	0
0	0	0	0	0	-\[CapitalGamma]	0	0

];
MatrixPlot[AA]*)


AA//MatrixPlot


{r,c}=MinimumBandwidthOrdering[AA,Method->"Sloan"]


|


AA[[r,c]]//MatrixPlot


?*MinimumBandwidthOrdering*


Needs["GraphUtilities`"]


?MinimumBandwidthOrdering


?*Minimum*


?Save





}


{
 {\[Omega] massU, 0, 0, -1/\[Omega] S},
 {0, \[Omega] massU, S, 0},
 {0, ConjugateTranspose[S], \[Omega] massQ, -massQ},
 {-ConjugateTranspose[S], 0, massQ, \[Omega] massQ}
}


Clear[S,massU,massQ]


bb[0]={0,0,-ConjugateTranspose[S]}
BB[0] = {{\[Omega] massU,S,0},{S,\[Omega] massQ,-massQ},{0,massQ,\[Omega] massQ}}


BB[0]-S


BB[0]-Outer[Times,{0,0,-S},{0,0,-ConjugateTranspose[S]}]


?Outer


BB[0]


Clear[bb]


Clear


BB[0]


B[


Outer[Times,bb,bb]


?Outer


AA





\[Omega]=1


Clear[\[Omega]]


$Assumptions=
{
massU\[Element]Matrices[{3,3},Reals,Symmetric[{1,2}]],
S\[Element]Matrices[{3,12},Reals],
\[CapitalGamma]\[Element]Matrices[{2,3},Reals],
\[CapitalXi]\[Element]Matrices[{2,12},Reals],
massQ\[Element]Matrices[{12,12},Reals,Symmetric[{1,2}] ],
\[Omega]\[Element] Reals
}


AA={
 { massU, 0, -1/\[Omega] S, 0},
 {0,  massU, 0, 1/\[Omega] S},
 {0, -Transpose[S], massQ, -\[Omega] massQ},
 {-Transpose[S], 0, \[Omega] massQ, massQ}
}


A[2]//MatrixPlot


Diagonal[A[2]]



myDot[{},{}] := 0
myDot[x_List,y_List] := Dot[First[x],First[y]]+myDot[Rest[x],Rest[y]]
myMatMat[x_,y_] := Table[myDot[x[[i,All]],y[[All,j]]],{i,1,8},{j,1,8}]
Clear[L,A,Ahalf,U]

A[0]=AA;
L[i_] := L[i]=Table[Boole[ii==jj]+If[jj==i,If[ii>i,-A[i-1][[ii,jj]].Inverse[A[i-1][[i,i]]],0],0],{ii,1,8},{jj,1,8}]
U[i_] := U[i]=
Table[
  Boole[ii==jj]
   +If[ii==i,
      If[jj>i,-Inverse[Ahalf[i][[i,i]]].Ahalf[i][[ii,jj]],0]
    ,0]
,{ii,1,8},{jj,1,8}]
Ahalf[i_] := Ahalf[i]=TensorReduce[myMatMat[L[i],A[i-1]]]
A[i_]:=A[i]=TensorReduce[myMatMat[Ahalf[i],U[i]]]


(*A[2]=myMatMat[L[2],A[1]];
A[3]=myMatMat[L[3],A[2]];
A[i_] := A[i]=myMatMat[L[i],A[i-1]];*)

Unprotect[Dot]

Dot[0,x_] := 0
Dot[x_,0] := 0
Dot[x_,1] := x
Dot[1,x_] := x
Dot[x_,-y_] := -Dot[x,y]
Dot[-x_,y_] := -Dot[x,y]
(*Dot[x_,inv[x_] ] := 1
Dot[inv[x_],x_ ] := 1*)
Dot[\[Omega],x_] := \[Omega] x
Dot[x_,(\[Omega])^-1 y_] := (\[Omega])^-1 Dot[x,y]
Dot[(\[Omega])^-1 x_, y_] := (\[Omega])^-1 Dot[x,y]
Dot[\[Omega] x_,y_] := \[Omega] Dot[x,y]




A[0]
Ahalf[1]
A[1]
Ahalf[2]
A[2]
Ahalf[3]
A[3]
Ahalf[4]
A[4]//MatrixPlot


A[5]//Diagonal












(* ::InheritFromParent:: *)
(*{{massU,0,-massU.Inverse[massU].(TensorTranspose[\[CapitalGamma],{2,1}]/\[Omega])+TensorTranspose[\[CapitalGamma],{2,1}]/\[Omega],0,S/\[Omega]-massU.Inverse[massU].(S/\[Omega]),0,0,0},{0,massU,0,-(TensorTranspose[\[CapitalGamma],{2,1}]/\[Omega]),0,-(S/\[Omega]),0,0},{-\[CapitalGamma]+\[CapitalGamma].Inverse[massU].massU,0,\[CapitalGamma].Inverse[massU].(TensorTranspose[\[CapitalGamma],{2,1}]/\[Omega])-(-\[CapitalGamma]+\[CapitalGamma].Inverse[massU].massU).Inverse[massU].(TensorTranspose[\[CapitalGamma],{2,1}]/\[Omega]),0,\[CapitalGamma].Inverse[massU].(S/\[Omega])-(-\[CapitalGamma]+\[CapitalGamma].Inverse[massU].massU).Inverse[massU].(S/\[Omega]),0,0,0},{0,\[CapitalGamma],0,0,0,0,0,0},{0,-TensorTranspose[S,{2,1}],0,0,massQ,massQ \[Omega],TensorTranspose[\[CapitalXi],{2,1}],0},{-TensorTranspose[S,{2,1}].Inverse[massU].massU+TensorTranspose[S,{2,1}],0,-TensorTranspose[S,{2,1}].Inverse[massU].(TensorTranspose[\[CapitalGamma],{2,1}]/\[Omega])-(-TensorTranspose[S,{2,1}].Inverse[massU].massU+TensorTranspose[S,{2,1}]).Inverse[massU].(TensorTranspose[\[CapitalGamma],{2,1}]/\[Omega]),0,massQ \[Omega]-TensorTranspose[S,{2,1}].Inverse[massU].(S/\[Omega])-(-TensorTranspose[S,{2,1}].Inverse[massU].massU+TensorTranspose[S,{2,1}]).Inverse[massU].(S/\[Omega]),-massQ,0,TensorTranspose[\[CapitalXi],{2,1}]},{0,0,0,0,\[CapitalXi],0,0,0},{0,0,0,0,0,\[CapitalXi],0,0}}*)


1


A[5]//MatrixPlot



(* ::InheritFromParent:: *)
(*{{massU,0,0,0},{0,massU,0,0},{0,0,massQ,-massQ \[Omega]+TensorTranspose[S,{2,1}].MatrixPower[massU,-1].S/\[Omega]},{0,0,massQ \[Omega]-Transpose[S,{2,1}].MatrixPower[massU,-1].S/\[Omega],massQ}}*)


-(-massQ \[Omega]+TensorTranspose[S,{2,1}].MatrixPower[massU,-1].S/\[Omega]).MatrixPower[massQ,-1].(-massQ \[Omega]+TensorTranspose[S,{2,1}].MatrixPower[massU,-1].S/\[Omega])//TensorReduce


A[3]


(* ::InheritFromParent:: *)
(**)


-S.TensorTranspose[S,{2,1}]-\[Omega]^2 ((1+\[Omega]^2) MatrixPower[S.TensorTranspose[S,{2,1}],-1]-2 MatrixPower[S.TensorTranspose[S,{2,1}],0])


-s^2-\[Omega]^2 ((1+\[Omega]^2) s^-2-2)//Expand//Factor


A[0]//MatrixPlot
(*A[1]//MatrixPlot
A[2]//MatrixPlot
A[3]//MatrixPlot*)


(A[3]//Diagonal)/. \[Omega]->1


TensorReduce[
myMatMat[
L[3],TensorReduce[
myMatMat[L[2],Simplify@TensorReduce[
myMatMat[L[1],TensorReduce[myMatMat[
TensorReduce[myMatMat[TensorReduce[myMatMat[AA,U[1]]],U[2]]],U[3]]]]]]]]]





L[4]


A[0]//MatrixPlot
Ahalf[1];
A[1]//MatrixPlot
Ahalf[2];
A[2]//MatrixPlot
Ahalf[3];
A[3]//MatrixPlot
Ahalf[4];
A[4]//MatrixPlot



subs={S.MatrixPower[massQ,-1].TensorTranspose[S,{2,1}]:> ZZ,}


A[3]//TensorReduce


(* ::InheritFromParent:: *)
(*{massQ,-massQ (1+\[Omega]^2),S.MatrixPower[massQ,-1].TensorTranspose[S,{2,1}]/(1+\[Omega]^2),2 massU \[Omega]^2-(\[Omega]^2+\[Omega]^4) massU.MatrixPower[S.MatrixPower[massQ,-1].TensorTranspose[S,{2,1}],-1].massU-S.MatrixPower[massQ,-1].TensorTranspose[S,{2,1}],-(1+\[Omega]^2) \[CapitalGamma].MatrixPower[S.MatrixPower[massQ,-1].TensorTranspose[S,{2,1}],-1].TensorTranspose[\[CapitalGamma],{2,1}],0,(1/(1+\[Omega]^2))(\[CapitalXi].MatrixPower[massQ,-1].TensorTranspose[\[CapitalXi],{2,1}]+\[Omega]^2 \[CapitalXi].MatrixPower[massQ,-1].TensorTranspose[S,{2,1}].MatrixPower[S.MatrixPower[massQ,-1].TensorTranspose[S,{2,1}],-1].S.MatrixPower[massQ,-1].TensorTranspose[\[CapitalXi],{2,1}]),(1/(1+\[Omega]^2))(-\[CapitalXi].MatrixPower[massQ,-1].TensorTranspose[\[CapitalXi],{2,1}]+\[CapitalXi].MatrixPower[massQ,-1].TensorTranspose[S,{2,1}].MatrixPower[S.MatrixPower[massQ,-1].TensorTranspose[S,{2,1}],-1].S.MatrixPower[massQ,-1].TensorTranspose[\[CapitalXi],{2,1}])}*)


(* ::InheritFromParent:: *)
(*(\[Omega]+\[Omega]^3) massQ.MatrixPower[-massQ (1+\[Omega]^2),-1].TensorTranspose[S,{2,1}]+\[Omega] TensorTranspose[S,{2,1}]*)


(* ::InheritFromParent:: *)
(*(\[Omega]+\[Omega]^3) massQ.MatrixPower[-massQ (1+\[Omega]^2),-1].TensorTranspose[S,{2,1}]+\[Omega] TensorTranspose[S,{2,1}]//TensorReduce//Simplify*)


-((\[Omega]^2 S.MatrixPower[massQ,-1].TensorTranspose[S,{2,1}])/(1+\[Omega]^2))+S.MatrixPower[massQ,-1].Transpose[S,{2,1}]//Together


Map[TensorReduce,Ahalf[1],{2}]


Transpose[S]//TensorReduce


A[0]//MatrixPlot
Ahalf[1]//MatrixPlot
A[1]//MatrixPlot
Ahalf[2]//MatrixPlot;
A[2]//MatrixPlot
Ahalf[3]//MatrixPlot
A[3]//MatrixPlot



A[1][[1,2]]


Ahalf[4]//MatrixPlot
A[4]//MatrixPlot
Ahalf[5]//MatrixPlot
A[5]//MatrixPlot


A[0]//MatrixPlot



A


(* ::InheritFromParent:: *)
(*{{massQ,massQ \[Omega],-ConjugateTranspose[S],0,0,0,ConjugateTranspose[\[CapitalXi]],0},{0,massQ+massQ \[Omega]^2,-\[Omega] massQ.Inverse[massQ].ConjugateTranspose[S],-ConjugateTranspose[S],0,0,\[Omega] massQ.Inverse[massQ].ConjugateTranspose[\[CapitalXi]],-ConjugateTranspose[\[CapitalXi]]},{0,-S \[Omega],S.Inverse[massQ].ConjugateTranspose[S],massU \[Omega],ConjugateTranspose[\[CapitalGamma]],0,-S.Inverse[massQ].ConjugateTranspose[\[CapitalXi]],0},{0,-S,massU \[Omega],0,0,-ConjugateTranspose[\[CapitalGamma]],0,0},{0,0,\[CapitalGamma],0,0,0,0,0},{0,0,0,-\[CapitalGamma],0,0,0,0},{0,-\[CapitalXi] \[Omega],\[CapitalXi].Inverse[massQ].ConjugateTranspose[S],0,0,0,-\[CapitalXi].Inverse[massQ].ConjugateTranspose[\[CapitalXi]],0},{0,\[CapitalXi],0,0,0,0,0,0}}*)


(* ::InheritFromParent:: *)
(*{{massQ,massQ \[Omega],-ConjugateTranspose[S],0,0,0,ConjugateTranspose[\[CapitalXi]],0},{-massQ \[Omega]+\[Omega] massQ.Inverse[massQ].massQ,massQ+\[Omega]^2 massQ.Inverse[massQ].massQ,-\[Omega] massQ.Inverse[massQ].ConjugateTranspose[S],-ConjugateTranspose[S],0,0,\[Omega] massQ.Inverse[massQ].ConjugateTranspose[\[CapitalXi]],-ConjugateTranspose[\[CapitalXi]]},{S-S.Inverse[massQ].massQ,-\[Omega] S.Inverse[massQ].massQ,S.Inverse[massQ].ConjugateTranspose[S],massU \[Omega],ConjugateTranspose[\[CapitalGamma]],0,-S.Inverse[massQ].ConjugateTranspose[\[CapitalXi]],0},{0,-S,massU \[Omega],0,0,-ConjugateTranspose[\[CapitalGamma]],0,0},{0,0,\[CapitalGamma],0,0,0,0,0},{0,0,0,-\[CapitalGamma],0,0,0,0},{\[CapitalXi]-\[CapitalXi].Inverse[massQ].massQ,-\[Omega] \[CapitalXi].Inverse[massQ].massQ,\[CapitalXi].Inverse[massQ].ConjugateTranspose[S],0,0,0,-\[CapitalXi].Inverse[massQ].ConjugateTranspose[\[CapitalXi]],0},{0,\[CapitalXi],0,0,0,0,0,0}}*)


-\[Omega] \[CapitalXi].Inverse[massQ].massQ//TensorReduce


S**inv[massQ]**ConjugateTranspose[S]-(\[Omega]^2 S**inv[massQ]**ConjugateTranspose[S])/(1+\[Omega]^2)//Together


?*Symm*




A[3]


(* ::InheritFromParent:: *)
(*{massQ+massQ \[Omega]^2:>(1+\[Omega]^2) massQ,inv[massQ (1+\[Omega]^2)]:>inv[massQ]}*)


inv[massQ (1+\[Omega]^2)]


A[0]//MatrixPlot
Ahalf[1]//MatrixPlot
A[1]//MatrixPlot
Ahalf[2]//MatrixPlot
A[2]//MatrixPlot
(*Ahalf[3]//MatrixPlot
A[3]//MatrixPlot
Ahalf[4]//MatrixPlot
A[4]//MatrixPlot
Ahalf[5]//MatrixPlot
A[5]//MatrixPlot
Ahalf[6]//MatrixPlot
A[6]//MatrixPlot
Ahalf[7]//MatrixPlot
A[7]//MatrixPlot*)


A[0];
A[1];
A[2];
Ahalf[3]


A[2][[7,7]]



(* ::InheritFromParent:: *)
(*{{massQ,0,0,0,0,0,0,0},{0,massQ (1+\[Omega]^2),-\[Omega] ConjugateTranspose[S],-ConjugateTranspose[S],0,0,\[Omega] ConjugateTranspose[\[CapitalXi]],-ConjugateTranspose[\[CapitalXi]]},{0,0,S**inv[massQ]**ConjugateTranspose[S]-\[Omega]^2 S**inv[massQ (1+\[Omega]^2)]**ConjugateTranspose[S],massU \[Omega]-\[Omega] S**inv[massQ (1+\[Omega]^2)]**ConjugateTranspose[S],ConjugateTranspose[\[CapitalGamma]],0,-S**inv[massQ]**ConjugateTranspose[\[CapitalXi]]+\[Omega]^2 S**inv[massQ (1+\[Omega]^2)]**ConjugateTranspose[\[CapitalXi]],-\[Omega] S**inv[massQ (1+\[Omega]^2)]**ConjugateTranspose[\[CapitalXi]]},{0,0,massU \[Omega]-\[Omega] S**inv[massQ (1+\[Omega]^2)]**ConjugateTranspose[S],-S**inv[massQ (1+\[Omega]^2)]**ConjugateTranspose[S],0,-ConjugateTranspose[\[CapitalGamma]],\[Omega] S**inv[massQ (1+\[Omega]^2)]**ConjugateTranspose[\[CapitalXi]],-S**inv[massQ (1+\[Omega]^2)]**ConjugateTranspose[\[CapitalXi]]},{0,0,\[CapitalGamma],0,0,0,0,0},{0,0,0,-\[CapitalGamma],0,0,0,0},{0,0,\[CapitalXi]**inv[massQ]**ConjugateTranspose[S]-\[Omega]^2 \[CapitalXi]**inv[massQ (1+\[Omega]^2)]**ConjugateTranspose[S],-\[Omega] \[CapitalXi]**inv[massQ (1+\[Omega]^2)]**ConjugateTranspose[S],0,0,-\[CapitalXi]**inv[massQ]**ConjugateTranspose[\[CapitalXi]]+\[Omega]^2 \[CapitalXi]**inv[massQ (1+\[Omega]^2)]**ConjugateTranspose[\[CapitalXi]],-\[Omega] \[CapitalXi]**inv[massQ (1+\[Omega]^2)]**ConjugateTranspose[\[CapitalXi]]},{0,0,\[Omega] \[CapitalXi]**inv[massQ (1+\[Omega]^2)]**ConjugateTranspose[S],\[CapitalXi]**inv[massQ (1+\[Omega]^2)]**ConjugateTranspose[S],0,0,-\[Omega] \[CapitalXi]**inv[massQ (1+\[Omega]^2)]**ConjugateTranspose[\[CapitalXi]],\[CapitalXi]**inv[massQ (1+\[Omega]^2)]**ConjugateTranspose[\[CapitalXi]]}}*)


(* ::InheritFromParent:: *)
(*{{massQ,massQ \[Omega],-ConjugateTranspose[S],0,0,0,ConjugateTranspose[\[CapitalXi]],0},{0,massQ+massQ \[Omega]^2,-\[Omega] ConjugateTranspose[S],-ConjugateTranspose[S],0,0,\[Omega] ConjugateTranspose[\[CapitalXi]],-ConjugateTranspose[\[CapitalXi]]},{0,-S \[Omega],S**inv[massQ]**ConjugateTranspose[S],massU \[Omega],ConjugateTranspose[\[CapitalGamma]],0,-S**inv[massQ]**ConjugateTranspose[\[CapitalXi]],0},{0,-S,massU \[Omega],0,0,-ConjugateTranspose[\[CapitalGamma]],0,0},{0,0,\[CapitalGamma],0,0,0,0,0},{0,0,0,-\[CapitalGamma],0,0,0,0},{0,-\[CapitalXi] \[Omega],\[CapitalXi]**inv[massQ]**ConjugateTranspose[S],0,0,0,-\[CapitalXi]**inv[massQ]**ConjugateTranspose[\[CapitalXi]],0},{0,\[CapitalXi],0,0,0,0,0,0}}*)


massQ+massQ \[Omega]^2:> (1+\[Omega]^2)massQ


L[1]
Ahalf[1]
U[1]


A[7][[1]]


L[1]//MatrixPlot


U[1]


A[1/2]


L[1]


A[0]//MatrixPlot
A[1]//MatrixPlot
A[2]//MatrixPlot
A[3]//MatrixPlot
A[4]//MatrixPlot
A[5]//MatrixPlot
A[6]//MatrixPlot
A[7]//MatrixPlot


A[7]//MatrixPlot


\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(7\)]\(L[i]\)\)//MatrixPlot


A[1]
A[2]


A[3]


A


CC=CholeskyDecomposition[massU];


A[0] := {{MM,SS},{-ConjugateTranspose[SS],0}};


L[1] := {{CCC,0},{-inv[CCC]**ConjugateTranspose[SS],1}}


A





ConjugateTranspose[CC].CC-massU


A[2][[3,3]]





L[2]//MatrixPlot



Unprotect[NonCommutativeMultiply]
NonCommutativeMultiply[0,x_] := 0
NonCommutativeMultiply[x_,0] := 0
NonCommutativeMultiply[x_,1] := x
NonCommutativeMultiply[1,x_] := x
NonCommutativeMultiply[x_,-y_] := -NonCommutativeMultiply[x,y]
NonCommutativeMultiply[-x_,y_] := -NonCommutativeMultiply[x,y]
NonCommutativeMultiply[x_,inv[x_] ] := 1
NonCommutativeMultiply[x_,\[Omega] y_]:= \[Omega] NonCommutativeMultiply[x, y]


L[1]//MatrixPlot





myMatMat[L[1],A[0]]


Table[myDot[L[1][[i,All]],A[0][[All,j]]],{i,1,8},{j,1,8}]//Expand


myDot[{},{}] := 0
myDot[x_List,y_List] := NonCommutativeMultiply[First[x],First[y]]+myDot[Rest[x],Rest[y]]
myMatMat[x_,y_] := Table[myDot[x[[i,All]],y[[All,j]]],{i,1,8},{j,1,8}]//Expand


myDot[x_List


myDot[


myDot[{a,b,c},{d,e,f}]


m


Outer[L[1],A[0],Times]


?Dot


?Outer


?Outer





L[1]//MatrixPlot


?*Mult*


L[1]


AA//Normal


CholeskyDecomposition[massU]


\[DoubleStruckCapitalA]=SparseArray@ArrayFlatten[{
 {\[Omega] massU, 0, 0, -S},
 {0, \[Omega] massU, S, 0},
 {0, ConjugateTranspose[S], \[Omega] massQ, -massQ},
 {-ConjugateTranspose[S], 0, massQ, \[Omega] massQ}
}];
\[DoubleStruckCapitalB]=SparseArray@ArrayFlatten[{
 {0, -ConjugateTranspose[\[CapitalGamma]], 0, 0},
 {ConjugateTranspose[\[CapitalGamma]], 0, 0, 0},
 {0, 0, 0, -ConjugateTranspose[\[CapitalXi]]},
 {0, 0, ConjugateTranspose[\[CapitalXi]], 0}
}];


RowReduce[SparseArray@ArrayFlatten[{
 {\[Omega] massU, 0, 0, -S, \[ScriptCapitalI]u, 0, 0, 0},
 {0, \[Omega] massU, S, 0, 0, \[ScriptCapitalI]u, 0, 0},
 {0, ConjugateTranspose[S], \[Omega] massQ, -massQ, 0, 0, \[ScriptCapitalI]q, 0},
 {-ConjugateTranspose[S], 0, massQ, \[Omega] massQ, 0, 0, 0, \[ScriptCapitalI]q}
}]]//MatrixPlot

RowReduce[SparseArray@ArrayFlatten[{
 {\[ScriptCapitalI]u, 0, 0, -1/\[Omega] inverseMassU.S, 1/\[Omega] inverseMassU, 0, 0, 0},
 {0, \[Omega] massU, S, 0, 0, \[ScriptCapitalI]u, 0, 0},
 {0, ConjugateTranspose[S], \[Omega] massQ, -massQ, 0, 0, \[ScriptCapitalI]q, 0},
 {-ConjugateTranspose[S], 0, massQ, \[Omega] massQ, 0, 0, 0, \[ScriptCapitalI]q}
}]]//MatrixPlot

RowReduce[SparseArray@ArrayFlatten[{
 {\[ScriptCapitalI]u, 0, 0, -1/\[Omega] inverseMassU.S, 1/\[Omega] inverseMassU, 0, 0, 0},
 {0, \[ScriptCapitalI]u, 1/\[Omega] inverseMassU.S, 0, 0, 1/\[Omega] inverseMassU, 0, 0},
 {0, ConjugateTranspose[S], \[Omega] massQ, -massQ, 0, 0, \[ScriptCapitalI]q, 0},
 {0, 0, massQ, -1/\[Omega] ConjugateTranspose[S].inverseMassU.S+\[Omega] massQ, 1/\[Omega] ConjugateTranspose[S].inverseMassU, 0, 0, \[ScriptCapitalI]q}
}]]//MatrixPlot

RowReduce[SparseArray@ArrayFlatten[{
 {\[ScriptCapitalI]u, 0, 0, -1/\[Omega] inverseMassU.S, 1/\[Omega] inverseMassU, 0, 0, 0},
 {0, \[ScriptCapitalI]u, 1/\[Omega] inverseMassU.S, 0, 0, 1/\[Omega] inverseMassU, 0, 0},
 {0, 0, 1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ, massQ, 0, 1/\[Omega] ConjugateTranspose[S].inverseMassU, -\[ScriptCapitalI]q, 0},
 {0, 0, massQ, -1/\[Omega] ConjugateTranspose[S].inverseMassU.S+\[Omega] massQ, 1/\[Omega] ConjugateTranspose[S].inverseMassU, 0, 0, \[ScriptCapitalI]q}
}]]//MatrixPlot

RowReduce[SparseArray@ArrayFlatten[{
 {\[ScriptCapitalI]u, 0, 0, -1/\[Omega] inverseMassU.S, 1/\[Omega] inverseMassU, 0, 0, 0},
 {0, \[ScriptCapitalI]u, 1/\[Omega] inverseMassU.S, 0, 0, 1/\[Omega] inverseMassU, 0, 0},
 {0, 0, 1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ, massQ, 0, 1/\[Omega] ConjugateTranspose[S].inverseMassU, -\[ScriptCapitalI]q, 0},
 {0, 0, \[ScriptCapitalI]q, inverseMassQ.(-1/\[Omega] ConjugateTranspose[S].inverseMassU.S+\[Omega] massQ), 1/\[Omega] inverseMassQ.ConjugateTranspose[S].inverseMassU, 0, 0, inverseMassQ}
}]]//MatrixPlot

RowReduce[SparseArray@ArrayFlatten[{
 {\[ScriptCapitalI]u, 0, 0, -1/\[Omega] inverseMassU.S, 1/\[Omega] inverseMassU, 0, 0, 0},
 {0, \[ScriptCapitalI]u, 1/\[Omega] inverseMassU.S, 0, 0, 1/\[Omega] inverseMassU, 0, 0},
 {0, 0, \[ScriptCapitalI]q, -inverseMassQ.(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ), 1/\[Omega] inverseMassQ.ConjugateTranspose[S].inverseMassU, 0, 0, inverseMassQ},
 {0, 0, 1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ, massQ, 0, 1/\[Omega] ConjugateTranspose[S].inverseMassU, -\[ScriptCapitalI]q, 0}
}]]//MatrixPlot

RowReduce[SparseArray@ArrayFlatten[{
 {\[ScriptCapitalI]u, 0, 0, -1/\[Omega] inverseMassU.S, 1/\[Omega] inverseMassU, 0, 0, 0},
 {0, \[ScriptCapitalI]u, 1/\[Omega] inverseMassU.S, 0, 0, 1/\[Omega] inverseMassU, 0, 0},
 {0, 0, \[ScriptCapitalI]q, -inverseMassQ.(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ), 1/\[Omega] inverseMassQ.ConjugateTranspose[S].inverseMassU, 0, 0, inverseMassQ},
 {0, 0, 0, (1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ).inverseMassQ.(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ)+massQ, -(1/\[Omega])(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ).inverseMassQ.ConjugateTranspose[S].inverseMassU, 1/\[Omega] ConjugateTranspose[S].inverseMassU, -\[ScriptCapitalI]q, -(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ).inverseMassQ}
}]]//MatrixPlot

\[Kappa]=(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ).inverseMassQ.(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ)+massQ;
inv\[Kappa]=Inverse[\[Kappa]];

RowReduce[SparseArray@ArrayFlatten[{
 {\[ScriptCapitalI]u, 0, 0, -1/\[Omega] inverseMassU.S, 1/\[Omega] inverseMassU, 0, 0, 0},
 {0, \[ScriptCapitalI]u, 1/\[Omega] inverseMassU.S, 0, 0, 1/\[Omega] inverseMassU, 0, 0},
 {0, 0, \[ScriptCapitalI]q, -inverseMassQ.(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ), 1/\[Omega] inverseMassQ.ConjugateTranspose[S].inverseMassU, 0, 0, inverseMassQ},
 {0, 0, 0, \[Kappa], -(1/\[Omega])(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ).inverseMassQ.ConjugateTranspose[S].inverseMassU, 1/\[Omega] ConjugateTranspose[S].inverseMassU, -\[ScriptCapitalI]q, -(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ).inverseMassQ}
}]]//MatrixPlot

RowReduce[SparseArray@ArrayFlatten[{
 {\[ScriptCapitalI]u, 0, 0, -1/\[Omega] inverseMassU.S, 1/\[Omega] inverseMassU, 0, 0, 0},
 {0, \[ScriptCapitalI]u, 1/\[Omega] inverseMassU.S, 0, 0, 1/\[Omega] inverseMassU, 0, 0},
 {0, 0, \[ScriptCapitalI]q, -inverseMassQ.(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ), 1/\[Omega] inverseMassQ.ConjugateTranspose[S].inverseMassU, 0, 0, inverseMassQ},
 {0, 0, 0, \[ScriptCapitalI]q, inv\[Kappa].(-(1/\[Omega])(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ).inverseMassQ.ConjugateTranspose[S].inverseMassU), 1/\[Omega] inv\[Kappa].ConjugateTranspose[S].inverseMassU, -inv\[Kappa], -inv\[Kappa].(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ).inverseMassQ}
}]]//MatrixPlot



RowReduce[SparseArray@ArrayFlatten[{
 {\[ScriptCapitalI]u, 0, 0, -1/\[Omega] inverseMassU.S, 1/\[Omega] inverseMassU, 0, 0, 0},
 {0, \[ScriptCapitalI]u, 1/\[Omega] inverseMassU.S, 0, 0, 1/\[Omega] inverseMassU, 0, 0},
 {0, 0, \[ScriptCapitalI]q, 0, 1/\[Omega] inverseMassQ.ConjugateTranspose[S].inverseMassU+(inverseMassQ.(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ)).(inv\[Kappa].(-(1/\[Omega])(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ).inverseMassQ.ConjugateTranspose[S].inverseMassU)), (inverseMassQ.(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ)).(1/\[Omega] inv\[Kappa].ConjugateTranspose[S].inverseMassU), (inverseMassQ.(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ)).(-inv\[Kappa]), inverseMassQ+(inverseMassQ.(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ)).(-inv\[Kappa].(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ).inverseMassQ)},
 {0, 0, 0, \[ScriptCapitalI]q, inv\[Kappa].(-(1/\[Omega])(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ).inverseMassQ.ConjugateTranspose[S].inverseMassU), 1/\[Omega] inv\[Kappa].ConjugateTranspose[S].inverseMassU, -inv\[Kappa], -inv\[Kappa].(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ).inverseMassQ}
}]]//MatrixPlot


RowReduce[SparseArray@ArrayFlatten[{
 {\[ScriptCapitalI]u, 0, 0, 0, 1/\[Omega] inverseMassU+(1/\[Omega] inverseMassU.S).(inv\[Kappa].(-(1/\[Omega])(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ).inverseMassQ.ConjugateTranspose[S].inverseMassU)), +(1/\[Omega] inverseMassU.S).(1/\[Omega] inv\[Kappa].ConjugateTranspose[S].inverseMassU), +(1/\[Omega] inverseMassU.S).(-inv\[Kappa]), (1/\[Omega] inverseMassU.S).(-inv\[Kappa].(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ).inverseMassQ)},
 {0, \[ScriptCapitalI]u, 0, 0, -(1/\[Omega] inverseMassU.S).(1/\[Omega] inverseMassQ.ConjugateTranspose[S].inverseMassU+(inverseMassQ.(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ)).(inv\[Kappa].(-(1/\[Omega])(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ).inverseMassQ.ConjugateTranspose[S].inverseMassU))), 1/\[Omega] inverseMassU-(1/\[Omega] inverseMassU.S).((inverseMassQ.(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ)).(1/\[Omega] inv\[Kappa].ConjugateTranspose[S].inverseMassU)), -(1/\[Omega] inverseMassU.S).((inverseMassQ.(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ)).(-inv\[Kappa])), -(1/\[Omega] inverseMassU.S).(inverseMassQ+(inverseMassQ.(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ)).(-inv\[Kappa].(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ).inverseMassQ))},
 {0, 0, \[ScriptCapitalI]q, 0, 1/\[Omega] inverseMassQ.ConjugateTranspose[S].inverseMassU+(inverseMassQ.(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ)).(inv\[Kappa].(-(1/\[Omega])(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ).inverseMassQ.ConjugateTranspose[S].inverseMassU)), (inverseMassQ.(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ)).(1/\[Omega] inv\[Kappa].ConjugateTranspose[S].inverseMassU), (inverseMassQ.(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ)).(-inv\[Kappa]), inverseMassQ+(inverseMassQ.(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ)).(-inv\[Kappa].(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ).inverseMassQ)},
 {0, 0, 0, \[ScriptCapitalI]q, inv\[Kappa].(-(1/\[Omega])(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ).inverseMassQ.ConjugateTranspose[S].inverseMassU), 1/\[Omega] inv\[Kappa].ConjugateTranspose[S].inverseMassU, -inv\[Kappa], -inv\[Kappa].(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ).inverseMassQ}
}]]//MatrixPlot







1/\[Omega] inverseMassU+(1/\[Omega] inverseMassU.S).(inv\[Kappa].(-(1/\[Omega])(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ).inverseMassQ.ConjugateTranspose[S].inverseMassU))


(1/\[Omega] inverseMassU.S)


\[Omega]=10


{
 {\[Omega] massU, 0, 0, -S},
 {0, \[Omega] massU, S, 0},
 {0, ConjugateTranspose[S], \[Omega] massQ, -massQ},
 {-ConjugateTranspose[S], 0, massQ, \[Omega] massQ}
}


AA=SparseArray@ArrayFlatten[{
 {massQ, -\[Omega] massQ, 0, -ConjugateTranspose[S]},
 {\[Omega] massQ, massQ, -ConjugateTranspose[S], 0},
 {0, S, 0, \[Omega] massU},
 {-S, 0, \[Omega] massU, 0}
}];


$Assumptions=
{
massU\[Element]Matrices[{3,3},Reals,Symmetric[{1,2}]],
S\[Element]Matrices[{3,12},Reals],
\[CapitalGamma]\[Element]Matrices[{2,3},Reals],
\[CapitalXi]\[Element]Matrices[{2,12},Reals],
massQ\[Element]Matrices[{12,12},Reals,Symmetric[{1,2}] ],
invMassQ\[Element]Matrices[{12,12},Reals,Symmetric[{1,2}] ],
\[Omega]\[Element] Reals

}


A=ArrayFlatten[{
 {massQ, -\[Omega] massQ},
 {\[Omega] massQ, massQ}
}]
B=ArrayFlatten[{
 {0, -Transpose[S]},
 {-Transpose[S], 0}
}]
CC=ArrayFlatten[{
 {0, S},
 {-S, 0}
}]
DD=ArrayFlatten[{
 {0, \[Omega] massU},
 {\[Omega] massU, 0}
}]
invA=ArrayFlatten[Inverse[{{1,-\[Omega]},{\[Omega],1}}]invMassQ]



CCinvAdotBB={{-1,\[Omega]},{\[Omega],1}} 1/(1+\[Omega]^2) S invMassQ Transpose[S];
CCinvA={{-\[Omega],1},{-1,-\[Omega]}} 1/(1+\[Omega]^2) S invMassQ
DDmCCinvADotBB=DD-CCinvAdotBB
invAdotBB={{-\[Omega],-1},{-1,\[Omega]}} 1/(1+\[Omega]^2) invMassQ Transpose[S]








(* ::InheritFromParent:: *)
(*{{-((invMassQ \[Omega] Transpose[S])/(1+\[Omega]^2)),-((invMassQ Transpose[S])/(1+\[Omega]^2))},{-((invMassQ Transpose[S])/(1+\[Omega]^2)),(invMassQ \[Omega] Transpose[S])/(1+\[Omega]^2)}}*)


{{1,0,-((invMassQ \[Omega] Transpose[S])/(1+\[Omega]^2)),-((invMassQ Transpose[S])/(1+\[Omega]^2))},{0,1,-((invMassQ Transpose[S])/(1+\[Omega]^2)),(invMassQ \[Omega] Transpose[S])/(1+\[Omega]^2)},{0,0,1,0},{0,0,0,1}}//Inverse


invA.B//Expand//Hash


CC.invA.B//Expand//Hash


({{-\[Omega],1},{-1,-\[Omega]}} 1/(1+\[Omega]^2) S invMassQ).{{0,Transpose[S]},{-S\[Transpose],0}}//Expand


invA.CC//Expand//Hash


A//InputForm





?*Inverse*


Inverse[A]


invA


myDot[


CC.invA.B


CC.invA


invA=Inverse[{{1,\[Omega]},{\[Omega],-1}}]invMassQ;


AA//MatrixPlot


AA=SparseArray@ArrayFlatten[{
 {\[Omega] massU, 0, 0, -S, 0, -ConjugateTranspose[\[CapitalGamma]], 0, 0},
 {0, \[Omega] massU, S, 0, ConjugateTranspose[\[CapitalGamma]], 0, 0, 0},
 {0, ConjugateTranspose[S], \[Omega] massQ, -massQ, 0, 0, 0, -ConjugateTranspose[\[CapitalXi]]},
 {-ConjugateTranspose[S], 0, massQ, \[Omega] massQ, 0, 0, ConjugateTranspose[\[CapitalXi]], 0},
 {0, -\[CapitalGamma], 0, 0, 0, 0, 0, 0},
 {\[CapitalGamma], 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, -\[CapitalXi], 0, 0, 0, 0},
 {0, 0, \[CapitalXi], 0, 0, 0, 0, 0}
}];
L1=ArrayFlatten[{
 {IdentityMatrix[Length[ massU]], 0, 0, 0, 0, 0, 0, 0},
 {0, IdentityMatrix[Length[massU]], 0, 0, 0, 0, 0, 0},
 {0, 0,  IdentityMatrix[Length[massQ]], 0, 0, 0, 0, 0},
 {1/\[Omega] ConjugateTranspose[S].inverseMassU, 0, 0,  IdentityMatrix[Length[massQ]], 0, 0, 0, 0},
 {0, 0, 0, 0, IdentityMatrix[Length[\[CapitalGamma]]], 0, 0, 0},
 {-(1/\[Omega])\[CapitalGamma].inverseMassU, 0, 0, 0, 0, IdentityMatrix[Length[\[CapitalGamma]]], 0, 0},
 {0, 0, 0, 0, 0, 0, IdentityMatrix[Length[\[CapitalXi]]], 0},
 {0, 0, 0, 0, 0, 0, 0, IdentityMatrix[Length[\[CapitalXi]]]}
}];
A1=ArrayFlatten[{
 {IdentityMatrix[Length[ massU]], 0, 0, 0, 0, 0, 0, 0},
 {0, IdentityMatrix[Length[massU]], 0, 0, 0, 0, 0, 0},
 {0, 0,  IdentityMatrix[Length[massQ]], 0, 0, 0, 0, 0},
 {1/\[Omega] ConjugateTranspose[S].inverseMassU, 0, 0,  IdentityMatrix[Length[massQ]], 0, 0, 0, 0},
 {0, 0, 0, 0, IdentityMatrix[Length[\[CapitalGamma]]], 0, 0, 0},
 {-(1/\[Omega])\[CapitalGamma].inverseMassU, 0, 0, 0, 0, IdentityMatrix[Length[\[CapitalGamma]]], 0, 0},
 {0, 0, 0, 0, 0, 0, IdentityMatrix[Length[\[CapitalXi]]], 0},
 {0, 0, 0, 0, 0, 0, 0, IdentityMatrix[Length[\[CapitalXi]]]}
}].SparseArray@ArrayFlatten[{
 {\[Omega] massU, 0, 0, -S, 0, -ConjugateTranspose[\[CapitalGamma]], 0, 0},
 {0, \[Omega] massU, S, 0, ConjugateTranspose[\[CapitalGamma]], 0, 0, 0},
 {0, ConjugateTranspose[S], \[Omega] massQ, -massQ, 0, 0, 0, -ConjugateTranspose[\[CapitalXi]]},
 {-ConjugateTranspose[S], 0, massQ, \[Omega] massQ, 0, 0, ConjugateTranspose[\[CapitalXi]], 0},
 {0, -\[CapitalGamma], 0, 0, 0, 0, 0, 0},
 {\[CapitalGamma], 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, -\[CapitalXi], 0, 0, 0, 0},
 {0, 0, \[CapitalXi], 0, 0, 0, 0, 0}
}]
A11=SparseArray@ArrayFlatten[{
 {\[Omega] massU, 0, 0, \[Placeholder], 0, 0, 0, 0},
 {\[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder]},
 {\[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder]},
 {\[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder]},
 {\[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder]},
 {\[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder]},
 {\[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder]},
 {\[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder]}
}]
L2=ArrayFlatten[{
 {IdentityMatrix[Length[ massU]], 0, 0, 0, 0, 0, 0, 0},
 {0, IdentityMatrix[Length[massU]], 0, 0, 0, 0, 0, 0},
 {0, 0,  IdentityMatrix[Length[massQ]], 0, 0, 0, 0, 0},
 {0, 0, 0,  IdentityMatrix[Length[massQ]], 0, 0, 0, 0},
 {0, 0, 0, 0, IdentityMatrix[Length[\[CapitalGamma]]], 0, 0, 0},
 {0, 0, 0, 0, 0, IdentityMatrix[Length[\[CapitalGamma]]], 0, 0},
 {0, 0, 0, 0, 0, 0, IdentityMatrix[Length[\[CapitalXi]]], 0},
 {0, 0, 0, 0, 0, 0, 0, IdentityMatrix[Length[\[CapitalXi]]]}
}];



Inverse[BBB]


L1.AA//MatrixPlot


inverseMassU.ConjugateTranspose[S]


?ArrayFlatten


(inverseMassQ.(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ))


inverseMassQ.(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ)


\[Kappa]//MatrixPlot
inv\[Kappa]=Inverse[\[Kappa]];=


inv\[Kappa]//MatrixPlot


\[Kappa]=(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ).inverseMassQ.(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ)+massQ;


Eigenvalues[(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ).inverseMassQ.(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ)+massQ]


(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ)


(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ)


(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ)


(1/\[Omega] ConjugateTranspose[S].inverseMassU.S-\[Omega] massQ)


\[ScriptCapitalI]u
\[ScriptCapitalI]q


Inverse[\[DoubleStruckCapitalA]]


ArrayFlatten[{{\[DoubleStruckCapitalA],\[DoubleStruckCapitalB]},{-ConjugateTranspose[\[DoubleStruckCapitalB]],0}}]





AA=SparseArray@ArrayFlatten[{
 {\[Omega] massU, 0, 0, -S, 0, -ConjugateTranspose[\[CapitalGamma]], 0, 0},
 {0, \[Omega] massU, S, 0, ConjugateTranspose[\[CapitalGamma]], 0, 0, 0},
 {0, ConjugateTranspose[S], \[Omega] massQ, -massQ, 0, 0, 0, -ConjugateTranspose[\[CapitalXi]]},
 {-ConjugateTranspose[S], 0, massQ, \[Omega] massQ, 0, 0, ConjugateTranspose[\[CapitalXi]], 0},
 {0, -\[CapitalGamma], 0, 0, 0, 0, 0, 0},
 {\[CapitalGamma], 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, -\[CapitalXi], 0, 0, 0, 0},
 {0, 0, \[CapitalXi], 0, 0, 0, 0, 0}
}];
rhsBig=Join[0 rhs1,rhs1,0rhs2,rhs2,0rhs3,rhs3,0rhs4,rhs4];
LinearSolve[AA,rhsBig]-Join[uR,uI,qR,qI,\[Lambda]R,\[Lambda]I,\[Mu]R,\[Mu]I]//Chop//Norm


AA=SparseArray@ArrayFlatten[{
 {\[Omega] massU, 0, 0, -S, 0, -ConjugateTranspose[\[CapitalGamma]], 0, 0},
 {0, \[Omega] massU, S, 0, ConjugateTranspose[\[CapitalGamma]], 0, 0, 0},
 {0, ConjugateTranspose[S], \[Omega] massQ, -massQ, 0, 0, 0, -ConjugateTranspose[\[CapitalXi]]},
 {0, 0, \[Omega] massQ, -ConjugateTranspose[S].inverseMassU.S +\[Omega]^2 massQ, 0, -ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], \[Omega] ConjugateTranspose[\[CapitalXi]], 0},
 {0, -\[CapitalGamma], 0, 0, 0, 0, 0, 0},
 {\[CapitalGamma], 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, -\[CapitalXi], 0, 0, 0, 0},
 {0, 0, \[CapitalXi], 0, 0, 0, 0, 0}
}];
rhsBig=Join[0 rhs1,rhs1,0rhs2,\[Omega] rhs2,0rhs3,rhs3,0rhs4,rhs4];
LinearSolve[AA,rhsBig]-Join[uR,uI,qR,qI,\[Lambda]R,\[Lambda]I,\[Mu]R,\[Mu]I]//Chop//Norm


AA=SparseArray@ArrayFlatten[{
 {\[Omega] massU, 0, 0, -S, 0, -ConjugateTranspose[\[CapitalGamma]], 0, 0},
 {0, \[Omega] massU, S, 0, ConjugateTranspose[\[CapitalGamma]], 0, 0, 0},
 {0, 0, ConjugateTranspose[S].inverseMassU.S -\[Omega]^2 massQ, \[Omega] massQ, ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, \[Omega] ConjugateTranspose[\[CapitalXi]]},
 {0, 0, \[Omega] massQ, -ConjugateTranspose[S].inverseMassU.S +\[Omega]^2 massQ, 0, -ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], \[Omega] ConjugateTranspose[\[CapitalXi]], 0},
 {0, -\[CapitalGamma], 0, 0, 0, 0, 0, 0},
 {\[CapitalGamma], 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, -\[CapitalXi], 0, 0, 0, 0},
 {0, 0, \[CapitalXi], 0, 0, 0, 0, 0}
}];
rhsBig=Join[0 rhs1,rhs1,ConjugateTranspose[S].inverseMassU.rhs1,\[Omega] rhs2,0rhs3,rhs3,0rhs4,rhs4];
LinearSolve[AA,rhsBig]-Join[uR,uI,qR,qI,\[Lambda]R,\[Lambda]I,\[Mu]R,\[Mu]I]//Chop//Norm

AA=SparseArray@ArrayFlatten[{
 {\[Omega] massU, 0, 0, -S, 0, -ConjugateTranspose[\[CapitalGamma]], 0, 0},
 {0, \[Omega] massU, S, 0, ConjugateTranspose[\[CapitalGamma]], 0, 0, 0},
 {0, 0, ConjugateTranspose[S].inverseMassU.S -\[Omega]^2 massQ, \[Omega] massQ, ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, \[Omega] ConjugateTranspose[\[CapitalXi]]},
 {0, 0, \[Omega] massQ, -ConjugateTranspose[S].inverseMassU.S +\[Omega]^2 massQ, 0, -ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], \[Omega] ConjugateTranspose[\[CapitalXi]], 0},
 {0, -\[Omega] \[CapitalGamma], 0, 0, 0, 0, 0, 0},
 {0, 0, 0, -\[CapitalGamma].inverseMassU.S, 0, -\[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0},
 {0, 0, 0, -\[CapitalXi], 0, 0, 0, 0},
 {0, 0, \[CapitalXi], 0, 0, 0, 0, 0}
}];
rhsBig=Join[0 rhs1,rhs1,ConjugateTranspose[S].inverseMassU.rhs1,\[Omega] rhs2,0rhs3,-\[Omega] rhs3,0rhs4,rhs4];
LinearSolve[AA,rhsBig]-Join[uR,uI,qR,qI,\[Lambda]R,\[Lambda]I,\[Mu]R,\[Mu]I]//Chop//Norm



AA=SparseArray@ArrayFlatten[{
 {\[Omega] massU, 0, 0, -S, 0, -ConjugateTranspose[\[CapitalGamma]], 0, 0},
 {0, \[Omega] massU, S, 0, ConjugateTranspose[\[CapitalGamma]], 0, 0, 0},
 {0, 0, ConjugateTranspose[S].inverseMassU.S -\[Omega]^2 massQ, \[Omega] massQ, ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, \[Omega] ConjugateTranspose[\[CapitalXi]]},
 {0, 0, \[Omega] massQ, -(ConjugateTranspose[S].inverseMassU.S -\[Omega]^2 massQ), 0, -ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], \[Omega] ConjugateTranspose[\[CapitalXi]], 0},
 {0, 0, \[CapitalGamma].inverseMassU.S, 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0},
 {0, 0, 0, -\[CapitalGamma].inverseMassU.S, 0, -\[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0},
 {0, 0, 0, -\[CapitalXi], 0, 0, 0, 0},
 {0, 0, \[CapitalXi], 0, 0, 0, 0, 0}
}];
rhsBig=Join[0 rhs1,rhs1,ConjugateTranspose[S].inverseMassU.rhs1,\[Omega] rhs2,\[CapitalGamma].inverseMassU.rhs1,-\[Omega] rhs3,0rhs4,rhs4];
LinearSolve[AA,rhsBig]-Join[uR,uI,qR,qI,\[Lambda]R,\[Lambda]I,\[Mu]R,\[Mu]I]//Chop//Norm

AA=SparseArray@ArrayFlatten[{
 {\[Omega] massU, 0, 0, -S, 0, -ConjugateTranspose[\[CapitalGamma]], 0, 0},
 {0, \[Omega] massU, S, 0, ConjugateTranspose[\[CapitalGamma]], 0, 0, 0},
 {0, 0, ConjugateTranspose[S].inverseMassU.S -\[Omega]^2 massQ, \[Omega] massQ, ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, \[Omega] ConjugateTranspose[\[CapitalXi]]},
 {0, 0, \[Omega] massQ, -(ConjugateTranspose[S].inverseMassU.S -\[Omega]^2 massQ), 0, -ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], \[Omega] ConjugateTranspose[\[CapitalXi]], 0},
 {0, 0, \[CapitalGamma].inverseMassU.S, 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0},
 {0, 0, 0, -\[CapitalGamma].inverseMassU.S, 0, -\[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0},
 {0, 0, 0, -\[CapitalXi], 0, 0, 0, 0},
 {0, 0, 0, -\[CapitalXi].inverseMassQ.(ConjugateTranspose[S].inverseMassU.S -\[Omega]^2 massQ), 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], \[Omega] \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0}
}];
rhsBig=Join[0 rhs1,rhs1,ConjugateTranspose[S].inverseMassU.rhs1,\[Omega] rhs2,\[CapitalGamma].inverseMassU.rhs1,-\[Omega] rhs3,0rhs4,\[Omega] \[CapitalXi].inverseMassQ.rhs2-\[Omega] rhs4];
LinearSolve[AA,rhsBig]-Join[uR,uI,qR,qI,\[Lambda]R,\[Lambda]I,\[Mu]R,\[Mu]I]//Chop//Norm

AA=SparseArray@ArrayFlatten[{
 {\[Omega] massU, 0, 0, -S, 0, -ConjugateTranspose[\[CapitalGamma]], 0, 0},
 {0, \[Omega] massU, S, 0, ConjugateTranspose[\[CapitalGamma]], 0, 0, 0},
 {0, 0, ConjugateTranspose[S].inverseMassU.S -\[Omega]^2 massQ, \[Omega] massQ, ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, \[Omega] ConjugateTranspose[\[CapitalXi]]},
 {0, 0, \[Omega] massQ, -(ConjugateTranspose[S].inverseMassU.S -\[Omega]^2 massQ), 0, -ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], \[Omega] ConjugateTranspose[\[CapitalXi]], 0},
 {0, 0, \[CapitalGamma].inverseMassU.S, 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0},
 {0, 0, 0, -\[CapitalGamma].inverseMassU.S, 0, -\[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0},
 {0, 0, \[CapitalXi].inverseMassQ.(ConjugateTranspose[S].inverseMassU.S -\[Omega]^2 massQ), 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, \[Omega] \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]]},
 {0, 0, 0, -\[CapitalXi].inverseMassQ.(ConjugateTranspose[S].inverseMassU.S -\[Omega]^2 massQ), 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.ConjugateTranspose[\[CapitalGamma]], \[Omega] \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0}
}];
rhsBig=Join[0 rhs1,rhs1,ConjugateTranspose[S].inverseMassU.rhs1,\[Omega] rhs2,\[CapitalGamma].inverseMassU.rhs1,-\[Omega] rhs3,\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.rhs1,\[Omega] \[CapitalXi].inverseMassQ.rhs2-\[Omega] rhs4];
LinearSolve[AA,rhsBig]-Join[uR,uI,qR,qI,\[Lambda]R,\[Lambda]I,\[Mu]R,\[Mu]I]//Chop//Norm




AA=SparseArray@ArrayFlatten[{
 {\[Omega] massU, 0, 0, -S, 0, -ConjugateTranspose[\[CapitalGamma]], 0, 0},
 {0, \[Omega] massU, S, 0, ConjugateTranspose[\[CapitalGamma]], 0, 0, 0},
 {0, ConjugateTranspose[S], \[Omega] massQ, -massQ, 0, 0, 0, -ConjugateTranspose[\[CapitalXi]]},
 {-ConjugateTranspose[S], 0, massQ, \[Omega] massQ, 0, 0, ConjugateTranspose[\[CapitalXi]], 0},
 {0, -\[CapitalGamma], 0, 0, 0, 0, 0, 0},
 {0, 0, 0, -\[CapitalGamma].inverseMassU.S, 0, -\[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0},
 {0, 0, 0, -\[CapitalXi], 0, 0, 0, 0},
 {0, 0, \[CapitalXi], 0, 0, 0, 0, 0}
}];
rhsBig=Join[0 rhs1,rhs1,0rhs2,rhs2,0rhs3,-\[Omega] rhs3,0rhs4,rhs4];
LinearSolve[AA,rhsBig]-Join[uR,uI,qR,qI,\[Lambda]R,\[Lambda]I,\[Mu]R,\[Mu]I]//Chop//Norm

AA=SparseArray@ArrayFlatten[{
 {\[Omega] massU, 0, 0, -S, 0, -ConjugateTranspose[\[CapitalGamma]], 0, 0},
 {0, \[Omega] massU, S, 0, ConjugateTranspose[\[CapitalGamma]], 0, 0, 0},
 {0, ConjugateTranspose[S], \[Omega] massQ, -massQ, 0, 0, 0, -ConjugateTranspose[\[CapitalXi]]},
 {-ConjugateTranspose[S], 0, massQ, \[Omega] massQ, 0, 0, ConjugateTranspose[\[CapitalXi]], 0},
 {0, -\[CapitalGamma], \[CapitalGamma].inverseMassU.S, 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0},
 {0, 0, 0, -\[CapitalGamma].inverseMassU.S, 0, -\[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0},
 {0, 0, 0, -\[CapitalXi], 0, 0, 0, 0},
 {0, 0, \[CapitalXi], 0, 0, 0, 0, 0}
}];
rhsBig=Join[0 rhs1,rhs1,0rhs2,rhs2,\[CapitalGamma].inverseMassU.rhs1+0\[Omega] rhs3,-\[Omega] rhs3,0rhs4,rhs4];
LinearSolve[AA,rhsBig]-Join[uR,uI,qR,qI,\[Lambda]R,\[Lambda]I,\[Mu]R,\[Mu]I]//Chop//Norm

AA=SparseArray@ArrayFlatten[{
 {\[Omega] massU, 0, 0, -S, 0, -ConjugateTranspose[\[CapitalGamma]], 0, 0},
 {0, \[Omega] massU, S, 0, ConjugateTranspose[\[CapitalGamma]], 0, 0, 0},
 {0, ConjugateTranspose[S], \[Omega] massQ, -massQ, 0, 0, 0, -ConjugateTranspose[\[CapitalXi]]},
 {-ConjugateTranspose[S], 0, massQ, \[Omega] massQ, 0, 0, ConjugateTranspose[\[CapitalXi]], 0},
 {0, 0, \[CapitalGamma].inverseMassU.S, 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0},
 {0, 0, 0, -\[CapitalGamma].inverseMassU.S, 0, -\[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0},
 {0, 0, 0, -\[CapitalXi], 0, 0, 0, 0},
 {0, 0, \[CapitalXi], 0, 0, 0, 0, 0}
}];
rhsBig=Join[0 rhs1,rhs1,0rhs2,rhs2,\[CapitalGamma].inverseMassU.rhs1+0\[Omega] rhs3,-\[Omega] rhs3,0rhs4,rhs4];
LinearSolve[AA,rhsBig]-Join[uR,uI,qR,qI,\[Lambda]R,\[Lambda]I,\[Mu]R,\[Mu]I]//Chop//Norm


AA=SparseArray@ArrayFlatten[{
 {\[Omega] massU, 0, 0, -S, 0, -ConjugateTranspose[\[CapitalGamma]], 0, 0},
 {0, \[Omega] massU, S, 0, ConjugateTranspose[\[CapitalGamma]], 0, 0, 0},
 {0, ConjugateTranspose[S], \[Omega] massQ, -massQ, 0, 0, 0, -ConjugateTranspose[\[CapitalXi]]},
 {-ConjugateTranspose[S], 0, massQ, \[Omega] massQ, 0, 0, ConjugateTranspose[\[CapitalXi]], 0},
 {0, 0, \[CapitalGamma].inverseMassU.S, 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0},
 {0, 0, 0, -\[CapitalGamma].inverseMassU.S, 0, -\[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0},
 {0, 0, 0, -\[CapitalXi], 0, 0, 0, 0},
 {-\[CapitalXi].inverseMassQ.ConjugateTranspose[S], 0, 0, \[Omega] \[CapitalXi].inverseMassQ.massQ, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0}
}];

rhsBig=Join[0 rhs1,rhs1,0rhs2,rhs2,\[CapitalGamma].inverseMassU.rhs1+0\[Omega] rhs3,-\[Omega] rhs3,0rhs4,\[CapitalXi].inverseMassQ.rhs2-rhs4];
LinearSolve[AA,rhsBig]-Join[uR,uI,qR,qI,\[Lambda]R,\[Lambda]I,\[Mu]R,\[Mu]I]//Chop//Norm

AA=SparseArray@ArrayFlatten[{
 {\[Omega] massU, 0, 0, -S, 0, -ConjugateTranspose[\[CapitalGamma]], 0, 0},
 {0, \[Omega] massU, S, 0, ConjugateTranspose[\[CapitalGamma]], 0, 0, 0},
 {0, ConjugateTranspose[S], \[Omega] massQ, -massQ, 0, 0, 0, -ConjugateTranspose[\[CapitalXi]]},
 {-ConjugateTranspose[S], 0, massQ, \[Omega] massQ, 0, 0, ConjugateTranspose[\[CapitalXi]], 0},
 {0, 0, \[CapitalGamma].inverseMassU.S, 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0},
 {0, 0, 0, -\[CapitalGamma].inverseMassU.S, 0, -\[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0},
 {0, \[CapitalXi].inverseMassQ.ConjugateTranspose[S], \[Omega] \[CapitalXi].inverseMassQ.massQ, 0, 0, 0, 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]]},
 {-\[CapitalXi].inverseMassQ.ConjugateTranspose[S], 0, 0, \[Omega] \[CapitalXi].inverseMassQ.massQ, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0}
}];

rhsBig=Join[0 rhs1,rhs1,0rhs2,rhs2,\[CapitalGamma].inverseMassU.rhs1+0\[Omega] rhs3,-\[Omega] rhs3,0 rhs4,\[CapitalXi].inverseMassQ.rhs2-rhs4];
LinearSolve[AA,rhsBig]-Join[uR,uI,qR,qI,\[Lambda]R,\[Lambda]I,\[Mu]R,\[Mu]I]//Chop//Norm

AA=SparseArray@ArrayFlatten[{
 {\[Omega] massU, 0, 0, -S, 0, -ConjugateTranspose[\[CapitalGamma]], 0, 0},
 {0, \[Omega] massU, S, 0, ConjugateTranspose[\[CapitalGamma]], 0, 0, 0},
 {0, ConjugateTranspose[S], \[Omega] massQ, -massQ, 0, 0, 0, -ConjugateTranspose[\[CapitalXi]]},
 {-ConjugateTranspose[S], 0, massQ, \[Omega] massQ, 0, 0, ConjugateTranspose[\[CapitalXi]], 0},
 {0, 0, \[CapitalGamma].inverseMassU.S, 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0},
 {0, 0, 0, -\[CapitalGamma].inverseMassU.S, 0, -\[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0},
 {0, \[CapitalXi].inverseMassQ.ConjugateTranspose[S], \[Omega] \[CapitalXi].inverseMassQ.massQ, 0, 0, 0, 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]]},
 {-\[CapitalXi].inverseMassQ.ConjugateTranspose[S], 0, 0, \[Omega] \[CapitalXi].inverseMassQ.massQ, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0}
}];

rhsBig=Join[0 rhs1,rhs1,0rhs2,rhs2,\[CapitalGamma].inverseMassU.rhs1+0\[Omega] rhs3,-\[Omega] rhs3,0 rhs4,\[CapitalXi].inverseMassQ.rhs2-rhs4];
LinearSolve[AA,rhsBig]-Join[uR,uI,qR,qI,\[Lambda]R,\[Lambda]I,\[Mu]R,\[Mu]I]//Chop//Norm


AA=SparseArray@ArrayFlatten[{
 {\[Omega] massU, 0, 0, -S, 0, -ConjugateTranspose[\[CapitalGamma]], 0, 0},
 {0, \[Omega] massU, S, 0, ConjugateTranspose[\[CapitalGamma]], 0, 0, 0},
 {0, ConjugateTranspose[S], \[Omega] massQ, -massQ, 0, 0, 0, -ConjugateTranspose[\[CapitalXi]]},
 {-ConjugateTranspose[S], 0, massQ, \[Omega] massQ, 0, 0, ConjugateTranspose[\[CapitalXi]], 0},
 {0, 0, \[CapitalGamma].inverseMassU.S, 0, \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0, 0},
 {0, 0, 0, -\[CapitalGamma].inverseMassU.S, 0, -\[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]], 0, 0},
 {0, \[CapitalXi].inverseMassQ.ConjugateTranspose[S], \[Omega] \[CapitalXi].inverseMassQ.massQ, 0, 0, 0, 0, -\[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]]},
 {-\[CapitalXi].inverseMassQ.ConjugateTranspose[S], 0, 0, \[Omega] \[CapitalXi].inverseMassQ.massQ, 0, 0, \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]], 0}
}];

rhsBig=Join[0 rhs1,rhs1,0rhs2,rhs2,\[CapitalGamma].inverseMassU.rhs1+0\[Omega] rhs3,-\[Omega] rhs3,0 rhs4,\[CapitalXi].inverseMassQ.rhs2-rhs4];
LinearSolve[AA,rhsBig]-Join[uR,uI,qR,qI,\[Lambda]R,\[Lambda]I,\[Mu]R,\[Mu]I]//Chop//Norm


\[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]]//MatrixPlot
\[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]]//MatrixPlot


AA//MatrixPlot


AA//MatrixPlot


\[Omega] \[CapitalXi].inverseMassQ.massQ.qI//Chop
\[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]].\[Mu]R
aa \[CapitalXi].inverseMassQ.rhs2-bb rhs4


massQ.qR+\[Omega] massQ.qI+ConjugateTranspose[\[CapitalXi]].\[Mu]R





\[CapitalXi].inverseMassQ.rhs2


AA//MatrixPlot


\[CapitalGamma].inverseMassU.rhs1


AA.Join[Join[uR,uI,qR,qI,\[Lambda]R,\[Lambda]I,\[Mu]R,\[Mu]I]]-rhsBig//Chop


MatrixPlot[AA]





\[CapitalXi].inverseMassQ.ConjugateTranspose[S].inverseMassU.rhs1


AA//MatrixPlot


\[CapitalGamma].inverseMassU.rhs1


\[ScriptCapitalI]u=IdentityMatrix[Length[massU]];
\[ScriptCapitalI]q=IdentityMatrix[Length[massQ]];








Join[0 rhs1 ]











128


rhs4//Length


uI


rhs1


r





\[CapitalXi]R
\[CapitalXi]I


\[Lambda]Soln


ConjugateTranspose[\[CapitalGamma]]//MatrixPlot


Length[ConjugateTranspose[\[CapitalGamma]]]


AA=ArrayFlatten[{
 {0.0massU, 0, 0, -1.0S, {
   {0, -ConjugateTranspose[\[CapitalGamma]], 0, 0},
   {ConjugateTranspose[\[CapitalGamma]], 0, 0, 0},
   {0, 0, 0, -ConjugateTranspose[\[CapitalXi]]},
   {0, 0, ConjugateTranspose[\[CapitalXi]], 0}
  }, \[Placeholder], \[Placeholder], \[Placeholder]},
 {0, 0.0massU, S, 0, \[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder]},
 {0, ConjugateTranspose[S], 0.0massQ, -massQ, \[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder]},
 {-ConjugateTranspose[S], 0, massQ, 0.0massQ, \[Placeholder], \[Placeholder], \[Placeholder], \[Placeholder]},
 {0, -\[CapitalGamma], 0, 0, 0, 0, 0, 0},
 {\[CapitalGamma], 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, \[CapitalXi], 0, 0, 0, 0},
 {0, 0, -\[CapitalXi], 0, 0, 0, 0, 0}
}];

BB=ArrayFlatten[{
 {massU, 0, 0, -0.0S, 0, 0.0ConjugateTranspose[\[CapitalGamma]], 0, 0},
 {0, massU, 0.0S, 0, 0.0ConjugateTranspose[\[CapitalGamma]], 0, 0, 0},
 {0, 0.0ConjugateTranspose[S], massQ, -0.0massQ, 0, 0, 0, -0.0 ConjugateTranspose[\[CapitalXi]]},
 {-0.0 ConjugateTranspose[S], 0, 0.0massQ, massQ, 0, 0, 0.0ConjugateTranspose[\[CapitalXi]], 0},
 {0, 0.0\[CapitalGamma], 0, 0, 0, 0, 0, 0},
 {0.0\[CapitalGamma], 0, 0, 0, 0, 0, 0, 0},
 {0, 0, 0, 0.0\[CapitalXi], 0, 0, 0, 0},
 {0, 0, 0.0\[CapitalXi], 0, 0, 0, 0, 0}
}];


\[DoubleStruckCapitalA]=SparseArray@ArrayFlatten[{
 {massU, 0, 0, -0.0S},
 {0, massU, 0.0S, 0},
 {0, 0.0ConjugateTranspose[S], massQ, -0.0massQ},
 {-0.0 ConjugateTranspose[S], 0, 0.0massQ, massQ}
}];
\[DoubleStruckCapitalB]=SparseArray@ArrayFlatten[{
 {0.0massU, 0, 0, -1.0S},
 {0, 0.0massU, 1.0S, 0},
 {0, 1.0ConjugateTranspose[S], 0.0massQ, -1.0massQ},
 {-1.0 ConjugateTranspose[S], 0, 1.0massQ, 0.0massQ}
}];
\[DoubleStruckCapitalC]=SparseArray@ArrayFlatten[{
 {0, -ConjugateTranspose[\[CapitalGamma]], 0, 0},
 {ConjugateTranspose[\[CapitalGamma]], 0, 0, 0},
 {0, 0, 0, -ConjugateTranspose[\[CapitalXi]]},
 {0, 0, ConjugateTranspose[\[CapitalXi]], 0}
}];




\[DoubleStruckCapitalD]=\[DoubleStruckCapitalA]+\[DoubleStruckCapitalB];


ConjugateTranspose[\[DoubleStruckCapitalC]].LinearSolve[\[DoubleStruckCapitalD],\[DoubleStruckCapitalC]]//MatrixPlot





\[DoubleStruckCapitalA]//MatrixP


//


\[DoubleStruckCapitalA]//MatrixPlot
\[DoubleStruckCapitalB]


AA//MatrixPlot


AA


AA//MatrixPlot


af={
 {1, 2},
 {3, 4}
}


AA//MatrixPlot
BB//MatrixPlot
CC//MatrixPlot


ZZ=ArrayFlatten[{{BB,CC},{ConjugateTranspose[CC],0}}];
ArrayFlatten[{{1.0AA,0.0CC},{0.0ConjugateTranspose[CC],0.0}}]


ZZ//SingularValueList






zzz=SparseArray@ArrayFlatten[{{massQ,S,ConjugateTranspose[\[CapitalXi]],0},{-ConjugateTranspose[R],0,0,ConjugateTranspose[\[CapitalGamma]]},{\[CapitalXi],0,0,0},{0,\[CapitalGamma],0,0}}];
MatrixPlot@zzz


Q=RandomInteger[{-10,10},Length[massQ]];
U=RandomInteger[{-10,10},Length[ConjugateTranspose[R]]];
\[Mu]=RandomInteger[{-10,10},Length[\[CapitalXi]]];
\[Lambda]=RandomInteger[{-10,10},Length[\[CapitalGamma]]];

f=massQ.Q+S.U+ConjugateTranspose[\[CapitalXi]].\[Mu];
g=-ConjugateTranspose[R].Q+ ConjugateTranspose[\[CapitalGamma]].\[Lambda];
h=\[CapitalXi].Q;
j=\[CapitalGamma].U;




\[DoubleStruckCapitalY]=SparseArray@Chop@ConjugateTranspose[NullSpace[\[CapitalXi]]];
\[DoubleStruckCapitalZ]=ConjugateTranspose[(SparseArray@Chop@NullSpace[\[CapitalGamma]])];
a=LinearSolve[ConjugateTranspose[\[DoubleStruckCapitalY]].massQ.\[DoubleStruckCapitalY],ConjugateTranspose[\[DoubleStruckCapitalY]].massQ.Q];

\[ScriptL]=LinearSolve[\[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]],\[CapitalXi].Q];

b=LinearSolve[ConjugateTranspose[\[DoubleStruckCapitalZ]].massU.\[DoubleStruckCapitalZ],ConjugateTranspose[\[DoubleStruckCapitalZ]].massU.U];
\[ScriptM]=LinearSolve[\[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]],\[CapitalGamma].U];


\[DoubleStruckCapitalY][[1;;4,105;;108]]//MatrixPlot


inverseMassU.ConjugateTranspose[\[CapitalGamma]].\[ScriptM]+\[DoubleStruckCapitalZ].b-U//Chop



ConjugateTranspose[\[DoubleStruckCapitalY]].(massQ.Q+S.U+ConjugateTranspose[\[CapitalXi]].\[Mu]-f)//Chop
ConjugateTranspose[\[DoubleStruckCapitalY]].(massQ.(inverseMassQ.ConjugateTranspose[\[CapitalXi]].\[ScriptL]+\[DoubleStruckCapitalY].a)+S.U+ConjugateTranspose[\[CapitalXi]].\[Mu]-f)//Chop


hh=SparseArray@LinearSolve[ConjugateTranspose[\[DoubleStruckCapitalZ]].massU.\[DoubleStruckCapitalZ],SparseArray@(ConjugateTranspose[\[DoubleStruckCapitalZ]].ConjugateTranspose[R].\[DoubleStruckCapitalY].SparseArray@LinearSolve[ConjugateTranspose[\[DoubleStruckCapitalY]].massQ.\[DoubleStruckCapitalY],ConjugateTranspose[\[DoubleStruckCapitalY]].S.\[DoubleStruckCapitalZ]])];


SparseArray@LinearSolve[ConjugateTranspose[\[DoubleStruckCapitalZ]].massU.\[DoubleStruckCapitalZ],SparseArray[ConjugateTranspose[\[DoubleStruckCapitalZ]].ConjugateTranspose[S].\[DoubleStruckCapitalY].LinearSolve[ConjugateTranspose[\[DoubleStruckCapitalY]].massQ.\[DoubleStruckCapitalY],ConjugateTranspose[\[DoubleStruckCapitalY]].S.\[DoubleStruckCapitalZ]]]]//Chop//MatrixPlot


\[DoubleStruckCapitalZ]//MatrixPlot
massU//MatrixPlot


Inverse[ConjugateTranspose[\[DoubleStruckCapitalZ]].massU.\[DoubleStruckCapitalZ]]-ConjugateTranspose[\[DoubleStruckCapitalZ]].inverseMassU.\[DoubleStruckCapitalZ]


Inverse[ConjugateTranspose[\[DoubleStruckCapitalZ]].massU.\[DoubleStruckCapitalZ]]//Norm


massU//MatrixPlot


64-48


(1/#& /@ Eigenvalues[{ConjugateTranspose[\[DoubleStruckCapitalZ]].massU.\[DoubleStruckCapitalZ],SparseArray@(ConjugateTranspose[\[DoubleStruckCapitalZ]].ConjugateTranspose[R].\[DoubleStruckCapitalY].SparseArray@LinearSolve[ConjugateTranspose[\[DoubleStruckCapitalY]].massQ.\[DoubleStruckCapitalY],ConjugateTranspose[\[DoubleStruckCapitalY]].S.\[DoubleStruckCapitalZ]])}])//ListPlot


hh//Chop


Eigenvalues[hh]


ConjugateTranspose[\[DoubleStruckCapitalY]].S.\[DoubleStruckCapitalZ]-ConjugateTranspose[\[DoubleStruckCapitalY]].R.\[DoubleStruckCapitalZ]//Chop//ArrayRules//Short


S+R//Chop


S


{omegaSquared,b}=(Eigensystem[#,-1]&@SparseArray@LinearSolve[-ConjugateTranspose[\[DoubleStruckCapitalZ]].massU.\[DoubleStruckCapitalZ],ConjugateTranspose[\[DoubleStruckCapitalZ]].ConjugateTranspose[R].\[DoubleStruckCapitalY].SparseArray@LinearSolve[ConjugateTranspose[\[DoubleStruckCapitalY]].massQ.\[DoubleStruckCapitalY],ConjugateTranspose[\[DoubleStruckCapitalY]].S.\[DoubleStruckCapitalZ]]]);
b=First[b];
omegaSquared=First[omegaSquared];
omega=omegaSquared^(1/2);


1/omega LinearSolve[ConjugateTranspose[\[DoubleStruckCapitalY]].massQ.\[DoubleStruckCapitalY],ConjugateTranspose[\[DoubleStruckCapitalY]].S.\[DoubleStruckCapitalZ].b]//Chop





\[GothicCapitalQ]=SparseArray@(ConjugateTranspose[\[CapitalGamma]].LinearSolve[\[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]],\[CapitalGamma].inverseMassU]//Chop)
\[GothicCapitalP]=SparseArray@(ConjugateTranspose[\[CapitalXi]].LinearSolve[\[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]],\[CapitalXi].inverseMassQ]//Chop)


\[ScriptCapitalI]2=SparseArray@IdentityMatrix[Length[Q]];
\[ScriptCapitalI]1=SparseArray@IdentityMatrix[Length[U]];


inverseMassU.(\[ScriptCapitalI]1-\[GothicCapitalQ]).ConjugateTranspose[R].inverseMassQ.(\[ScriptCapitalI]2-\[GothicCapitalP]).S


hh


244-208


Eigenvalues[inverseMassU.(\[ScriptCapitalI]1-\[GothicCapitalQ]).ConjugateTranspose[R].inverseMassQ.(\[ScriptCapitalI]2-\[GothicCapitalP]).S]//Chop
Eigenvalues[hh]


Chop[hh]
MatrixPlot[Chop[hh]]


Eigenvalues[ConjugateTranspose[\[DoubleStruckCapitalZ]].ConjugateTranspose[R].ConjugateTranspose[\[DoubleStruckCapitalZ]].massU.\[DoubleStruckCapitalZ]


(*\[DoubleStruckCapitalZ]^\[ConjugateTranspose].massQ.\[DoubleStruckCapitalZ].a+\[DoubleStruckCapitalZ]^\[ConjugateTranspose].massQ.inverseMassQ.\[CapitalXi]^\[ConjugateTranspose].\[ScriptL]+ \[DoubleStruckCapitalZ]^\[ConjugateTranspose].S.\[DoubleStruckCapitalY].b+\[DoubleStruckCapitalZ]^\[ConjugateTranspose].S.inverseMassU.\[CapitalGamma]^\[ConjugateTranspose].\[ScriptM]+\[DoubleStruckCapitalZ]^\[ConjugateTranspose].\[CapitalXi]^\[ConjugateTranspose].\[Mu]-\[DoubleStruckCapitalZ]^\[ConjugateTranspose].f//Chop//Norm
\[CapitalXi].inverseMassQ.\[CapitalXi]^\[ConjugateTranspose].\[ScriptL]+ \[CapitalXi].inverseMassQ.S.\[DoubleStruckCapitalY].b+\[CapitalXi].inverseMassQ.S.inverseMassU.\[CapitalGamma]^\[ConjugateTranspose].\[ScriptM]+ \[CapitalXi].inverseMassQ.\[CapitalXi]^\[ConjugateTranspose].\[Mu]-\[CapitalXi].inverseMassQ.f//Chop//Norm

-\[DoubleStruckCapitalY]^\[ConjugateTranspose].R^\[ConjugateTranspose].\[DoubleStruckCapitalZ].a-\[DoubleStruckCapitalY]^\[ConjugateTranspose].R^\[ConjugateTranspose].inverseMassQ.\[CapitalXi]^\[ConjugateTranspose].\[ScriptL]+ \[DoubleStruckCapitalY]^\[ConjugateTranspose].\[CapitalGamma]^\[ConjugateTranspose].\[Lambda]-\[DoubleStruckCapitalY]^\[ConjugateTranspose].g//Chop//Norm
-\[CapitalGamma].inverseMassU.R^\[ConjugateTranspose].\[DoubleStruckCapitalZ].a - \[CapitalGamma].inverseMassU.R^\[ConjugateTranspose].inverseMassQ.\[CapitalXi]^\[ConjugateTranspose].\[ScriptL]+ \[CapitalGamma].inverseMassU.\[CapitalGamma]^\[ConjugateTranspose].\[Lambda]-\[CapitalGamma].inverseMassU.g//Chop//Norm

\[CapitalXi].\[DoubleStruckCapitalZ].a+\[CapitalXi].inverseMassQ.\[CapitalXi]^\[ConjugateTranspose].\[ScriptL]-h//Chop//Norm

\[CapitalGamma].\[DoubleStruckCapitalY].b+\[CapitalGamma].inverseMassU.\[CapitalGamma]^\[ConjugateTranspose].\[ScriptM]-j//Chop//Norm*)


(*\[DoubleStruckCapitalZ]^\[ConjugateTranspose].massQ.\[DoubleStruckCapitalZ].a+0\[DoubleStruckCapitalZ]^\[ConjugateTranspose].massQ.inverseMassQ.\[CapitalXi]^\[ConjugateTranspose].\[ScriptL]+ \[DoubleStruckCapitalZ]^\[ConjugateTranspose].S.\[DoubleStruckCapitalY].b+\[DoubleStruckCapitalZ]^\[ConjugateTranspose].S.inverseMassU.\[CapitalGamma]^\[ConjugateTranspose].\[ScriptM]+0\[DoubleStruckCapitalZ]^\[ConjugateTranspose].\[CapitalXi]^\[ConjugateTranspose].\[Mu]-\[DoubleStruckCapitalZ]^\[ConjugateTranspose].f//Chop//Norm
\[CapitalXi].inverseMassQ.\[CapitalXi]^\[ConjugateTranspose].\[ScriptL]+ \[CapitalXi].inverseMassQ.S.\[DoubleStruckCapitalY].b+\[CapitalXi].inverseMassQ.S.inverseMassU.\[CapitalGamma]^\[ConjugateTranspose].\[ScriptM]+ \[CapitalXi].inverseMassQ.\[CapitalXi]^\[ConjugateTranspose].\[Mu]-\[CapitalXi].inverseMassQ.f//Chop//Norm

-\[DoubleStruckCapitalY]^\[ConjugateTranspose].R^\[ConjugateTranspose].\[DoubleStruckCapitalZ].a-\[DoubleStruckCapitalY]^\[ConjugateTranspose].R^\[ConjugateTranspose].inverseMassQ.\[CapitalXi]^\[ConjugateTranspose].\[ScriptL]+ 0\[DoubleStruckCapitalY]^\[ConjugateTranspose].\[CapitalGamma]^\[ConjugateTranspose].\[Lambda]-\[DoubleStruckCapitalY]^\[ConjugateTranspose].g//Chop//Norm
-\[CapitalGamma].inverseMassU.R^\[ConjugateTranspose].\[DoubleStruckCapitalZ].a - \[CapitalGamma].inverseMassU.R^\[ConjugateTranspose].inverseMassQ.\[CapitalXi]^\[ConjugateTranspose].\[ScriptL]+ \[CapitalGamma].inverseMassU.\[CapitalGamma]^\[ConjugateTranspose].\[Lambda]-\[CapitalGamma].inverseMassU.g//Chop//Norm

0\[CapitalXi].\[DoubleStruckCapitalZ].a+\[CapitalXi].inverseMassQ.\[CapitalXi]^\[ConjugateTranspose].\[ScriptL]-h//Chop//Norm

0\[CapitalGamma].\[DoubleStruckCapitalY].b+\[CapitalGamma].inverseMassU.\[CapitalGamma]^\[ConjugateTranspose].\[ScriptM]-j//Chop//Norm*)


(*\[DoubleStruckCapitalZ]^\[ConjugateTranspose].massQ.\[DoubleStruckCapitalZ].a+ \[DoubleStruckCapitalZ]^\[ConjugateTranspose].S.\[DoubleStruckCapitalY].b+\[DoubleStruckCapitalZ]^\[ConjugateTranspose].S.inverseMassU.\[CapitalGamma]^\[ConjugateTranspose].\[ScriptM]-\[DoubleStruckCapitalZ]^\[ConjugateTranspose].f//Chop//Norm
\[CapitalXi].inverseMassQ.\[CapitalXi]^\[ConjugateTranspose].\[ScriptL]+ \[CapitalXi].inverseMassQ.S.\[DoubleStruckCapitalY].b+\[CapitalXi].inverseMassQ.S.inverseMassU.\[CapitalGamma]^\[ConjugateTranspose].\[ScriptM]+ \[CapitalXi].inverseMassQ.\[CapitalXi]^\[ConjugateTranspose].\[Mu]-\[CapitalXi].inverseMassQ.f//Chop//Norm

-\[DoubleStruckCapitalY]^\[ConjugateTranspose].R^\[ConjugateTranspose].\[DoubleStruckCapitalZ].a-\[DoubleStruckCapitalY]^\[ConjugateTranspose].R^\[ConjugateTranspose].inverseMassQ.\[CapitalXi]^\[ConjugateTranspose].\[ScriptL]-\[DoubleStruckCapitalY]^\[ConjugateTranspose].g//Chop//Norm
-\[CapitalGamma].inverseMassU.R^\[ConjugateTranspose].\[DoubleStruckCapitalZ].a - \[CapitalGamma].inverseMassU.R^\[ConjugateTranspose].inverseMassQ.\[CapitalXi]^\[ConjugateTranspose].\[ScriptL]+ \[CapitalGamma].inverseMassU.\[CapitalGamma]^\[ConjugateTranspose].\[Lambda]-\[CapitalGamma].inverseMassU.g//Chop//Norm

\[CapitalXi].inverseMassQ.\[CapitalXi]^\[ConjugateTranspose].\[ScriptL]-h//Chop//Norm

\[CapitalGamma].inverseMassU.\[CapitalGamma]^\[ConjugateTranspose].\[ScriptM]-j//Chop//Norm\[DoubleStruckCapitalZ]^\[ConjugateTranspose].massQ.\[DoubleStruckCapitalZ].a+0\[DoubleStruckCapitalZ]^\[ConjugateTranspose].massQ.inverseMassQ.\[CapitalXi]^\[ConjugateTranspose].\[ScriptL]+ \[DoubleStruckCapitalZ]^\[ConjugateTranspose].S.\[DoubleStruckCapitalY].b+\[DoubleStruckCapitalZ]^\[ConjugateTranspose].S.inverseMassU.\[CapitalGamma]^\[ConjugateTranspose].\[ScriptM]+0\[DoubleStruckCapitalZ]^\[ConjugateTranspose].\[CapitalXi]^\[ConjugateTranspose].\[Mu]-\[DoubleStruckCapitalZ]^\[ConjugateTranspose].f//Chop//Norm*)

(*\[DoubleStruckCapitalZ]^\[ConjugateTranspose].massQ.\[DoubleStruckCapitalZ].a+ \[DoubleStruckCapitalZ]^\[ConjugateTranspose].S.\[DoubleStruckCapitalY].b+\[DoubleStruckCapitalZ]^\[ConjugateTranspose].S.inverseMassU.\[CapitalGamma]^\[ConjugateTranspose].\[ScriptM]-\[DoubleStruckCapitalZ]^\[ConjugateTranspose].f//Chop//Norm
\[CapitalXi].inverseMassQ.\[CapitalXi]^\[ConjugateTranspose].\[ScriptL]+ \[CapitalXi].inverseMassQ.S.\[DoubleStruckCapitalY].b+\[CapitalXi].inverseMassQ.S.inverseMassU.\[CapitalGamma]^\[ConjugateTranspose].\[ScriptM]+ \[CapitalXi].inverseMassQ.\[CapitalXi]^\[ConjugateTranspose].\[Mu]-\[CapitalXi].inverseMassQ.f//Chop//Norm

-\[DoubleStruckCapitalY]^\[ConjugateTranspose].R^\[ConjugateTranspose].\[DoubleStruckCapitalZ].a-\[DoubleStruckCapitalY]^\[ConjugateTranspose].R^\[ConjugateTranspose].inverseMassQ.\[CapitalXi]^\[ConjugateTranspose].\[ScriptL]-\[DoubleStruckCapitalY]^\[ConjugateTranspose].g//Chop//Norm
-\[CapitalGamma].inverseMassU.R^\[ConjugateTranspose].\[DoubleStruckCapitalZ].a - \[CapitalGamma].inverseMassU.R^\[ConjugateTranspose].inverseMassQ.\[CapitalXi]^\[ConjugateTranspose].\[ScriptL]+ \[CapitalGamma].inverseMassU.\[CapitalGamma]^\[ConjugateTranspose].\[Lambda]-\[CapitalGamma].inverseMassU.g//Chop//Norm

\[CapitalXi].inverseMassQ.\[CapitalXi]^\[ConjugateTranspose].\[ScriptL]-h//Chop//Norm

\[CapitalGamma].inverseMassU.\[CapitalGamma]^\[ConjugateTranspose].\[ScriptM]-j//Chop//Norm*)


ConjugateTranspose[\[DoubleStruckCapitalZ]].massQ.\[DoubleStruckCapitalZ].a+ ConjugateTranspose[\[DoubleStruckCapitalZ]].S.\[DoubleStruckCapitalY].b+ConjugateTranspose[\[DoubleStruckCapitalZ]].S.inverseMassU.ConjugateTranspose[\[CapitalGamma]].LinearSolve[\[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]],j]-ConjugateTranspose[\[DoubleStruckCapitalZ]].f//Chop//Norm
-ConjugateTranspose[\[DoubleStruckCapitalY]].ConjugateTranspose[R].\[DoubleStruckCapitalZ].a-ConjugateTranspose[\[DoubleStruckCapitalY]].ConjugateTranspose[R].inverseMassQ.ConjugateTranspose[\[CapitalXi]].LinearSolve[\[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]],h]-ConjugateTranspose[\[DoubleStruckCapitalY]].g//Chop//Norm
\[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]].LinearSolve[\[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]],h]+ \[CapitalXi].inverseMassQ.S.\[DoubleStruckCapitalY].b+\[CapitalXi].inverseMassQ.S.inverseMassU.ConjugateTranspose[\[CapitalGamma]].LinearSolve[\[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]],j]+ \[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]].\[Mu]-\[CapitalXi].inverseMassQ.f//Chop//Norm
-\[CapitalGamma].inverseMassU.ConjugateTranspose[R].\[DoubleStruckCapitalZ].a - \[CapitalGamma].inverseMassU.ConjugateTranspose[R].inverseMassQ.ConjugateTranspose[\[CapitalXi]].LinearSolve[\[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]],h]+ \[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]].\[Lambda]-\[CapitalGamma].inverseMassU.g//Chop//Norm





FF=(-ConjugateTranspose[\[DoubleStruckCapitalZ]].S.inverseMassU.ConjugateTranspose[\[CapitalGamma]].LinearSolve[\[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]],j]+ConjugateTranspose[\[DoubleStruckCapitalZ]].f);
GG=(ConjugateTranspose[\[DoubleStruckCapitalY]].ConjugateTranspose[R].inverseMassQ.ConjugateTranspose[\[CapitalXi]].LinearSolve[\[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]],h]+ConjugateTranspose[\[DoubleStruckCapitalY]].g);


ConjugateTranspose[\[DoubleStruckCapitalY]].ConjugateTranspose[R].\[DoubleStruckCapitalZ].LinearSolve[ConjugateTranspose[\[DoubleStruckCapitalZ]].massQ.\[DoubleStruckCapitalZ],ConjugateTranspose[\[DoubleStruckCapitalZ]].massQ.\[DoubleStruckCapitalZ].a+ConjugateTranspose[\[DoubleStruckCapitalZ]].S.\[DoubleStruckCapitalY].b-FF]//Chop//Norm
Plus[
 ConjugateTranspose[\[DoubleStruckCapitalY]].ConjugateTranspose[R].\[DoubleStruckCapitalZ].LinearSolve[ConjugateTranspose[\[DoubleStruckCapitalZ]].massQ.\[DoubleStruckCapitalZ],ConjugateTranspose[\[DoubleStruckCapitalZ]].massQ.\[DoubleStruckCapitalZ].a+ ConjugateTranspose[\[DoubleStruckCapitalZ]].S.\[DoubleStruckCapitalY].b-FF]
]
Plus[
,
 
]


LinearSolve[


ConjugateTranspose[\[DoubleStruckCapitalZ]].massQ.\[DoubleStruckCapitalZ].a+ ConjugateTranspose[\[DoubleStruckCapitalZ]].S.\[DoubleStruckCapitalY].b+ConjugateTranspose[\[DoubleStruckCapitalZ]].S.inverseMassU.ConjugateTranspose[\[CapitalGamma]].LinearSolve[\[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]],j]-ConjugateTranspose[\[DoubleStruckCapitalZ]].f//Chop//Norm


LinearSolve[ConjugateTranspose[\[DoubleStruckCapitalY]].ConjugateTranspose[R].\[DoubleStruckCapitalZ].LinearSolve[ConjugateTranspose[\[DoubleStruckCapitalZ]].massQ.\[DoubleStruckCapitalZ], ConjugateTranspose[\[DoubleStruckCapitalZ]].S.\[DoubleStruckCapitalY]],ConjugateTranspose[\[DoubleStruckCapitalY]].ConjugateTranspose[R].\[DoubleStruckCapitalZ].LinearSolve[ConjugateTranspose[\[DoubleStruckCapitalZ]].massQ.\[DoubleStruckCapitalZ], FF]+GG]-b//Chop
LinearSolve[
 ConjugateTranspose[\[DoubleStruckCapitalY]].ConjugateTranspose[R].\[DoubleStruckCapitalZ].LinearSolve[ConjugateTranspose[\[DoubleStruckCapitalZ]].massQ.\[DoubleStruckCapitalZ], ConjugateTranspose[\[DoubleStruckCapitalZ]].S.\[DoubleStruckCapitalY]],
 ConjugateTranspose[\[DoubleStruckCapitalY]].ConjugateTranspose[R].\[DoubleStruckCapitalZ].LinearSolve[ConjugateTranspose[\[DoubleStruckCapitalZ]].massQ.\[DoubleStruckCapitalZ], -ConjugateTranspose[\[DoubleStruckCapitalZ]].S.inverseMassU.ConjugateTranspose[\[CapitalGamma]].LinearSolve[\[CapitalGamma].inverseMassU.ConjugateTranspose[\[CapitalGamma]],j]+ConjugateTranspose[\[DoubleStruckCapitalZ]].f]+(ConjugateTranspose[\[DoubleStruckCapitalY]].ConjugateTranspose[R].inverseMassQ.ConjugateTranspose[\[CapitalXi]].LinearSolve[\[CapitalXi].inverseMassQ.ConjugateTranspose[\[CapitalXi]],h]+ConjugateTranspose[\[DoubleStruckCapitalY]].g)
]-b//Chop//Norm


LinearSolve[ConjugateTranspose[\[DoubleStruckCapitalZ]].massQ.\[DoubleStruckCapitalZ],-ConjugateTranspose[\[DoubleStruckCapitalZ]].S.\[DoubleStruckCapitalY].b+FF]-a//Chop





b


SparseArray@(ConjugateTranspose[\[DoubleStruckCapitalY]].ConjugateTranspose[R].\[DoubleStruckCapitalZ].LinearSolve[ConjugateTranspose[\[DoubleStruckCapitalZ]].massQ.\[DoubleStruckCapitalZ], ConjugateTranspose[\[DoubleStruckCapitalZ]].S.\[DoubleStruckCapitalY]])//MatrixPlot


-ConjugateTranspose[\[DoubleStruckCapitalY]].ConjugateTranspose[R].\[DoubleStruckCapitalZ].a-GG//Chop//Norm


p






