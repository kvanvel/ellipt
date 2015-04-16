(* ::Package:: *)

ClearAll["Global`*"]


ClearAll["Global`*"]
SetDirectory["~/Workspace/HEAT/src/step-1"];
massU = Import["MassU.mtx"];
massQ = Import["MassQ.mtx"];
inverseMassU = Import["InverseMassU.mtx"];
inverseMassQ = Import["InverseMassQ.mtx"];
traceUdir = Import["TraceUDir.mtx"];
traceQNeu = Import["TraceQNeu.mtx"]; 
totalUFromQ=Import["TotalUFromQ.mtx"];
totalQFromU=Import["TotalQFromU.mtx"];
fluxBoundaryUFromQ=Import["FluxBoundaryUFromQ.mtx"];
fluxBoundaryQFromU=Import["FluxBoundaryQFromU.mtx"];
massElec=Import["MassElec.mtx"];

stiff = Import["StiffQFromU.mtx"];

S=totalQFromU;
minusSstar= totalUFromQ;
R=ConjugateTranspose[totalUFromQ];
\[CapitalGamma]=traceUdir;
\[CapitalXi]=traceQNeu;

Needs["GraphUtilities`"]
{r,c} = MinimumBandwidthOrdering[\[CapitalGamma]];

permutationMatrix[p_List]:=SparseArray@IdentityMatrix[Length[p]][[p]]

zz=SparseArray[Transpose[SparseArray[Chop[NullSpace[\[CapitalGamma][[r,c]]]]].permutationMatrix[c]]];


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






