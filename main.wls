#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Title:: *)
(*Epidemic modelling framework*)


(* ::Subsection::Closed:: *)
(*Input*)


susceptible="S";
transitions={
	{"E","I", \[Alpha]},
	{"I","R", \[Gamma]},
	{"S","V",\[Nu]},
	{"R","S",\[Omega]}
};
transmissions={
	{"I","S","E",\[Beta]},
	{"I","V","E",\[Tau]*\[Beta]}
};
chains=<|
	"E"->3,
	"I"->2
|>;


ic=<|
	"S"->99,
	"I1"->1
|>;
parameters=<|
	\[Beta]->0.5, 
	\[Alpha]->1/5.0, 
	\[Gamma]->1/4.0,
	\[Tau]->0.5,
	\[Omega]->1/100,
	\[Nu]->1/10
|>;


(* ::Subsection::Closed:: *)
(*Generate simplified graph (without chains)*)


vertices=
DeleteDuplicates[Join[
	DeleteDuplicates[Flatten[transitions[[All,1;;2]]]],
	DeleteDuplicates[Flatten[transmissions[[All,1;;3]]]]
]];


edges=
Join[
	DeleteDuplicates[Map[Function[v,v[[1]]->v[[2]]],transitions]],
	DeleteDuplicates[Map[Function[v,v[[2]]->v[[3]]],transmissions]]
];


edgeWeights=
Join[
	DeleteDuplicates[Flatten[transitions[[All,3]]]],
	DeleteDuplicates[Flatten[transmissions[[All,4]]]]
];


simplifiedTransitionGraph=Graph[
	vertices,
	edges,
	EdgeWeight->edgeWeights,
	VertexLabels->"Name",
	EdgeLabels->"EdgeWeight"
]


(* ::Subsection::Closed:: *)
(*Generate graph with chains*)


outVertices=
Function[
	x,
	VertexOutComponent[simplifiedTransitionGraph,x,1][[2;;]]
];
inVertices=
Function[
	x,
	VertexInComponent[simplifiedTransitionGraph,x,1][[2;;]]
];


verticesFull=
Flatten[ReplacePart[vertices,
	Thread[
		Map[Function[v, Position[vertices,v][[1,1]]],Keys[chains]]->
		Table[Table[comp<>ToString[idx],{idx,Range[chains[comp]]}],{comp, Keys[chains]}]
	]
]];


positionsToChange=
Join[
	Map[
		Function[e, {Position[edges,e][[1,1]],2}],
		Flatten[Map[
			Function[comp,If[Length[inVertices[comp]]>0,Thread[inVertices[comp] -> comp],Nothing]],
			Keys[chains]
		]]
	],
	Map[
		Function[e, {Position[edges,e][[1,1]],1}],
		Flatten[Map[
			Function[comp,If[Length[outVertices[comp]]>0,Thread[comp -> outVertices[comp]],Nothing]],
			Keys[chains]
		]]
	]
];
edgesFull=
Flatten[Join[
	ReplacePart[edges,
		Thread[
			positionsToChange->
			Flatten[Join[
				Table[comp<>ToString[1],{comp, Keys[chains]},{x,Range[Length[inVertices[comp]]]}],
				Table[comp<>ToString[chains[comp]],{comp, Keys[chains]},{x,Range[Length[outVertices[comp]]]}]
			]]
		]
	],
	Flatten[Table[comp<>ToString[i-1]->comp<>ToString[i],{comp,Keys[chains]},{i,Range[2,chains[comp]]}]]
]];


positionsToChange=
Map[
	Function[x, {Position[edges,x][[1,1]]}],
	Flatten[Map[
		Function[comp,If[Length[outVertices[comp]]>0, Thread[comp -> outVertices[comp]],Nothing]],
		Keys[chains]
	]]
];
edgeWeightsFull=
Flatten[Join[
	ReplacePart[edgeWeights,
		Thread[
			positionsToChange->
			Flatten[Table[chains[comp]*edgeWeights[[Position[vertices,comp][[1,1]]]],{comp, Keys[chains]},{x,Range[Length[outVertices[comp]]]}]]]
	],
	Flatten[Table[chains[comp]*edgeWeights[[Position[vertices,comp][[1,1]]]],{comp,Keys[chains]},{i,Range[2,chains[comp]]}]]
]];


fullTransitionGraph=
Graph[
	verticesFull,
	edgesFull,
	EdgeWeight->edgeWeightsFull,
	VertexLabels->"Name",
	EdgeLabels->"EdgeWeight"
]


(* ::Subsection::Closed:: *)
(*Generate model equations*)


transmissionsFull=
Flatten[Table[
	If[MemberQ[Keys[chains],tr[[1]]],Table[{tr[[1]]<>ToString[i],tr[[2]],tr[[3]]<>ToString[1],tr[[4]]},{i,Range[chains[tr[[1]]]]}],tr],
	{tr,transmissions}
],1];


p=
Association[
	Table[edgesFull[[eIdx]]->edgeWeightsFull[[eIdx]],{eIdx,Length[edgesFull]}]
];
Table[
	p[transmissionsFull[[eIdx,2]]->transmissionsFull[[eIdx,3]]]="",
	{eIdx,Length[transmissionsFull]}];
Table[
	p[transmissionsFull[[eIdx,2]]->transmissionsFull[[eIdx,3]]]=p[transmissionsFull[[eIdx,2]]->transmissionsFull[[eIdx,3]]]<>"+"<>
		ToString[transmissionsFull[[eIdx,4]]]<>"*"<>ToString[transmissionsFull[[eIdx,1]]]<>"[t]",
	{eIdx,Length[transmissionsFull]}];


AllComp=
Association[Table[comp->{Symbol[comp]},{comp,verticesFull}]];


c=
Association[
	Table[comp->Map[Function[x, x[t]],AllComp[comp]],{comp,Keys[AllComp]}]
];


vars=Flatten[Values[AllComp]];


icLHS=
Association[
	Table[
		comp->Map[Function[x, x[ t/;t<=0]],AllComp[comp]],
		{comp,Keys[AllComp]}]
];


icLHSList=Flatten[Values[icLHS]];


LHS=
Association[
	Table[
		comp->Map[Function[x, x'[t]],AllComp[comp]],
		{comp,Keys[AllComp]}]
];


LHSList=Flatten[Values[LHS]];


outVerticesFull=
Function[
	x,
	VertexOutComponent[fullTransitionGraph,x,1][[2;;]]
];
inVerticesFull=
Function[
	x,
	VertexInComponent[fullTransitionGraph,x,1][[2;;]]
];


ode:=
Module[{},
	RHS=Association[Table[comp->0,{comp,Keys[AllComp]}]];
	IC=Association[Table[comp->0,{comp,Keys[icLHS]}]];
	Table[
		Module[{},
			RHS[comp]=RHS[comp]-Sum[If[p[comp->to][[0]]===String,ToExpression[p[comp->to]],p[comp->to]]*c[comp],{to,outVerticesFull[comp]}];
			RHS[comp]=RHS[comp]+Sum[If[p[from->comp][[0]]===String,ToExpression[p[from->comp]],p[from->comp]]*c[from],{from,inVerticesFull[comp]}];
			If[MemberQ[Keys[ic],comp],IC[comp]=ic[comp],IC[comp]=0];
		],
		{comp,verticesFull}];
	RHSList=Flatten[Values[RHS]];
	equations=Thread[LHSList==RHSList];
	icRHSList=Flatten[Values[IC]];
	initval=Thread[icLHSList==icRHSList];
	Join[equations,initval]
]


(* ::Subsection:: *)
(*Example*)


ode//TableForm


paramPattern=Thread[Keys[parameters]->Values[parameters]];
Plot[
	Evaluate[{S[t],I1[t]}/.NDSolve[ode/.paramPattern,vars,{t,0,10}][[1]]],
	{t,0,10}
]
