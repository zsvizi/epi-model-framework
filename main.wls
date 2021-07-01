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


(*
Function to get vertices connected with outward edges
from the simplified transition graph
*)
outVertices=
Function[
	x,
	VertexOutComponent[simplifiedTransitionGraph,x,1][[2;;]]
];
(*
Function to get vertices connected with inward edges
from the simplified transition graph
*)
inVertices=
Function[
	x,
	VertexInComponent[simplifiedTransitionGraph,x,1][[2;;]]
];


(*
Extend list of verticies for compartments, 
where linear chain trick is used
*)
verticesFull=
Flatten[ReplacePart[vertices,
	Thread[
		Map[Function[v, Position[vertices,v][[1,1]]],Keys[chains]]->
		Table[Table[comp<>ToString[idx],{idx,Range[chains[comp]]}],{comp, Keys[chains]}]
	]
]];


(*
The following changes have to be applied for edge list of simplified transition graph
due to linear chain trick:
- inward edge ? -> X is changed to ? -> X1
- outward edge X -> ? is changed to Xn -> ? (n is the length of the chain for component X)
*)
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
(*
1. Apply previously determined changes to the edge list.
2. Add new directed edges of the chains
*)
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


(*
Find edges, where weights has to be changed
*)
positionsToChange=
Map[
	Function[x, {Position[edges,x][[1,1]]}],
	Flatten[Map[
		Function[comp,If[Length[outVertices[comp]]>0, Thread[comp -> outVertices[comp]],Nothing]],
		Keys[chains]
	]]
(*
1. Apply changes for previously determined weights
2. Add weight for newly added directed edges
*)
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


(* ::Subsection:: *)
(*Generate model equations*)


(*
Create list for transmissions in the full transition graph (i.e. all components of a chain affect the transmission) 
*)
transmissionsFull=
Flatten[Table[
	If[MemberQ[Keys[chains],tr[[1]]],Table[{tr[[1]]<>ToString[i],tr[[2]],tr[[3]]<>ToString[1],tr[[4]]},{i,Range[chains[tr[[1]]]]}],tr],
	{tr,transmissions}
],1];


(*
Create association, where 
- keys are the directed edges, 
- values are the edge weights
*)
p=
Association[
	Table[edgesFull[[eIdx]]->edgeWeightsFull[[eIdx]],{eIdx,Length[edgesFull]}]
];
(*
Loop for initializing edge weights for transmissions
*)
Table[
	p[transmissionsFull[[eIdx,2]]->transmissionsFull[[eIdx,3]]]="",
	{eIdx,Length[transmissionsFull]}];
(*
Loop for determining edge weights for transmissions
*)
Table[
	p[transmissionsFull[[eIdx,2]]->transmissionsFull[[eIdx,3]]]=
		p[transmissionsFull[[eIdx,2]]->transmissionsFull[[eIdx,3]]] <> "+" <> 
		ToString[transmissionsFull[[eIdx,4]]] <> "*" <> (*transmission rate*)
		ToString[transmissionsFull[[eIdx,1]]] <> "[t]", (*infective compartment*)
	{eIdx,Length[transmissionsFull]}];


(*
Association:
- key: component name (from verticesFull)
- value: component symbol
*)
AllComp=
Association[Table[comp->{Symbol[comp]},{comp,verticesFull}]];
(*
Association:
- key: component name (from Keys[AllComp])
- value: component time function
*)
c=
Association[
	Table[comp->Map[Function[x, x[t]],AllComp[comp]],{comp,Keys[AllComp]}]
];
(*
List of component time functions = variables
*)
vars=Flatten[Values[AllComp]];


(*
Association
- key: component name
- value: LHS of initial condition equality
*)
icLHS=
Association[
	Table[
		comp->Map[Function[x, x[ t/;t<=0]],AllComp[comp]],
		{comp,Keys[AllComp]}]
];
icLHSList=Flatten[Values[icLHS]];


(*
Association:
- key: component name
- value: derivative of the component function (symbolic)
*)
LHS=
Association[
	Table[
		comp->Map[Function[x, x'[t]],AllComp[comp]],
		{comp,Keys[AllComp]}]
];
LHSList=Flatten[Values[LHS]];


(*
Function for collecting vertices of inward/outward edges for a component 
*)
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


(*
Function return with merged list of differential equations and initial condition equalities
*)
ode:=
Module[{},
	RHS=Association[Table[comp->0,{comp,Keys[AllComp]}]];
	IC=Association[Table[comp->0,{comp,Keys[icLHS]}]];
	(*
	Loop for filling RHS of equations and IC equalities
	- transmissions "rates" are stored as String (e.g. "\[Beta]*I[t]")
	- simple transitions are stored as Expression (e.g. 3*\[Alpha])
	*)
	Table[
		Module[{},
			RHS[comp]-=Sum[If[p[comp->to][[0]]===String,ToExpression[p[comp->to]],p[comp->to]]*c[comp],{to,outVerticesFull[comp]}];
			RHS[comp]+=Sum[If[p[from->comp][[0]]===String,ToExpression[p[from->comp]],p[from->comp]]*c[from],{from,inVerticesFull[comp]}];
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


fullTransitionGraph


ode//TableForm


paramPattern=Thread[Keys[parameters]->Values[parameters]];
Plot[
	Evaluate[{S[t],I1[t]}/.NDSolve[ode/.paramPattern,vars,{t,0,10}][[1]]],
	{t,0,10}
]
