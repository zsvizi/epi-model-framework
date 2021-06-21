#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Title:: *)
(*Epidemic modelling framework*)


(* ::Subsection:: *)
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


(* ::Subsection:: *)
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


(* ::Subsection:: *)
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
		Table[Table[StringJoin[comp, ToString[idx]],{idx,Range[chains[comp]]}],{comp, Keys[chains]}]
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
				Table[StringJoin[comp,ToString[1]],{comp, Keys[chains]},{x,Range[Length[inVertices[comp]]]}],
				Table[StringJoin[comp,ToString[chains[comp]]],{comp, Keys[chains]},{x,Range[Length[outVertices[comp]]]}]
			]]
		]
	],
	Flatten[Table[StringJoin[comp,ToString[i-1]]->StringJoin[comp,ToString[i]],{comp,Keys[chains]},{i,Range[2,chains[comp]]}]]
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
