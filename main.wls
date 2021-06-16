#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Title:: *)
(*Epidemic modelling framework*)


transitions={
	{"E","I", \[Alpha]},
	{"I","R", \[Gamma]}
};
transmissions={
	{"I","S","E",\[Beta]}
};


vertices=DeleteDuplicates[Join[
	DeleteDuplicates[Flatten[transitions[[All,1;;2]]]],
	DeleteDuplicates[Flatten[transmissions[[All,1;;3]]]]
]];
edges=Join[
	DeleteDuplicates[Map[Function[x,x[[1]]->x[[2]]],transitions]],
	DeleteDuplicates[Map[Function[x,x[[2]]->x[[3]]],transmissions]]
];
edgeWeights=Join[
	DeleteDuplicates[Flatten[transitions[[All,3]]]],
	DeleteDuplicates[Flatten[transmissions[[All,4]]]]
];


graph=Graph[
	vertices,
	edges,
	EdgeWeight->edgeWeights,
	VertexLabels->"Name",
	EdgeLabels->"EdgeWeight"
]
