(* ::Package:: *)

(* ::Subtitle:: *)
(*QuantumAgebra 2.0*)


(* ::Subsection:: *)
(*General Information*)


(* :Title:  Quantum Algebra *)

(* :Context:  QuantumAlgebra`QuantumAlgebra` *)

(* :Author:  C\[EAcute]sar Guerra
             cguerra@pucp.edu.pe
     
	         Secci\[OAcute]n F\[IAcute]sica, Of. 104        
             Pontificia Universidad Catolica del Peru 
             Av. Universitaria 1801, San Miguel
             Lima-32 PERU. Telf. (511) - 62620000 Anexo: 4108 *)

(* :Summary:  Quantum Algebra is a package to perform Quantum Calculations
              using non commutative algebra. Dirac Notation is implemented. *)

(* :Mathematica Versions: 6.0, 7.0 *)

(* :Package Version:  2.0 *)

(* :Keywords: Quantum Mechanics, Non Commutative Algebra, 
              Kets, Bras, Brakets, Dirac Notation *)

(* :Requirements: To visualize Dirac notations a notebook-based front end 
                  is required. If this is not available, notations can not 
                  be loaded but the package can still be used by using the
                  full form of each object. *)

(* :Warning:     This version of Quantum Algebra doesn't include warning 
                 and error messages. Only usage messages are available. *)

(* :Discussion:  See QuantumAlgebra.nb demo. Several features and some 
                 examples are discussed there. *)



(* ::Subsection:: *)
(* Clear system definitions of Ket, Bra and Braket *)


If[And @@ (NameQ /@ {"System`Bra", "System`Ket", "System`BraKet", "System`TensorProduct"}),
	Unprotect["System`Ket","System`Bra","System`BraKet", "System`TensorProduct"];
	(*ClearAll["System`Ket","System`Bra","System`BraKet"];*)
	Remove["System`Ket","System`Bra","System`BraKet", "System`TensorProduct"]
]


(*Begin and Usage Messages*)
BeginPackage["QuantumAlgebra`QuantumAlgebra`"]

Unprotect @  "QuantumAlgebra`QuantumAlgebra`*";
ClearAll @@ Complement[Names @ "QuantumAlgebra`QuantumAlgebra`*", {"AutoLoadNotationPalette", AutoLoadNavigatorPalette}]
ClearAll @@ Names["QuantumAlgebra`QuantumAlgebra`Private`*"];

BraQ::usage =
"BraQ[expr] gives True if expr is a bra and False otherwise."

Commutator::usage =
"Commutator[A,B] represents the commutator between operators A and B.
 
ALIASES
[ESC] com [ESC] : commutator."

ConstantQ::usage = 
"ConstantQ[expr] gives True if expr is not a quantum object, i.e. it's not an operator, a bra or a ket."

Ket::usage =
"Ket[A,...] is the internal representation of ket with labels A,...
Ket[A,...,{x,...},{}] ket with labels A,... and subscripts x,...
Ket[A,...,{},{y,...}] ket with labels A,... and superscripts y,...
Ket[A,...,{x,...},{y,...}] ket with labels A,...; subscripts x,... and superscripts y,...

ALIASES 
[ESC] ket [ESC]: ket
[ESC] ketd [ESC]: ket with subscript
[ESC] ketu [ESC]: ket with superscript
[ESC] ketdu [ESC]: ket with both scripts"

KetQ::usage =
"KetQ[expr] gives True if expr is a ket and False otherwise."

Operator::usage =
"Operator[A] is the internal representation of operator A. 
Operator[A,{x,...},{}] represents operator A with subscripts x,... 
Operator[A,{},{y,...}] represents operator A with superscripts y...
Operator[A,{x,...},{y,...}] represents operator A with with subscripts x,... and superscripts y,... If y is an integer this is interpreted as a power. 
Operator[...][t,...] represents an operator with aditional parameters t,...
Operator[A_,\"H\"] represents the associated hermitian operator of operator A.
Operator[A_,\"H\",{x,...},{y,...}] subscripts and superscripts can be added as well. 
Operator[1] represents the identity operator. 
Operator[0] represents the null operator.

ALIASES
[ESC] op [ESC]: operator.
[ESC] opd [ESC]: operator with subscript
[ESC] opu [ESC]: operator with superscript
[ESC] opdu [ESC]: operator with both scripts 
[ESC] hop [ESC]: hermitian operator."

OperatorQ::usage = "OperatorQ[expr] gives True if expr is an operator and False otherwise."

QuantumObjectQ::usage = "QuantumObjectQ[expr] gives True if expr is a quantum object, i.e contains operators, bras or kets."

QATrace::usage = "QATrace[op, base] computes the trace of op in basis base";

QAMatrix::usage = "QAMatrix[op, base] computes matrix representation of op in basis base";

QADirac::usage = "QADirac[op, base] computes Dirac representation of op in basis base";

QACompute::usage = "";

(******* OLD *********)

Bra::usage =
"Bra is the head for a Bra object. Some full forms that 
\ncan be interpreted are:
\n\t- Bra[A__]
\n\t- Bra[A__,{x__},{}], 
\n\t- Bra[A__,{},{y__}],
\n\t- Bra[A__,{x__},{y__}]
\nwhere A labels the Ket and x and y are subsripts and
\nsuperscripts respectively.
\nALIASES: 
\n[ESC] bra [ESC]   : bra
\n[ESC] dbra [ESC]  : bra with subscript
\n[ESC] ubra [ESC]  : bra with superscript
\n[ESC] dubra [ESC] : bra with both scripts.
\nBra is also used as a tag for TaxBoxes that guide the 
\ninterpretation of the internal structure of bra notation."

Braket::usage =
"Braket is the head for a Braket object. Some full forms that 
\ncan be interpreted are:
\n\t- Braket[{A__},{B__}]
\n\t- Braket[{A__,{x__},{}},{B__}], 
\n\t- Braket[{A__,{},{y__}},{B__}],
\n\t- Braket[{A__,{x__},{y__}},{B__,{x__},{Y__}}]
\nwhere A and B label the Braket and x and y are subsripts 
\nand superscripts respectively. Of course more combinations
\nare posible.
\nALIASES: 
\n[ESC] braket [ESC]     : braket
\n[ESC] braketd [ESC]    : braket with subscript
\n[ESC] ubraket [ESC]    : braket with superscript
\n[ESC] dubraketdu [ESC] : braket with both scripts.
\nMore combinations are possible using u(up) and d(down).
\nBraket is also used as a tag for TaxBoxes that guide the 
\ninterpretation of the internal structure of braket notation."

Hermitian::usage = 
"Hermitian[expr] performs the hermitian operation on expr.
\nA expression with a dagger as a superscript is interpreted
\nas the hermiatian operation on expression"

QAPowerExpand::usage =
"QAPowerExpand[expr] expands all powers in expr as non 
\ncommutative products."

QAPowerContract::usage =
"QAPowerContract[expr] is intented to contract terms like: 
\nA^n ** A^m to A^(n+m) where A is an operator and n and m
\nsimple exponents."

QAProductExpand::usage =
"QAProductExpand[expr] expands all non commutative products of
\nconstant terms and quantum objects in expr as simple products
\nof constant terms and non commutative products of quantum
\nobjects."

QACommutatorExpand::usage =
"QACommutatorExpand[expr] expands all commutators in expr."

QACommutatorContract::usage =
"QACommutatorContract[expr] contracts all expanded commutators 
\nin expr. 
\nQACommutatorContract[expr,{A,B}] contracts all expanded form 
\nof Commutator[A,B] in expr. 
\nQACommutatorContract[expr,{{A,B},...}] performs several
\ncommutators contracts."

QAExpandAll::usage =
"QAExpandAll[expr] performs in expr the following commands:
\nQAPowerExpand, QAProductExpand, and QACommutatorExpand in that 
\norder"

QACollect::usage =
" "

QASeries::usage =
"QASeries[expr,{A,B,n}] expands the operator expr depending 
\nonly on A, in power series of A around the operator B up to
\nn order. The operator B must be of the form c Operator[1] or
\nOperator[0] (cero). 
\nQASeries[Exp[A] ** B ** Exp[-A],n] expands the first argument
\nas a series expansion in terms of A and B operators up to n 
\nterms.
\nQASeries does not return a series object, but a polinomial 
\non operaotor A."

QAOrdering::usage =
"QAOrdering[expr_, list_] transform expr which is formed by 
\nproduct of operators to an equivalent expression in which the
\norder of the products is changed. The new order will be as 
\ngiven in list from left to right."

BCHRelation::usage =
"BCHRelation[expr] use the Baker Campbell Hausdorff relation
\nto expand all operators of the form Exp[A+B] in expr where 
\nA and B are operators. BCHRelation automatically checks if 
\nA and B commute with their commutator."

GramSchmidtVectors::usage =
" "

QASet::usage =
"QASet[ expr1 = expr2] and QASet[expr1 := expr2] where expr1 and
\nexpr2 involve quantum objects like operators, bras, kets...,
\nfirst evaluate the Set and SetDelayed operations as given in 
\nits argument. After that, QASet try to set up definitions that
\ncan be deduced from the first relation. 
\nFor example if we evaluate QASet[ Commutator[A,B] = C ] then 
\nthe argument is evaluated in the usual way. But from this set 
\nit can be deduced that Commutator[B,A] = -C, so QASet also 
\ndoes this assignation.
\nQASet can handle the following definitions: 
\n\t- QASet[ Com[A,B] = expr ]
\n\t- QASet[ Com[A,Com[A,B]] = expr ]
\n\t- QASet[ Com[Com[A,B],Com[C,D]] = expr ]
\n\t- QASet[ Com[A,Com[A,B]] = 0, Com[B,Com[A,B]] = 0 ]
\n\t- QASet[ A ** Ket[n] = a Ket[n] ]
\nwhere A, B, C, D are operators and Com means Commutator. Also,
\nit can be used := instead of = if neccesary. For a full dicussion
\nsee QuantumAlgebra.nb Section 8. "

QAClear::usage =
"QAClear[ lhs =. ] removes all definitions made by QASet associated 
\nwith the definition of lhs. See QASet. "

SuperDagger::usage =
"Perform same operationa as Hermitian."

SuperStar::usage =
"Perform same operation as conjugate"

CircleTimes::usage =
" "

BraArgs::usage =
"BraArgs is a tag for TaxBoxes that guide the interpretation 
\nof the internal structure of bra notation."

KetArgs::usage =
"KetArgs is a tag for TaxBoxes that guide the interpretation 
\nof the internal structure of ket notation."

BraKetArgs::usage =
"BraKetArgs is a tag for TaxBoxes that guide the interpretation 
\nof the internal structure of braket notation."

UDScript::usage =
"UDScripts is a tag for TaxBoxes that guide the interpretation  
\nof the internal structure of bras, kets, and brakets notation."

TensorProduct::usage = "TensorProduct[opts__] defines the tensor product of operators opts"

DoTensorProduct::usage = ""

QACommute::usage = ""

QAClearCommute::usage = ""

ExtendSystem::usage = ""

AutoLoadNotationPalette::usage = "Flag to load the notation palette"

AutoLoadNavigatorPalette::usage = "Flag to load the navigator palette"


Begin["`Private`"]

SetAttributes[{Hermitian, OperatorQ, KetQ, BraQ, ConstantQ, QuantumObjectQ, 
	QAProductExpand, QAPowerExpand, QACommutatorExpand, QAExpandAll}, Listable]

SetAttributes[DoTensorProduct,{HoldAll}]


protected = Unprotect[NonCommutativeMultiply, Exp, Times, Power, Conjugate, SubscriptBox, RowBox, FullForm];

(* ::Subsection:: *)
(*Defining the type of objects*)


ConstantQ[expr_] := FreeQ[expr,Operator] && FreeQ[expr,Ket] && FreeQ[expr,Bra]

QuantumObjectQ[expr_] := !ConstantQ[expr]

OperatorQ[expr_] := !FreeQ[expr, Operator]      

KetQ[expr_] := !FreeQ[expr, Ket]

BraQ[expr_] := !FreeQ[expr, Bra]


(* ::Subsection:: *)
(*Algebra for NCM*)


P_ ** M_List := P ** # & /@ M /; ArrayDepth[M] == 1

P_ ** M_List := {#1, P ** #2} & @@@ M /; ArrayDepth[M] == 2

0 ** P_ = 0

P_ ** 0 = 0

Operator[0] ** P_ := 0
 
P_ ** Operator[0] := 0
 
Operator/: Operator[0] + P_ := P
 
Operator/: Operator[0] a_ := 0
 
Operator[1] ** P_ := P 
 
P_ ** Operator[1] := P 
 
P___ ** (Q_Plus) ** R___ := Map[(P ** # ** R)&, Q]

P___ ** (c_? ConstantQ Q_) ** R___ := c P ** Q ** R /; ({P,R}=!={})



(* ::Subsection:: *)
(*Algebra for Tensor Product*)


Ket[x__Subscript] := Ket[Sequence@@Sort[{x},OrderedQ[{#1[[2]],#2[[2]]}]&]]/;
	(!OrderedQ[{x},OrderedQ[{#1[[2]],#2[[2]]}]&])

Bra[x__Subscript] := Bra[Sequence@@Sort[{x},OrderedQ[{#1[[2]],#2[[2]]}]&]]/;
	(!OrderedQ[{x},OrderedQ[{#1[[2]],#2[[2]]}]&])

Unprotect[TensorProduct];	

SetAttributes[TensorProduct, {Flat, OneIdentity}]
	
TensorProduct[P___, Q_Plus, R___] := Map[TensorProduct[P, #, R]&, Q]

TensorProduct[P___, (c_? ConstantQ Q_),  R___] := c TensorProduct[P, Q, R] /; ({P,R}=!={})

TensorProduct[lista__List]:= Module[{prod},
  	prod = Flatten[Outer[TensorProduct, lista]]
	]/; Length[{lista}] >= 2

TensorProduct[ g:(Ket[__Subscript]..)] := Ket @@(Sequence@@#&/@{g}) /; Length[{g}]>1

TensorProduct[ g:(Ket[__Subscript,{x___},{y___}]..)] := 
	Ket[Sequence@@Flatten[Drop[#,-2]&/@(List@@#&/@{g})],{x},{y}] /; Length[{g}]>1
	
TensorProduct[ g:(Bra[__Subscript]..)] := Bra @@(Sequence@@#&/@{g}) /; Length[{g}]>1

TensorProduct[ g:(Bra[__Subscript,{x___},{y___}]..)] := 
	Bra[Sequence@@Flatten[Drop[#,-2]&/@(List@@#&/@{g})],{x},{y}] /; Length[{g}]>1
	
TensorProduct[ g:((Bra|Ket)[__Subscript,{__},{___}]..)] := TensorProduct @@ (TensorProduct @@ # & /@ GatherBy[{g}, #[[2]] &]) /; Length[{g}]>1 && DeleteDuplicates[{g}[[All,2]]] != {g}[[All,2]]

TensorProduct[NonCommutativeMultiply[x_Ket,y_Bra],
	NonCommutativeMultiply[u_Ket,v_Bra]]:=
		NonCommutativeMultiply[TensorProduct[x,u],TensorProduct[y,v]];

TensorProduct[NonCommutativeMultiply[TensorProduct[x__Ket],TensorProduct[y__Bra]],
	NonCommutativeMultiply[u_Ket,v_Bra]]:=
		NonCommutativeMultiply[TensorProduct[x,u],TensorProduct[y,v]];

TensorProduct[NonCommutativeMultiply[x_Ket,y_Bra],
	NonCommutativeMultiply[TensorProduct[u__Ket],TensorProduct[v__Bra]]]:=
		NonCommutativeMultiply[TensorProduct[x,u],TensorProduct[y,v]];

TensorProduct[NonCommutativeMultiply[TensorProduct[x__Ket],TensorProduct[y__Bra]],
	NonCommutativeMultiply[TensorProduct[u__Ket],TensorProduct[v__Bra]]]:=
		NonCommutativeMultiply[TensorProduct[x,u],TensorProduct[y,v]];

(* SetAttributes[TensorProduct, {Orderless}] *)

(*TensorProduct[a_Ket] := (Print["1"];a)

TensorProduct[a_Bra] := a*)

(* ::Subsection:: *)
(*Tensor Product Generation*)


DoTensorProduct[expr_,{n_}] := TensorProduct @@ Table[expr,{n}]

DoTensorProduct[expr_,{i_,n_}] := TensorProduct@@Table[expr,{i,n}]

DoTensorProduct[expr_,{i_,n_,n_}] := First@Table[expr,{i,n,n}]

DoTensorProduct[expr_,{i_,n_,m_}] := TensorProduct@@Table[expr,{i,n,m}]

DoTensorProduct[expr_,{i_,n_,m_,d_}] := TensorProduct@@Table[expr,{i,n,m,d}]

DoTensorProduct[expr_,{i_,{j__}}] := TensorProduct@@Table[expr,{i,{j}}]


(* ::Subsection:: *)
(*Algebra for braketing*)


Bra /: y__ ** Bra[x___, Subscript[a_, n_], z___] ** Ket[u___, Subscript[b_, n_], v___] ** w__:= 
  Which[
	Length[{x, z}] >= 1 && Length[{u, v}] >= 1, 
			Braket[{Subscript[a, n]}, {Subscript[b, n]}] * y ** Bra[x, z] ** Ket[u, v] ** w, 
	Length[{x, z}] == 0 && Length[{u, v}] >= 1, 
			Braket[{Subscript[a, n]},{Subscript[b, n]}] * y ** Ket[u, v] ** w, 
	Length[{x, z}] >= 1 && Length[{u, v}] == 0, 
			Braket[{Subscript[a, n]}, {Subscript[b, n]}] * y ** Bra[x, z] ** w,
	Length[{x, z}] == 0 && Length[{u, v}] == 0, 
			Braket[{Subscript[a, n]}, {Subscript[b, n]}] y ** w]  /;
		And@@(FreeQ[#,List]&/@{Sequence[z,v]})

Bra /: Bra[x___, Subscript[a_, n_], z___] ** Ket[u___, Subscript[b_, n_], v___] ** y_:= 
  Which[
	Length[{x, z}] >= 1 && Length[{u, v}] >= 1, 
			Braket[{Subscript[a, n]}, {Subscript[b, n]}]* Bra[x, z] ** Ket[u, v] ** y, 
	Length[{x, z}] == 0 && Length[{u, v}] >= 1, 
			Braket[{Subscript[a, n]},{Subscript[b, n]}] * Ket[u, v] ** y, 
	Length[{x, z}] >= 1 && Length[{u, v}] == 0, 
			Braket[{Subscript[a, n]}, {Subscript[b, n]}] * Bra[x, z] ** y,
	Length[{x, z}] == 0 && Length[{u, v}] == 0, 
			Braket[{Subscript[a, n]}, {Subscript[b, n]}] * y] /;
		True && And@@(FreeQ[#,List]&/@{Sequence[z,v]})

Bra /: y_ ** Bra[x___, Subscript[a_, n_], z___] ** Ket[u___, Subscript[b_, n_], v___] := 
  Which[
	Length[{x, z}] >= 1 && Length[{u, v}] >= 1, 
			Braket[{Subscript[a, n]}, {Subscript[b, n]}] * y ** Bra[x, z] ** Ket[u, v], 
	Length[{x, z}] == 0 && Length[{u, v}] >= 1, 
			Braket[{Subscript[a, n]},{Subscript[b, n]}] * y ** Ket[u, v], 
	Length[{x, z}] >= 1 && Length[{u, v}] == 0, 
			Braket[{Subscript[a, n]}, {Subscript[b, n]}] * y ** Bra[x, z],
	Length[{x, z}] == 0 && Length[{u, v}] == 0, 
			Braket[{Subscript[a, n]}, {Subscript[b, n]}] * y] /;
		True && And@@(FreeQ[#,List]&/@{Sequence[z,v]})

Bra /: Bra[x___, Subscript[a_, n_], z___] ** Ket[u___, Subscript[b_, n_], v___] := 
  Which[
	Length[{x, z}] >= 1 && Length[{u, v}] >= 1, 
			Braket[{Subscript[a, n]}, {Subscript[b, n]}] * Bra[x, z] ** Ket[u, v], 
	Length[{x, z}] == 0 && Length[{u, v}] >= 1, 
			Braket[{Subscript[a, n]},{Subscript[b, n]}] * Ket[u, v], 
	Length[{x, z}] >= 1 && Length[{u, v}] == 0, 
			Braket[{Subscript[a, n]}, {Subscript[b, n]}] * Bra[x, z],
	Length[{x, z}] == 0 && Length[{u, v}] == 0, 
			Braket[{Subscript[a, n]}, {Subscript[b, n]}]] /;
		And@@(FreeQ[#,List]&/@{Sequence[z,v]})

Bra /: y_ ** Bra[x___, Subscript[a_, n_], z___, {p___}, {q___}] ** Ket[u___, Subscript[b_, n_], v___, {p___}, {q___}] ** w_:= 
  Which[
	Length[{x, z}] >= 1 && Length[{u, v}] >= 1, 
			Braket[{Subscript[a, n],{p},{q}}, {Subscript[b, n], {p}, {q}} ] * y ** Bra[x, z, {p}, {q}] ** Ket[u, v, {p}, {q}] ** w, 
	Length[{x, z}] == 0 && Length[{u, v}] >= 1, 
			Braket[{Subscript[a, n], {p}, {q}},{Subscript[b, n], {p}, {q}}] * y ** Ket[u, v, {p}, {q}] ** w, 
	Length[{x, z}] >= 1 && Length[{u, v}] == 0, 
			Braket[{Subscript[a, n],{p},{q}}, {Subscript[b, n], {p}, {q}}] * y ** Bra[x, z, {p}, {q}] ** w,
	Length[{x, z}] == 0 && Length[{u, v}] == 0, 
			Braket[{Subscript[a, n],{p},{q}}, {Subscript[b, n], {p}, {q}}] * y ** w] /; True

Bra /: Bra[x___, Subscript[a_, n_], z___, {p___}, {q___}] ** Ket[u___, Subscript[b_, n_], v___, {p___}, {q___}] ** y_ := 
  Which[
	Length[{x, z}] >= 1 && Length[{u, v}] >= 1, 
			Braket[{Subscript[a, n],{p},{q}}, {Subscript[b, n], {p}, {q}} ] * Bra[x, z, {p}, {q}] ** Ket[u, v, {p}, {q}] ** y, 
	Length[{x, z}] == 0 && Length[{u, v}] >= 1, 
			Braket[{Subscript[a, n], {p}, {q}},{Subscript[b, n], {p}, {q}}] * Ket[u, v, {p}, {q}] ** y, 
	Length[{x, z}] >= 1 && Length[{u, v}] == 0, 
			Braket[{Subscript[a, n],{p},{q}}, {Subscript[b, n], {p}, {q}}] * Bra[x, z, {p}, {q}] ** y,
	Length[{x, z}] == 0 && Length[{u, v}] == 0, 
			Braket[{Subscript[a, n],{p},{q}}, {Subscript[b, n], {p}, {q}}] * y ] /; True

Bra /: y_ ** Bra[x___, Subscript[a_, n_], z___, {p___}, {q___}] ** Ket[u___, Subscript[b_, n_], v___, {p___}, {q___}] := 
  Which[
	Length[{x, z}] >= 1 && Length[{u, v}] >= 1, 
			Braket[{Subscript[a, n],{p},{q}}, {Subscript[b, n],{p},{q}}] * y ** Bra[x, z,{p},{q}] ** Ket[u, v, {p}, {q}], 
	Length[{x, z}] == 0 && Length[{u, v}] >= 1, 
			Braket[{Subscript[a, n],{p},{q}},{Subscript[b, n],{p},{q}}] * y ** Ket[u, v, {p}, {q}], 
	Length[{x, z}] >= 1 && Length[{u, v}] == 0, 
			Braket[{Subscript[a, n],{p},{q}}, {Subscript[b, n],{p},{q}}] * y ** Bra[x, z, {p}, {q}],
	Length[{x, z}] == 0 && Length[{u, v}] == 0, 
			Braket[{Subscript[a, n],{p},{q}}, {Subscript[b, n], {p}, {q}}] * y ]/; True

Bra /: Bra[x___, Subscript[a_, n_], z___, {p___}, {q___}] ** Ket[u___, Subscript[b_, n_], v___, {p___}, {q___}] := 
  Which[
	Length[{x, z}] >= 1 && Length[{u, v}] >= 1, 
			Braket[{Subscript[a, n],{p},{q}}, {Subscript[b, n],{p},{q}}] * Bra[x, z,{p},{q}] ** Ket[u, v,{p},{q}], 
	Length[{x, z}] == 0 && Length[{u, v}] >= 1, 
			Braket[{Subscript[a, n],{p},{q}},{Subscript[b, n],{p},{q}}] * Ket[u, v,{p},{q}], 
	Length[{x, z}] >= 1 && Length[{u, v}] == 0, 
			Braket[{Subscript[a, n],{p},{q}}, {Subscript[b, n],{p},{q}}] * Bra[x, z, {p}, {q}],
	Length[{x, z}] == 0 && Length[{u, v}] == 0, 
			Braket[{Subscript[a, n],{p},{q}}, {Subscript[b, n], {p}, {q}}] ]

Bra /: Bra[x___, Subscript[a_, n_], z___] ** (f_Ket) ** Ket[u___, Subscript[b_, n_], v___] := 
  Which[
	Length[{x, z}] >= 1 && Length[{u, v}] >= 1, 
			Braket[{Subscript[a, n]}, {Subscript[b, n]}] * Bra[x, z] ** f ** Ket[u, v], 
	Length[{x, z}] == 0 && Length[{u, v}] >= 1, 
			Braket[{Subscript[a, n]},{Subscript[b, n]}] * f ** Ket[u, v], 
	Length[{x, z}] >= 1 && Length[{u, v}] == 0, 
			Braket[{Subscript[a, n]}, {Subscript[b, n]}] * Bra[x, z] ** f,
	Length[{x, z}] == 0 && Length[{u, v}] == 0, 
			Braket[{Subscript[a, n]}, {Subscript[b, n]}] * f] /;
		And@@(FreeQ[#,List]&/@{Sequence[z,v]})

Bra /: Bra[x___, Subscript[a_, n_], z___] ** (f__Ket) ** Ket[u___, Subscript[b_, n_], v___] := 
  Which[
	Length[{x, z}] >= 1 && Length[{u, v}] >= 1, 
			Braket[{Subscript[a, n]}, {Subscript[b, n]}] * Bra[x, z] ** f ** Ket[u, v], 
	Length[{x, z}] == 0 && Length[{u, v}] >= 1, 
			Braket[{Subscript[a, n]},{Subscript[b, n]}] * f ** Ket[u, v], 
	Length[{x, z}] >= 1 && Length[{u, v}] == 0, 
			Braket[{Subscript[a, n]}, {Subscript[b, n]}] * Bra[x, z] ** f,
	Length[{x, z}] == 0 && Length[{u, v}] == 0, 
			Braket[{Subscript[a, n]}, {Subscript[b, n]}] * NonCommutativeMultiply[f]] /;
		And@@(FreeQ[#,List]&/@{Sequence[z,v]})

Bra /: Bra[x___, Subscript[a_, n_], z___,{c___},{d___}] ** (f_Ket) ** Ket[u___, Subscript[b_, n_], v___,{c___},{d___}] := 
  Which[
	Length[{x, z}] >= 1 && Length[{u, v}] >= 1, 
			Braket[{Subscript[a, n],{c},{d}}, {Subscript[b, n],{c},{d}}]* Bra[x, z,{c},{d}] ** f ** Ket[u, v,{c},{d}], 
	Length[{x, z}] == 0 && Length[{u, v}] >= 1, 
			Braket[{Subscript[a, n],{c},{d}},{Subscript[b, n],{c},{d}}] * f ** Ket[u, v, {c}, {d}], 
	Length[{x, z}] >= 1 && Length[{u, v}] == 0, 
			Braket[{Subscript[a, n],{c},{d}}, {Subscript[b, n],{c},{d}}] * Bra[x, z, {c}, {d}] ** f,
	Length[{x, z}] == 0 && Length[{u, v}] == 0, 
			Braket[{Subscript[a, n],{c},{d}}, {Subscript[b, n], {c}, {d}}] f ]
 

Bra /: Bra[x___, Subscript[a_, n_], z___,{c___},{d___}] ** (f__Ket) ** Ket[u___, Subscript[b_, n_], v___,{c___},{d___}] := 
  Which[
	Length[{x, z}] >= 1 && Length[{u, v}] >= 1, 
			Braket[{Subscript[a, n],{c},{d}}, {Subscript[b, n],{c},{d}}]* Bra[x, z,{c},{d}] ** Ket[u, v,{c},{d}], 
	Length[{x, z}] == 0 && Length[{u, v}] >= 1, 
			Braket[{Subscript[a, n],{c},{d}},{Subscript[b, n],{c},{d}}] * Ket[u, v, {c}, {d}], 
	Length[{x, z}] >= 1 && Length[{u, v}] == 0, 
			Braket[{Subscript[a, n],{c},{d}}, {Subscript[b, n],{c},{d}}] * Bra[x, z, {c}, {d}],
	Length[{x, z}] == 0 && Length[{u, v}] == 0, 
			Braket[{Subscript[a, n],{c},{d}}, {Subscript[b, n], {c}, {d}}] * NonCommutativeMultiply[f] ]

Bra /: Bra[x___, Subscript[a_, n_], z___] ** (f_Bra) ** Ket[u___, Subscript[b_, n_], v___] := 
  Which[
	Length[{x, z}] >= 1 && Length[{u, v}] >= 1, 
			Braket[{Subscript[a, n]}, {Subscript[b, n]}] * Bra[x, z] ** f ** Ket[u, v], 
	Length[{x, z}] == 0 && Length[{u, v}] >= 1, 
			Braket[{Subscript[a, n]},{Subscript[b, n]}] * f ** Ket[u, v], 
	Length[{x, z}] >= 1 && Length[{u, v}] == 0, 
			Braket[{Subscript[a, n]}, {Subscript[b, n]}] * Bra[x, z] ** f,
	Length[{x, z}] == 0 && Length[{u, v}] == 0, 
			Braket[{Subscript[a, n]}, {Subscript[b, n]}] * f ]  /;
		And@@(FreeQ[#,List]&/@{Sequence[z,v]})

Bra /: Bra[x___, Subscript[a_, n_], z___] ** (f__Bra) ** Ket[u___, Subscript[b_, n_], v___] := 
  Which[
	Length[{x, z}] >= 1 && Length[{u, v}] >= 1, 
			Braket[{Subscript[a, n]}, {Subscript[b, n]}]* Bra[x, z] ** f ** Ket[u, v], 
	Length[{x, z}] == 0 && Length[{u, v}] >= 1, 
			Braket[{Subscript[a, n]},{Subscript[b, n]}] * f ** Ket[u, v], 
	Length[{x, z}] >= 1 && Length[{u, v}] == 0, 
			Braket[{Subscript[a, n]}, {Subscript[b, n]}] * Bra[x, z] ** f,
	Length[{x, z}] == 0 && Length[{u, v}] == 0, 
			Braket[{Subscript[a, n]}, {Subscript[b, n]}] * NonCommutativeMultiply[f]] /;
		And@@(FreeQ[#,List]&/@{Sequence[z,v]})

Bra /: Bra[x___, Subscript[a_, n_], z___,{c___},{d___}] ** (f_Bra) ** Ket[u___, Subscript[b_, n_], v___,{c___},{d___}] := 
  Which[
	Length[{x, z}] >= 1 && Length[{u, v}] >= 1, 
			Braket[{Subscript[a, n],{c},{d}}, {Subscript[b, n],{a},{b}}] * Bra[x, z,{c},{d}] ** f ** Ket[u, v,{c},{d}], 
	Length[{x, z}] == 0 && Length[{u, v}] >= 1, 
			Braket[{Subscript[a, n],{c},{d}},{Subscript[b, n],{c},{d}}] * f ** Ket[u, v, {c},{d}], 
	Length[{x, z}] >= 1 && Length[{u, v}] == 0, 
			Braket[{Subscript[a, n],{c},{d}}, {Subscript[b, n],{c},{d}}] * Bra[x, z,{c},{d}] ** f,
	Length[{x, z}] == 0 && Length[{u, v}] == 0, 
			Braket[{Subscript[a, n],{c},{d}}, {Subscript[b, n], {c}, {d}}] * f]
 

Bra /: Bra[x___, Subscript[a_, n_], z___,{c___},{d___}] ** (f__Bra) ** Ket[u___, Subscript[b_, n_], v___,{c___},{d___}] := 
  Which[
	Length[{x, z}] >= 1 && Length[{u, v}] >= 1, 
			Braket[{Subscript[a, n],{c},{d}}, {Subscript[b, n],{c},{d}}]* Bra[x, z,{c},{d}] ** f ** Ket[u, v,{c},{d}], 
	Length[{x, z}] == 0 && Length[{u, v}] >= 1, 
			Braket[{Subscript[a, n],{c},{d}},{Subscript[b, n],{c},{d}}] * f ** Ket[u, v,{c},{d}], 
	Length[{x, z}] >= 1 && Length[{u, v}] == 0, 
			Braket[{Subscript[a, n],{c},{d}}, {Subscript[b, n],{c},{d}}] * Bra[x, z,{c},{d}] ** f,
	Length[{x, z}] == 0 && Length[{u, v}] == 0, 
			Braket[{Subscript[a, n],{c},{d}}, {Subscript[b, n], {c}, {d}}] * NonCommutativeMultiply[f] ]
 

Bra /: Bra[x___, Subscript[a_, n_], z___] ** (f__Bra) ** (g__Ket) ** Ket[u___, Subscript[b_, n_], v___] := 
  Which[
	Length[{x, z}] >= 1 && Length[{u, v}] >= 1, 
			Braket[{Subscript[a, n]}, {Subscript[b, n]}] * Bra[x, z] ** f ** g ** Ket[u, v], 
	Length[{x, z}] == 0 && Length[{u, v}] >= 1, 
			Braket[{Subscript[a, n]},{Subscript[b, n]}] * f ** g ** Ket[u, v], 
	Length[{x, z}] >= 1 && Length[{u, v}] == 0, 
			Braket[{Subscript[a, n]},{Subscript[b, n]}] * Bra[x, z] ** f ** g, 
	Length[{x, z}] == 0 && Length[{u, v}] == 0, 
			Braket[{Subscript[a, n]}, {Subscript[b, n]}] * f ** g]  /;
		And@@(FreeQ[#,List]&/@{Sequence[z,v]})

Bra /: Bra[x___, Subscript[a_, n_], z___,{c___},{d___}] ** (f__Bra) ** (g__Ket) ** Ket[u___, Subscript[b_, n_], v___,{c___},{d___}] := 
  Which[
	Length[{x, z}] >= 1 && Length[{u, v}] >= 1, 
			Braket[{Subscript[a, n],{c},{d}}, {Subscript[b, n],{c},{d}}]* Bra[x, z,{c},{d}] ** f ** g ** Ket[u, v,{c},{d}], 
	Length[{x, z}] == 0 && Length[{u, v}] >= 1, 
			Braket[{Subscript[a, n],{c},{d}},{Subscript[b, n],{c},{d}}] * f ** g ** Ket[u, v,{c},{d}], 
	Length[{x, z}] >= 1 && Length[{u, v}] == 0, 
			Braket[{Subscript[a, n],{c},{d}}, {Subscript[b, n],{c},{d}}] * Bra[x, z,{c},{d}] ** f ** g,
	Length[{x, z}] == 0 && Length[{u, v}] == 0, 
			Braket[{Subscript[a, n],{c},{d}}, {Subscript[b, n], {c}, {d}}] f ** g]


Braket/: Conjugate[Braket[{m__},{n__}]] := Braket[{n},{m}];

Bra /: y__ ** Bra[u__] ** Ket[v__] ** z__ := Braket[{u},{v}] y ** z /; 
		 And@@(FreeQ[#,List]&/@{Sequence[u,v]}) && FreeQ[{u,v},Subscript]

Bra /: Bra[u__] ** Ket[v__] ** y_ :=  y Braket[{u},{v}] /; 
		True  && And@@(FreeQ[#,List]&/@{Sequence[u,v]}) && FreeQ[{u,v},Subscript]

Bra /: y_ ** Bra[u__] ** Ket[v__] := y Braket[{u},{v}] /; 
		True && And@@(FreeQ[#,List]&/@{Sequence[u,v]}) && FreeQ[{u,v},Subscript]

Bra /: Bra[u__] ** Ket[v__] :=  Braket[{u},{v}] /; 
		 And@@(FreeQ[#,List]&/@{Sequence[u,v]}) && FreeQ[{u,v},Subscript]

Bra /: y__ ** Bra[u__,{a___},{b___}] ** Ket[v__,{a___},{b___}] ** z__ := 
		Braket[{u,{a},{b}},{v,{a},{b}}] y ** z /; True  && FreeQ[{u,v},Subscript]

Bra /: Bra[u__,{a___},{b___}] ** Ket[v__,{a___},{b___}] ** y_ := 
		y Braket[{u,{a},{b}},{v,{a},{b}}] /; True  && FreeQ[{u,v},Subscript]

Bra /: y_ ** Bra[u__,{a___},{b___}] ** Ket[v__,{a___},{b___}] := 
		y Braket[{u,{a},{b}},{v,{a},{b}}] /; True  && FreeQ[{u,v},Subscript]

Bra /: Bra[u__,{a___},{b___}] ** Ket[v__,{a___},{b___}] := 
		Braket[{u,{a},{b}},{v,{a},{b}}] /; FreeQ[{u,v},Subscript]

Bra /: Bra[u__] ** (f_Ket) ** Ket[v__] := Braket[{u},{v}] * f /; 
		 And@@(FreeQ[#,List]&/@{Sequence[u,v]}) && FreeQ[{u,v},Subscript]

Bra /: Bra[u__] ** (f__Ket) ** Ket[v__] := Braket[{u},{v}] * NonCommutativeMultiply[f] /; 
		 And@@(FreeQ[#,List]&/@{Sequence[u,v]}) && FreeQ[{u,v},Subscript]

Bra /: Bra[u__,{a___},{b___}] ** (f_Ket) ** Ket[v__,{a___},{b___}] := 
		Braket[{u,{a},{b}},{v,{a},{b}}] * f /; FreeQ[{u,v},Subscript]

Bra /: Bra[u__,{a___},{b___}] ** (f__Ket) ** Ket[v__,{a___},{b___}] := 
		Braket[{u,{a},{b}},{v,{a},{b}}] * NonCommutativeMultiply[f] /; FreeQ[{u,v},Subscript]

Bra /: Bra[u__] ** (f_Bra) ** Ket[v__] := Braket[{u},{v}] * f /; 
		 And@@(FreeQ[#,List]&/@{Sequence[u,v]}) && FreeQ[{u,v},Subscript]

Bra /: Bra[u__] ** (f__Bra) ** Ket[v__] := Braket[{u},{v}] * NonCommutativeMultiply[f] /; 
		 And@@(FreeQ[#,List]&/@{Sequence[u,v]}) && FreeQ[{u,v},Subscript]

Bra /: Bra[u__,{a___},{b___}] ** (f_Bra) ** Ket[v__,{a___},{b___}] := 
		Braket[{u,{a},{b}},{v,{a},{b}}] * f /; FreeQ[{u,v},Subscript]

Bra /: Bra[u__,{a___},{b___}] ** (f__Bra) ** Ket[v__,{a___},{b___}] := 
		Braket[{u,{a},{b}},{v,{a},{b}}] * NonCommutativeMultiply[f] /; FreeQ[{u,v},Subscript]

Bra /: Bra[u__] ** (f__Bra) ** (g__Ket) ** Ket[v__] := Braket[{u},{v}] * 
		(f ** g) /;  And@@(FreeQ[#,List]&/@{Sequence[u,v]}) && FreeQ[{u,v},Subscript]

Bra /: Bra[u__,{a___},{b___}] ** (f__Bra) ** (g__Ket) ** Ket[v__,{a___},{b___}] := 
		(f ** g) Braket[{u,{a},{b}},{v,{a},{b}}] /; FreeQ[{u,v},Subscript]



commonIndicesQ[bras_,kets_]:= Module[
	{indBras, indKets, takeInd},

	takeInd[x__,{y___},{z___}]:=If[Head[#]===Subscript,{#[[2]],y,z},{y,z}]&/@{x};
	takeInd[x__]:={#[[2]]}& /@ Cases[{x},Subscript[__]];

	indBras = Flatten[Map[takeInd@@#&,Map[List@@#&,bras]],1];	
	indKets = Flatten[Map[takeInd@@#&,Map[List@@#&,kets]],1];

	Intersection[indBras,indKets]!={}]



NonCommutativeMultiply[a__, TensorProduct[b__Bra], TensorProduct[c__Ket], d__]:=
	(NonCommutativeMultiply[a,b,c,d] /. TensorProduct -> NonCommutativeMultiply /. 
	NonCommutativeMultiply[x__Ket, y__Bra, u__Ket, v__Bra] :> NonCommutativeMultiply[x, u, y, v] /.
	NonCommutativeMultiply[x__Bra] :> If[Length[{x}]>1, TensorProduct[x], Unevaluated[x]] /. 
	NonCommutativeMultiply[y__Ket] :> If[Length[{y}]>1, TensorProduct[y], Unevaluated[y]]) 

NonCommutativeMultiply[a__, b_Bra, TensorProduct[c__Ket], d__]:=
	(NonCommutativeMultiply[a,b,c,d] /. TensorProduct -> NonCommutativeMultiply /. 
	NonCommutativeMultiply[x__Ket, y__Bra, u__Ket, v__Bra] :> NonCommutativeMultiply[x, u, y, v] /.
	NonCommutativeMultiply[x__Bra] :> If[Length[{x}]>1, TensorProduct[x], Unevaluated[x]] /. 
	NonCommutativeMultiply[y__Ket] :> If[Length[{y}]>1, TensorProduct[y], Unevaluated[y]]) 

NonCommutativeMultiply[a__, TensorProduct[b__Bra], c_Ket, d__]:=
	(NonCommutativeMultiply[a,b,c,d] /. TensorProduct -> NonCommutativeMultiply /. 
		NonCommutativeMultiply[x__Ket, y__Bra, u__Ket, v__Bra] :> NonCommutativeMultiply[x, u, y, v] /.
		NonCommutativeMultiply[x__Bra] :> If[Length[{x}]>1, TensorProduct[x], Unevaluated[x]] /. 
		NonCommutativeMultiply[y__Ket] :> If[Length[{y}]>1, TensorProduct[y], Unevaluated[y]]) 

NonCommutativeMultiply[a___, TensorProduct[b__Bra], TensorProduct[c__Ket], d___]:=
	(NonCommutativeMultiply[a,b,c,d] /. TensorProduct -> NonCommutativeMultiply /. 
	NonCommutativeMultiply[x__Ket, y__Bra, u__Ket, v__Bra] :> NonCommutativeMultiply[x, u, y, v] /.
	NonCommutativeMultiply[x__Bra] :> If[Length[{x}]>1, TensorProduct[x], Unevaluated[x]] /. 
	NonCommutativeMultiply[y__Ket] :> If[Length[{y}]>1, TensorProduct[y], Unevaluated[y]]) /; 
		commonIndicesQ[{b},{c}]

NonCommutativeMultiply[a___, b_Bra, TensorProduct[c__Ket], d___]:=
	(Evaluate[NonCommutativeMultiply[a, b, c, d]] /. TensorProduct -> NonCommutativeMultiply /. 
		NonCommutativeMultiply[x__Ket, y__Bra, u__Ket, v__Bra] :> NonCommutativeMultiply[x, u, y, v] /.
	NonCommutativeMultiply[x__Bra] :> If[Length[{x}]>1, TensorProduct[x], Unevaluated[x]] /. 
	NonCommutativeMultiply[y__Ket] :> If[Length[{y}]>1, TensorProduct[y], Unevaluated[y]]) /; 
		commonIndicesQ[{b},{c}]

NonCommutativeMultiply[a___, TensorProduct[b__Bra], c_Ket, d___]:=
	(Evaluate[NonCommutativeMultiply[a,b,c,d]] /. TensorProduct -> NonCommutativeMultiply /. 
		NonCommutativeMultiply[x__Ket, y__Bra, u__Ket, v__Bra] :> NonCommutativeMultiply[x, u, y, v] /.
		NonCommutativeMultiply[x__Bra] :> If[Length[{x}]>1, TensorProduct[x], Unevaluated[x]] /. 
		NonCommutativeMultiply[y__Ket] :> If[Length[{y}]>1, TensorProduct[y], Unevaluated[y]]) /; 
			commonIndicesQ[{b},{c}]

NonCommutativeMultiply[x:(Ket[__Subscript,{___},{___}]|Ket[__Subscript]),
						y:(Bra[__Subscript,{___},{___}]|Bra[__Subscript]),
						u:(Ket[__Subscript,{___},{___}]|Ket[__Subscript]),
						v:(Bra[__Subscript,{___},{___}]|Bra[__Subscript])] := 
								NonCommutativeMultiply[TensorProduct[x, u],TensorProduct[y, v]]

NonCommutativeMultiply[x:(Ket[__Subscript,{___},{___}]|Ket[__Subscript]),
						u:(Ket[__Subscript,{___},{___}]|Ket[__Subscript])] := TensorProduct[x, u]

NonCommutativeMultiply[y:(Bra[__Subscript,{___},{___}]|Bra[__Subscript]),
						v:(Bra[__Subscript,{___},{___}]|Bra[__Subscript])] := TensorProduct[y, v]


(* ::Subsection:: *)
(*Power of Operators*)


Power[Operator[1], n_ ] := Operator[1]

Power[Operator[0], n_ ] := 0

Power[Operator[0],0] := Indeterminate

Power[expr_,0] := Operator[1] /; OperatorQ[expr]

Power/: P_^m_. ** P_^n_. := Operator[1] /; OperatorQ[P] && n==-m


(* ::Subsection:: *)
(*Hermitian Operation*)


Hermitian[Operator[0]] := Operator[0]

Hermitian[Operator[1]] := Operator[1]

Hermitian[Operator[P_,x___]] := Operator[P,"H",x]

Hermitian[Operator[P_,"H",x___]] := Operator[P,x]

Hermitian[c_] := Conjugate[c] /; ConstantQ[c]

Hermitian[c_ P_] := Conjugate[c] Hermitian[P] /; ConstantQ[c] && !ConstantQ[P]

Hermitian[ P_ + Q_] := Hermitian[P] + Hermitian[Q] /; !ConstantQ[P] && !ConstantQ[Q]

Hermitian[ P_ ** Q_ ] := Hermitian[Q] ** Hermitian[P]

Hermitian[Ket[a__]] := Bra[a]

Hermitian[Bra[a__]] := Ket[a]

Hermitian[Operator[P__][t__]] := Hermitian[Operator[P]][t]

Hermitian[Power[P_,n_]] := Power[Hermitian[P],n] /; OperatorQ[P]

Hermitian[Power[E,P_]] := Power[E,Hermitian[P]] /; OperatorQ[P]

Hermitian[TensorProduct[a__]]:=TensorProduct[Sequence@@(Hermitian[#]&/@{a})]


(* ::Subsection:: *)
(*Commutation Rules*)


Commutator[P_ , 0] := 0 

 Commutator[0, P_] := 0 

 Commutator[P_ ,Operator[0]] := 0

 Commutator[Operator[0],P_] := 0

 Commutator[P_ ,Operator[1]] := 0 /; OperatorQ[P]

 Commutator[Operator[1],P_] := 0 /; OperatorQ[P]

 Commutator[P_, P_] := 0 /; OperatorQ[P]

 Commutator[P_, n_ Q_] := n Commutator[P,Q] /; 
		!OperatorQ[n] && OperatorQ[P] && OperatorQ[Q] 
 
 Commutator[m_ P_, Q_] := m Commutator[P,Q] /; 
		!OperatorQ[m] && OperatorQ[P] && OperatorQ[Q]

 Commutator[P_, Q_+ R_] := Commutator[P,Q] + Commutator[P,R] /; 
		OperatorQ[P] && OperatorQ[Q] && OperatorQ[R]

 Commutator[P_+ Q_, R_] := Commutator[P,R] + Commutator[Q,R] /; 
		OperatorQ[P] && OperatorQ[Q] && OperatorQ[R]

 Commutator[P_, HoldPattern[Q_** R_]] := Commutator[P,Q] ** R + Q ** Commutator[P,R]/; 
		OperatorQ[P] && OperatorQ[Q] && OperatorQ[R]

 Commutator[HoldPattern[P_** Q_], R_] := Commutator[P,R] ** Q + P ** Commutator[Q,R]/; 
		OperatorQ[P] && OperatorQ[Q] && OperatorQ[R]

 Hermitian[HoldPattern[Commutator[P_, Q_]]] := Commutator[Hermitian[Q],Hermitian[P]]/; 
						OperatorQ[P] && OperatorQ[Q]

 Commutator[X:(Operator[__]|Operator[__][__]), expr_?OperatorQ]:= 0 /; checkOperator[X,expr]

 Commutator[expr_?OperatorQ,X:(Operator[__]|Operator[__][t__])] := 0 /; checkOperator[X,expr]


(* ::Subsection:: *)
(*QAMatrix, QATrace, QADirac and QACompute*)


Options[QADirac] = {Orthonormal -> True};

QADirac[{bra_?VectorQ}, base_] := Plus @@ (bra * Hermitian[base])

QADirac[ket:{{_} ..}, base_] := Plus @@ (ket[[All, 1]] * base)

QADirac[op_?MatrixQ, base_, OptionsPattern[]]:=
  Module[{outerMatrix,innerMatrix,tempMatrix},
	outerMatrix = Outer[NonCommutativeMultiply, base, Hermitian[base]];
	If[OptionValue[Orthonormal],
		Plus @@ Flatten[ op * outerMatrix ],
		innerMatrix = Outer[NonCommutativeMultiply, Hermitian[base], base];
		tempMatrix =  op . Inverse[innerMatrix];
		Plus @@ Flatten[ tempMatrix * outerMatrix ]
     ]
  ]

QADirac[op_, base_, OptionsPattern[]]:=
  Module[{opMatrix,outerMatrix,innerMatrix,tempMatrix},
	opMatrix = Outer[NonCommutativeMultiply[#1,op,#2]&, Hermitian[base], base];
	outerMatrix = Outer[NonCommutativeMultiply, base, Hermitian[base]];
	If[OptionValue[Orthonormal],
		Plus @@ Flatten[ opMatrix * outerMatrix ],
		innerMatrix = Outer[NonCommutativeMultiply, Hermitian[base], base];
		tempMatrix = Inverse[innerMatrix] . opMatrix . Inverse[innerMatrix];
		Plus @@ Flatten[ tempMatrix * outerMatrix ]
     ]
  ]


Options[QAMatrix] = {Orthonormal -> True};

QAMatrix[ket_, base_] := {# ** ket} & /@ Hermitian[base] /; KetQ[ket] && FreeQ[ket, NonCommutativeMultiply] && FreeQ[ket, Operator]

QAMatrix[bra_, base_] := {bra ** # & /@ base} /; BraQ[bra] && FreeQ[bra, NonCommutativeMultiply] && FreeQ[ket, Operator]

QAMatrix[Power[op_, n_], base_] := MatrixPower[QAMatrix[op, base], n]

QAMatrix[op_, base_, OptionsPattern[]] := 
  Module[{innerMatrix, opMatrix},	
	If[OptionValue[Orthonormal],
		Outer[NonCommutativeMultiply[#1,op,#2]&, Hermitian[base], base],
		innerMatrix = Outer[NonCommutativeMultiply, Hermitian[base], base];
		opMatrix = Outer[NonCommutativeMultiply[#1,op,#2]&, Hermitian[base], base];
		Inverse[innerMatrix] . opMatrix ] 
  ]


QATrace[op_, base_] := Sum[Hermitian[(base[[i]])] ** op ** base[[i]],
	{i, 1, Length[base]}]

SetAttributes[QACompute, HoldAll];
Options[QACompute] = {"ResultForm" -> "Dirac"};
	
QACompute[Power[op_, n_], base_, OptionsPattern[]] := 
	Switch[OptionValue["ResultForm"],
		"Dirac", QADirac[QAMatrix[Power[op, n], base], base],
		"Matrix", QAMatrix[Power[op, n], base]
		]

QACompute[CenterDot[factors__], base_, OptionsPattern[]] := 
	Module[{newfactors},
		newfactors = {factors} /. Power[x_, y_] /; OperatorQ[x] :> Sequence @@ Table[x, {y}];
		Switch[OptionValue["ResultForm"],
			"Dirac", QADirac[Dot@@(QAMatrix[#,base]&/@ newfactors), base],
			"Matrix", Dot @@ (QAMatrix[#, base] & /@ newfactors)
			]
		]

QACompute[ket_, base_, OptionsPattern[]] := 
	Switch[OptionValue["ResultForm"],
		"Dirac", QADirac[QAMatrix[ket, base], base],
		"Matrix", QAMatrix[ket, base]
		] /; KetQ[ket]

QACompute[bra_, base_, OptionsPattern[]] :=
	Switch[OptionValue["ResultForm"],
		"Dirac", QADirac[First@QAMatrix[bra, base], base],
		"Matrix", QAMatrix[bra, base]
		] /; BraQ[bra]

QACompute[op_, base_, OptionsPattern[]] := 
	Switch[OptionValue["ResultForm"],
		"Dirac",  QADirac[QAMatrix[op, base], base],
		"Matrix", QAMatrix[op, base]
		] /; OperatorQ[op]
 
 


(* ::Subsection:: *)
(*Gram Schmidt Procedure*)


GramSchmidtVectors[l_] := Module[{wk, t1, braketRule, conjRules}, 
	braketRule = HoldPattern[Bra[x__] ** Ket[y__]] :> Braket[{x},{y}];
	conjRules = {Conjugate[1/Sqrt[x_]] -> 1/Sqrt[Conjugate[x]], 
				 Conjugate[Sqrt[x_]] -> Sqrt[Conjugate[x]]}; 
    wk[1] = l[[1]]/Sqrt[Hermitian[l[[1]]]**l[[1]]]; 
    wk[n_] := wk[n] = (t1 = l[[n]] - Sum[Hermitian[wk[i]]**l[[n]]*wk[i], {i, 1, n - 1}] /. braketRule; 
       t1/Sqrt[Hermitian[t1]**t1] /. braketRule //. conjRules); Table[wk[j], {j, 1, Length[l]}]]



(* ::Subsection:: *)
(*QASet AND QAClear*)


SetAttributes[ QASet, {HoldAllComplete} ];
SetAttributes[ QAClear, {HoldAllComplete} ];
Options[QASet] = {ExtendSystem -> False};


divideKetOrBra[ketorbra_] := Flatten[ ketorbra /. TensorProduct -> List /. 
	{   Ket[x__,{y___},{z___}] :> Map[Ket[#,{y},{z}]&,{x}],
		Bra[x__,{y___},{z___}] :> Map[Bra[#,{y},{z}]&,{x}], 
		Ket[x__] :> Map[Ket[#]&,{x}],
		Bra[x__] :> Map[Bra[#]&,{x}]
}];

listOfPatterns[expr_] := Module[{p}, Cases[expr //. Pattern -> p, _p, \[Infinity]] //. p -> Pattern ];

deleteBlanks[expr_] := expr /. Verbatim[Pattern][t_, _] :> t;

filterOperators[expr_] := First[Union[Cases[expr,Operator[__],\[Infinity]]]];


QASet[expr__, options:OptionsPattern[]] := 
	Module[{setRules, invSetRules, set, setDelayed}, 
		setRules = {Set -> set, SetDelayed -> setDelayed};
		invSetRules = { set -> Set, setDelayed -> SetDelayed};
		QASet[#,options]&/@ ReleaseHold[Hold[{expr}] /. setRules] /. 
		invSetRules;
	] /; (Hold[{expr}] /. Hold[e_] :> Length[Unevaluated[e]]>=2) && FreeQ[Hold[{expr}],Rule]

QAClear[expr__] :=  Module[{}, 
    ReleaseHold[Map[QAClear,Hold[{expr}],{2}]]; ] /; (Hold[{expr}] /. 
		Hold[e_] :> Length[Unevaluated[e]]>=2)


applyQ[op_,ketorbra1_,ketorbra2_] := 
	Module[ {lista1, lista2, inter, 
				newinter, rules, afected},
 
		inter = Intersection[listOfPatterns[op],listOfPatterns[ketorbra1]];
		newinter = deleteBlanks[inter];

		rules = Thread[(Verbatim/@inter) -> newinter];
		
		lista1 = divideKetOrBra[ketorbra1 /. rules];
		lista2 = divideKetOrBra[ketorbra2];		

		afected = Flatten[Cases[lista2, #]& /@ lista1];

		(afected != {}) && Length[afected] == Length[lista1]
	]


applyQ[op_, x_, ketorbra1_, ketorbra2_] := Module[ {lista1, lista2, inter, 
				newinter, rules, afected, optemp, rules2},

		optemp = First[Union[Cases[{x},Operator[__],\[Infinity]]]];
		inter = Intersection[listOfPatterns[op],listOfPatterns[ketorbra1]];
		newinter = deleteBlanks[inter];

		rules = Thread[(Verbatim/@inter) -> newinter];
		rules2 = Thread[deleteBlanks/@Cases[Verbatim[op], Verbatim[Pattern][_,_],\[Infinity]] -> 
						Extract[optemp, Position[op,Verbatim[Pattern][_,_]]]];

		lista1 = divideKetOrBra[(ketorbra1 /. rules) /. rules2];
		lista2 = divideKetOrBra[ketorbra2];	
		
		afected = Flatten[Cases[lista2, #]& /@ lista1];

		(afected != {}) && Length[afected] == Length[lista1]
	]


QASet[Operator/: op_ ** ket1_ := ket2_, OptionsPattern[]]:=
	Module[{const,ket4},
	
	Off[RuleDelayed::rhs]; 
	Off[Condition::condp];
	Unprotect[NonCommutativeMultiply];

	Operator/: op ** ket1 := ket2;
	Operator/: Hermitian[ket1] ** Hermitian[op] := Hermitian[ket2];

	const = ket2 /. Ket[__]->1;
	ket4 = First[Cases[{ket2}, Ket[__], \[Infinity]]];
	
	If[ deleteBlanks[ket1] === ket4,
	
		Ket/: x_ ** ket1 /; (MatchQ[Union[Cases[{x},Operator[__],\[Infinity]]],{op} /. 
				Verbatim[Pattern][t_,z_]:>Pattern[Evaluate[ToExpression["$$$"~~ToString[t]~~"$$$"]],z]] && 
				Head[x] =!= NonCommutativeMultiply && Head[x] =!= Operator && applyQ[op, x, ket1, deleteBlanks[ket1]] ) := 
					Module[{c,op2},
						c = Evaluate[deleteBlanks[op]**deleteBlanks[ket1]] /. Ket[__] -> 1;
						op2 = op /. Verbatim[Pattern][t_,z_]:>Pattern[Evaluate[ToExpression["$$$"~~ToString[t]~~"$$$"]],z];
						QAProductExpand[( x deleteBlanks[ket1] /. op2 :> c )]
					];

		Bra/: Hermitian[ket1] ** x_ /; (MatchQ[Union[Cases[{x},Operator[__],\[Infinity]]],{Hermitian[op]}/. 
				Verbatim[Pattern][t_,z_]:>Pattern[Evaluate[ToExpression["$$$"~~ToString[t]~~"$$$"]],z]] && 
				Head[x] =!= NonCommutativeMultiply && Head[x] =!= Operator && applyQ[Hermitian[op], x, ket1, deleteBlanks[ket1]]) := 
					Module[{c,op2},
						c = deleteBlanks[op] ** deleteBlanks[ket1] /. Ket[__] -> 1;
						op2 = op /. Verbatim[Pattern][t_,z_]:>Pattern[Evaluate[ToExpression["$$$"~~ToString[t]~~"$$$"]],z];
						QAProductExpand[( x Hermitian[deleteBlanks[ket1]] /. Hermitian[op2] :> Conjugate[c] )]
					];
	];
	
	If[OptionValue[ExtendSystem],
		Operator/: op ** ket3_ /; applyQ[op, ket1, ket3] && KetQ[ket3] := 
			Module[{newket1, lista1, lista2, inter, 
						newinter, newop, rules, afected},
				
				inter = Intersection[listOfPatterns[op], listOfPatterns[ket1]];
				newinter = deleteBlanks[inter];
				rules = Thread[(Verbatim/@inter)->newinter];
				newop = deleteBlanks[op];
				lista1 = divideKetOrBra[ket1/.rules];
				lista2 = divideKetOrBra[ket3];
				afected = Flatten[Cases[lista2,#]&/@lista1];
				newket1 = If[Length[afected]==1, Sequence@@afected, TensorProduct@@afected];
				(newop ** newket1) \[CircleTimes] (TensorProduct@@Complement[lista2,afected])
			]  ;

		Operator/: bra_ ** Hermitian[op] /; applyQ[op, Hermitian[ket1], bra] && BraQ[bra] := 
			Module[{newbra, lista1, lista2, inter, 
						newinter, newop, rules, afected},
				
				inter = Intersection[listOfPatterns[op], listOfPatterns[ket1]];
				newinter = deleteBlanks[inter];
				rules = Thread[(Verbatim/@inter)->newinter];
				newop = deleteBlanks[Hermitian[op]];
				lista1 = divideKetOrBra[Hermitian[ket1]/.rules];
				lista2 = divideKetOrBra[bra];
				afected = Flatten[Cases[lista2,#]&/@lista1];
				newbra = If[Length[afected]==1, Sequence@@afected, TensorProduct@@afected];
				(TensorProduct@@Complement[lista2,afected]) \[CircleTimes] (newbra ** newop)
			]; 

		TensorProduct/: TensorProduct[op, y_] ** ket3_ /; 
							applyQ[op,ket1,ket3]&&KetQ[ket3] := 
									y ** (deleteBlanks[op] ** ket3) ;
	
		TensorProduct/: bra_ ** TensorProduct[Hermitian[op], y_] /; 
							applyQ[op,Hermitian[ket1],bra] && BraQ[bra] := 
									(bra ** deleteBlanks[Hermitian[op]]) ** y ;

		If[ deleteBlanks[ket1] === ket4,
	
				x_ ** ket3_ /; (MatchQ[Union[Cases[{x},Operator[__],\[Infinity]]], {op}] && 
					Head[x] =!= NonCommutativeMultiply && Head[x] =!= Operator &&
					applyQ[op, x, ket1, ket3] && KetQ[ket3]) := 
						Module[{newket1, lista1, lista2, inter, 
									newinter, newop, rules, afected, rules2, optemp},
							optemp = First[Union[Cases[{x},Operator[__],\[Infinity]]]];
							inter = Intersection[listOfPatterns[op], listOfPatterns[ket1]];
							newinter = deleteBlanks[inter];
							rules = Thread[(Verbatim/@inter)->newinter];
							rules2 = Thread[deleteBlanks/@Cases[Verbatim[op], Verbatim[Pattern][_,_],\[Infinity]] -> 
										Extract[optemp, Position[op,Verbatim[Pattern][_,_]]]];
							newop = deleteBlanks[op];
							lista1 = divideKetOrBra[(ket1/.rules)/.rules2];
							lista2 = divideKetOrBra[ket3];
							afected = Flatten[Cases[lista2,#]&/@lista1];
							newket1 = If[Length[afected]==1, Sequence@@afected, TensorProduct@@afected];
							(x ** newket1) \[CircleTimes] (TensorProduct@@Complement[lista2,afected])
						];

				TensorProduct/: TensorProduct[Y_, x_] ** ket3_ /; 
									(MatchQ[Union[Cases[{Y},Operator[__],\[Infinity]]], {op}] &&
									applyQ[op, First[Union[Cases[{Y},Operator[__],\[Infinity]]]], ket1, ket3] && KetQ[ket3]) := x ** (Y ** ket3);


				bra_ ** x_ /; (MatchQ[Union[Cases[{x},Operator[__],\[Infinity]]], {Hermitian[op]}] && 
					Head[x] =!= NonCommutativeMultiply && Head[x] =!= Operator && 
					applyQ[Hermitian[op], x, Hermitian[ket1], bra] && BraQ[bra]) := 
						Module[{newbra, lista1, lista2, inter, 
									newinter, newop, rules, afected, rules2, optemp},
							optemp = First[Union[Cases[{x},Operator[__],\[Infinity]]]];
							inter = Intersection[listOfPatterns[op], listOfPatterns[ket1]];
							newinter = deleteBlanks[inter];
							rules = Thread[(Verbatim /@ inter)->newinter];
							rules2 = Thread[deleteBlanks /@ Cases[Verbatim[op], Verbatim[Pattern][_,_],\[Infinity]] -> 
										Extract[optemp, Position[Hermitian[op],Verbatim[Pattern][_,_]]]];
							newop = deleteBlanks[Hermitian[op]];
							lista1 = divideKetOrBra[(Hermitian[ket1]/.rules)/.rules2];
							lista2 = divideKetOrBra[bra];
							afected = Flatten[Cases[lista2,#]&/@lista1];
							newbra = If[Length[afected]==1, Sequence@@afected, TensorProduct@@afected];
							(TensorProduct@@Complement[lista2,afected]) \[CircleTimes] (newbra ** x)
						];

				TensorProduct/: bra_ ** TensorProduct[Y_, x_] /; 
									(MatchQ[Union[Cases[{Y},Operator[__],\[Infinity]]], {Hermitian[op]}] &&
									applyQ[Hermitian[op], First[Union[Cases[{Y},Operator[__],\[Infinity]]]], Hermitian[ket1], bra] 
									&& BraQ[bra]) := (bra ** Y) ** x;

		];

	];
	
	Protect[NonCommutativeMultiply];
	On[RuleDelayed::rhs]; 
	On[Condition::condp];

	] /; OperatorQ[op] && KetQ[ket1] && KetQ[ket2]


QAClear[Operator/: op_ ** ket1_ =. ] :=
	Module[{},

	Unprotect[NonCommutativeMultiply];	
	Off[RuleDelayed::rhs]; 
	Off[Condition::condp];
	Off[TagUnset::norep];
	Off[Unset::norep];

	Operator/: op ** ket1 =. ;
	Operator/: Hermitian[ket1] ** Hermitian[op] =. ;
	
	Ket/: x_ ** ket1 /; (MatchQ[Union[Cases[{x},Operator[__],\[Infinity]]],{op} /. 
		Verbatim[Pattern][t_,z_]:>Pattern[Evaluate[ToExpression["$$$"~~ToString[t]~~"$$$"]],z]] && 
				Head[x] =!= NonCommutativeMultiply && Head[x] =!= Operator) =. ;

	Bra/: Hermitian[ket1] ** x_ /; (MatchQ[Union[Cases[{x},Operator[__],\[Infinity]]],{Hermitian[op]}/. 
		Verbatim[Pattern][t_,z_]:>Pattern[Evaluate[ToExpression["$$$"~~ToString[t]~~"$$$"]],z]] && 
				Head[x] =!= NonCommutativeMultiply && Head[x] =!= Operator) =. ;
	
	Operator/: op ** ket3_ /; applyQ[op, ket1, ket3] && KetQ[ket3] =. ;

	Operator/: bra_ ** Hermitian[op] /; applyQ[op, Hermitian[ket1], bra] && BraQ[bra] =.;

	TensorProduct/: TensorProduct[op, y_] ** ket3_ /; 
						applyQ[op,ket1,ket3]&&KetQ[ket3] =. ;
	
	TensorProduct/: bra_ ** TensorProduct[Hermitian[op], y_] /; 
						applyQ[op,Hermitian[ket1],bra] && BraQ[bra] =. ;
	
	x_ ** ket3_ /; (MatchQ[Union[Cases[{x},Operator[__],\[Infinity]]], {op}] && 
			Head[x] =!= NonCommutativeMultiply && Head[x] =!= Operator &&
					applyQ[op, x, ket1, ket3] && KetQ[ket3]) =. ;

	TensorProduct/: TensorProduct[Y_, x_] ** ket3_ /; 
			(MatchQ[Union[Cases[{Y},Operator[__],\[Infinity]]], {op}] &&
					applyQ[op, First[Union[Cases[{Y},Operator[__],\[Infinity]]]], ket1, ket3] && KetQ[ket3]) =. ;

	bra_ ** x_ /; (MatchQ[Union[Cases[{x},Operator[__],\[Infinity]]], {Hermitian[op]}] && 
			Head[x] =!= NonCommutativeMultiply && Head[x] =!= Operator && 
					applyQ[Hermitian[op], x, Hermitian[ket1], bra] && BraQ[bra]) =. ;

	TensorProduct/: bra_ ** TensorProduct[Y_, x_] /; 
			(MatchQ[Union[Cases[{Y},Operator[__],\[Infinity]]], {Hermitian[op]}] &&
					applyQ[Hermitian[op], First[Union[Cases[{Y},Operator[__],\[Infinity]]]], Hermitian[ket1], bra] && BraQ[bra]) =. ;
	
	Protect[NonCommutativeMultiply];	
	On[RuleDelayed::rhs]; 
	On[Condition::condp];
	On[TagUnset::norep];
	On[Unset::norep];

	] /; OperatorQ[op] && KetQ[ket1] 


QASet[ Set[ Commutator[Commutator[P_,Q_],Commutator[R_,S_]], T_ ], OptionsPattern[]] :=
		Module[ { },
                  Commutator[Commutator[P,Q],Commutator[R,S]] = T;
                  Commutator[Commutator[P,Q],Commutator[S,R]] = -T;
                  Commutator[Commutator[Q,P],Commutator[R,S]] = -T;
                  Commutator[Commutator[Q,P],Commutator[S,R]] = T;
                  Commutator[Commutator[R,S],Commutator[P,Q]] = -T;
                  Commutator[Commutator[R,S],Commutator[Q,P]] = T;
                  Commutator[Commutator[S,R],Commutator[P,Q]] = T;
                  Commutator[Commutator[S,R],Commutator[Q,P]] = -T; ]

QASet[ SetDelayed[ Commutator[Commutator[P_,Q_],Commutator[R_,S_]], T_ ], OptionsPattern[]] :=
		Module[ { },
                  Commutator[Commutator[P,Q],Commutator[R,S]] := T;
                  Commutator[Commutator[P,Q],Commutator[S,R]] := -T;
                  Commutator[Commutator[Q,P],Commutator[R,S]] := -T;
                  Commutator[Commutator[Q,P],Commutator[S,R]] := T;
                  Commutator[Commutator[R,S],Commutator[P,Q]] := -T;
                  Commutator[Commutator[R,S],Commutator[Q,P]] := T;
                  Commutator[Commutator[S,R],Commutator[P,Q]] := T;
                  Commutator[Commutator[S,R],Commutator[Q,P]] := -T; ]

QAClear[ Unset[ Commutator[Commutator[P_,Q_],Commutator[R_,S_]] ] ] :=
		Module[ { },
                  Commutator[Commutator[P,Q],Commutator[R,S]] =.;
                  Commutator[Commutator[P,Q],Commutator[S,R]] =.;
                  Commutator[Commutator[Q,P],Commutator[R,S]] =.;
                  Commutator[Commutator[Q,P],Commutator[S,R]] =.;
                  Commutator[Commutator[R,S],Commutator[P,Q]] =.;
                  Commutator[Commutator[R,S],Commutator[Q,P]] =.;
                  Commutator[Commutator[S,R],Commutator[P,Q]] =.;
                  Commutator[Commutator[S,R],Commutator[Q,P]] =.; ]


QASet[ Set[ ( Commutator[P_, Commutator[Q_,R_]] |
			  Commutator[Commutator[R_,Q_], P_] ), S_ ], 
		OptionsPattern[]] :=
	Module[ { }, 
		Commutator[P,Commutator[Q,R]] = S;
		Commutator[P,Commutator[R,Q]] = -S;
		Commutator[Commutator[Q,R],P] = -S;
		Commutator[Commutator[R,Q],P] = S; ]

QASet[ SetDelayed[ ( Commutator[P_, Commutator[Q_,R_]] |
					 Commutator[Commutator[R_,Q_], P_] ), S_ ],
		OptionsPattern[] ] :=
	Module[ { }, 
		Commutator[P,Commutator[Q,R]] := S;
		Commutator[P,Commutator[R,Q]] := -S;
		Commutator[Commutator[Q,R],P] := -S;
		Commutator[Commutator[R,Q],P] := S; ]

QAClear[ Unset[ ( Commutator[P_, Commutator[Q_,R_]] |
          Commutator[Commutator[R_,Q_], P_] )] ] :=
                Module[ { }, 
		    Commutator[P,Commutator[Q,R]] =.;
                  Commutator[P,Commutator[R,Q]] =.;
                  Commutator[Commutator[Q,R],P] =.;
                  Commutator[Commutator[R,Q],P] =.; ]


QASet[ Set[ Commutator[P_,Q_], R_] , OptionsPattern[]] := 
	Module[ {},

		Off[RuleDelayed::rhs]; 

		Commutator[P,Q] = R;
		Commutator[Q,P] = -R;           
        
		If[ Commutator[deleteBlanks[P], R] === Commutator[deleteBlanks[Q], R] === 0,
			
			Commutator[P^n_.,Q^m_.] /; IntegerQ[n] && IntegerQ[m] :=
			 Sum[Binomial[n,k] (m!/(m-k)!) (R^k) **
		              deleteBlanks[Q]^(m-k) ** deleteBlanks[P]^(n-k),{k,1,n}];

			Commutator[Q^n_.,P^m_.] /; IntegerQ[n] && IntegerQ[m] :=
			 Sum[Binomial[n,k] (m!/(m-k)!) ((-1)^k)(R^k) **
				deleteBlanks[P]^(m-k) ** deleteBlanks[Q]^(n-k),{k,1,n}];

			Commutator[P,x_] /; MatchQ[Union[Cases[{x},Operator[__],\[Infinity]]],{Q}] := 
				Module[
					{r = filterOperators[x], d = D[x,filterOperators[x]]/.Times->NonCommutativeMultiply},
					Commutator[deleteBlanks[P], r]  ** ( d ) // 
					MapAll[QAProductExpand,#]&] ;

			Commutator[x_,P] /; MatchQ[Union[Cases[{x},Operator[__],\[Infinity]]],{Q}] := 
				Module[
					{r = filterOperators[x], d = D[x,filterOperators[x]]/.Times->NonCommutativeMultiply},
					- Commutator[deleteBlanks[P], r]  ** ( d ) // 
					MapAll[QAProductExpand,#]&] ;

			Commutator[Q, x_] /; MatchQ[Union[Cases[{x},Operator[__],\[Infinity]]],{P}] := 
				Module[
					{r = filterOperators[x], d = D[x,filterOperators[x]]/.Times->NonCommutativeMultiply},
					- Commutator[r, deleteBlanks[Q]]  ** ( d ) // 
					MapAll[QAProductExpand,#]&];

			Commutator[x_,Q] /; MatchQ[Union[Cases[{x},Operator[__],\[Infinity]]],{P}] := 
				Module[
					{r = filterOperators[x], d = D[x,filterOperators[x]]/.Times->NonCommutativeMultiply},
					Commutator[r, deleteBlanks[Q]]  ** ( d ) // 
					MapAll[QAProductExpand,#]&] ;

			On[RuleDelayed::rhs]; 
		]; 
	] /; R =!= 0

QASet[ SetDelayed[ Commutator[P_,Q_], R_], OptionsPattern[]] := 
	Module[ {},
		Off[RuleDelayed::rhs]; 

		Commutator[P,Q] := R;
		Commutator[Q,P] := -R;           
        
		If[ Commutator[deleteBlanks[P], R] === Commutator[deleteBlanks[Q], R] === 0,
			
			Commutator[P^n_.,Q^m_.] /; IntegerQ[n] && IntegerQ[m] :=
			 Sum[Binomial[n,k] (m!/(m-k)!) (R^k) **
		              deleteBlanks[Q]^(m-k) ** deleteBlanks[P]^(n-k),{k,1,n}];

			Commutator[Q^n_.,P^m_.] /; IntegerQ[n] && IntegerQ[m] :=
			 Sum[Binomial[n,k] (m!/(m-k)!) ((-1)^k)(R^k) **
				deleteBlanks[P]^(m-k) ** deleteBlanks[Q]^(n-k),{k,1,n}];

			Commutator[P,x_] /; MatchQ[Union[Cases[{x},Operator[__],\[Infinity]]],{Q}] := Module[
				{r = filterOperators[x], d = D[x,filterOperators[x]]/.Times->NonCommutativeMultiply},
				Commutator[deleteBlanks[P], r]  ** ( d ) // 
				MapAll[QAProductExpand,#]&] ;

			Commutator[x_,P] /; MatchQ[Union[Cases[{x},Operator[__],\[Infinity]]],{Q}] := Module[
				{r = filterOperators[x], d = D[x,filterOperators[x]]/.Times->NonCommutativeMultiply},
				- Commutator[deleteBlanks[P], r]  ** ( d ) // 
				MapAll[QAProductExpand,#]&] ;

			Commutator[Q, x_] /; MatchQ[Union[Cases[{x},Operator[__],\[Infinity]]],{P}] := Module[
				{r = filterOperators[x], d = D[x,filterOperators[x]]/.Times->NonCommutativeMultiply},
				- Commutator[r, deleteBlanks[Q]]  ** ( d ) // 
				MapAll[QAProductExpand,#]&];

			Commutator[x_,Q] /; MatchQ[Union[Cases[{x},Operator[__],\[Infinity]]],{P}] := Module[
				{r = filterOperators[x], d = D[x,filterOperators[x]]/.Times->NonCommutativeMultiply},
				Commutator[r, deleteBlanks[Q]]  ** ( d ) // 
				MapAll[QAProductExpand,#]&] ;
			
			On[RuleDelayed::rhs]; 
		]; 
	] /; R =!= 0


QASet[ Set[ Commutator[P_,Q_], 0 ], OptionsPattern[]] :=
	Module[ {}, 
		Off[RuleDelayed::rhs];      
		Off[Condition::condp];

		Commutator[P,Q] = 0;
		Commutator[Q,P] = 0;
		
		Commutator[P,x_] /; MatchQ[Join[{(P /. Verbatim[Pattern][t_, _] :> t)},Cases[{x},Operator[__],\[Infinity]]],
			{P,Q}/.Verbatim[Pattern][t_,z_]:>Pattern[Evaluate[ToExpression["$$$"~~ToString[t]~~"$$$"]],z]] := 0  ;
	    
		Commutator[x_,P] /; MatchQ[Join[{(P /. Verbatim[Pattern][t_, _] :> t)},Cases[{x},Operator[__],\[Infinity]]],
			{P,Q}/.Verbatim[Pattern][t_,z_]:>Pattern[Evaluate[ToExpression["$$$"~~ToString[t]~~"$$$"]],z]]:= 0  ;
		
		Commutator[Q,x_] /; MatchQ[Join[{(Q /. Verbatim[Pattern][t_, _] :> t)},Cases[{x},Operator[__],\[Infinity]]],
			{Q,P}/.Verbatim[Pattern][t_,z_]:>Pattern[Evaluate[ToExpression["$$$"~~ToString[t]~~"$$$"]],z]] := 0  ;
		
		Commutator[x_,Q]  /; MatchQ[Join[{(Q /. Verbatim[Pattern][t_, _] :> t)},Cases[{x},Operator[__],\[Infinity]]],
		{Q,P}/.Verbatim[Pattern][t_,z_]:>Pattern[Evaluate[ToExpression["$$$"~~ToString[t]~~"$$$"]],z]] := 0  ; 
		
		Commutator[x_,y_]  /; MatchQ[Join[Cases[{x},Operator[__],\[Infinity]],Cases[{y},Operator[__],\[Infinity]]],
		{P,Q}/.Verbatim[Pattern][t_,z_]:>Pattern[Evaluate[ToExpression["$$$"~~ToString[t]~~"$$$"]],z]] := 0  ;
		
		Commutator[x_,y_] /; MatchQ[Join[Cases[{y},Operator[__],\[Infinity]],Cases[{x},Operator[__],\[Infinity]]],
		{P,Q}/.Verbatim[Pattern][t_,z_]:>Pattern[Evaluate[ToExpression["$$$"~~ToString[t]~~"$$$"]],z]] := 0  ; 
		
		On[RuleDelayed::rhs]; 
		On[Condition::condp];];

QASet[ SetDelayed[ Commutator[P_,Q_], 0 ], OptionsPattern[]] :=
	Module[ {}, 
		Off[RuleDelayed::rhs];      
		Off[Condition::condp];

		Commutator[P,Q] := 0;
		Commutator[Q,P] := 0;
		
		Commutator[P,x_] /; MatchQ[Join[{(P /. Verbatim[Pattern][t_, _] :> t)},Cases[{x},Operator[__],\[Infinity]]],
			{P,Q}/.Verbatim[Pattern][t_,z_]:>Pattern[Evaluate[ToExpression["$$$"~~ToString[t]~~"$$$"]],z]] := 0  ;
	    
		Commutator[x_,P] /; MatchQ[Join[{(P /. Verbatim[Pattern][t_, _] :> t)},Cases[{x},Operator[__],\[Infinity]]],
			{P,Q}/.Verbatim[Pattern][t_,z_]:>Pattern[Evaluate[ToExpression["$$$"~~ToString[t]~~"$$$"]],z]]:= 0  ;
		
		Commutator[Q,x_] /; MatchQ[Join[{(Q /. Verbatim[Pattern][t_, _] :> t)},Cases[{x},Operator[__],\[Infinity]]],
			{Q,P}/.Verbatim[Pattern][t_,z_]:>Pattern[Evaluate[ToExpression["$$$"~~ToString[t]~~"$$$"]],z]] := 0  ;
		
		Commutator[x_,Q]  /; MatchQ[Join[{(Q /. Verbatim[Pattern][t_, _] :> t)},Cases[{x},Operator[__],\[Infinity]]],
		{Q,P}/.Verbatim[Pattern][t_,z_]:>Pattern[Evaluate[ToExpression["$$$"~~ToString[t]~~"$$$"]],z]] := 0  ; 
		
		Commutator[x_,y_]  /; MatchQ[Join[Cases[{x},Operator[__],\[Infinity]],Cases[{y},Operator[__],\[Infinity]]],
		{P,Q}/.Verbatim[Pattern][t_,z_]:>Pattern[Evaluate[ToExpression["$$$"~~ToString[t]~~"$$$"]],z]] := 0  ;
		
		Commutator[x_,y_] /; MatchQ[Join[Cases[{y},Operator[__],\[Infinity]],Cases[{x},Operator[__],\[Infinity]]],
		{P,Q}/.Verbatim[Pattern][t_,z_]:>Pattern[Evaluate[ToExpression["$$$"~~ToString[t]~~"$$$"]],z]] := 0  ; 
		
		On[RuleDelayed::rhs]; 
		On[Condition::condp];];


QACommute[op_List] := Module[{}, 
		QASet[Commutator[#1,#2]=0]& @@@ Subsets[op,{2}];]

QAClearCommute[op_List] := Module[{}, 
		QAClear[Commutator[#1,#2]=.]& @@@ Subsets[op,{2}];]


QAClear[ Unset[ Commutator[P_,Q_] ] ] := 
	Module[ {   T = Commutator[deleteBlanks[P],deleteBlanks[Q]],
				R = Commutator[deleteBlanks[P],Commutator[deleteBlanks[P],deleteBlanks[Q]]], 
				S = Commutator[deleteBlanks[Q],Commutator[deleteBlanks[P],deleteBlanks[Q]]]},
		
		Off[Unset::norep];
		Off[Condition::condp];
		
		Commutator[P,Q] =. ;
		Commutator[Q,P] =. ;           
        
		If[ R===S===0 && T=!=0,
			
			Commutator[P^n_.,Q^m_.] /; IntegerQ[n] && IntegerQ[m] =.;
			Commutator[Q^n_.,P^m_.] /; IntegerQ[n] && IntegerQ[m] =.;
			Commutator[P,x_] /; MatchQ[Union[Cases[{x},Operator[__],\[Infinity]]],{Q}] =.;
			Commutator[x_,P] /; MatchQ[Union[Cases[{x},Operator[__],\[Infinity]]],{Q}] =.;
			Commutator[Q,x_] /; MatchQ[Union[Cases[{x},Operator[__],\[Infinity]]],{P}] =.;
			Commutator[x_,Q] /; MatchQ[Union[Cases[{x},Operator[__],\[Infinity]]],{P}] =.; ];

		If[T===0,

		Commutator[P,x_] /; MatchQ[Join[{(P /. Verbatim[Pattern][t_, _] :> t)},Cases[{x},Operator[__],\[Infinity]]],
			{P,Q}/.Verbatim[Pattern][t_,z_]:>Pattern[Evaluate[ToExpression["$$$"~~ToString[t]~~"$$$"]],z]] =. ;
	    
		Commutator[x_,P] /; MatchQ[Join[{(P /. Verbatim[Pattern][t_, _] :> t)},Cases[{x},Operator[__],\[Infinity]]],
			{P,Q}/.Verbatim[Pattern][t_,z_]:>Pattern[Evaluate[ToExpression["$$$"~~ToString[t]~~"$$$"]],z]] =.  ;
		
		Commutator[Q,x_] /; MatchQ[Join[{(Q /. Verbatim[Pattern][t_, _] :> t)},Cases[{x},Operator[__],\[Infinity]]],
			{Q,P}/.Verbatim[Pattern][t_,z_]:>Pattern[Evaluate[ToExpression["$$$"~~ToString[t]~~"$$$"]],z]] =.  ;
		
		Commutator[x_,Q]  /; MatchQ[Join[{(Q /. Verbatim[Pattern][t_, _] :> t)},Cases[{x},Operator[__],\[Infinity]]],
		{Q,P}/.Verbatim[Pattern][t_,z_]:>Pattern[Evaluate[ToExpression["$$$"~~ToString[t]~~"$$$"]],z]] =.  ; 
		
		Commutator[x_,y_]  /; MatchQ[Join[Cases[{x},Operator[__],\[Infinity]],Cases[{y},Operator[__],\[Infinity]]],
		{P,Q}/.Verbatim[Pattern][t_,z_]:>Pattern[Evaluate[ToExpression["$$$"~~ToString[t]~~"$$$"]],z]] =.  ;
		
		Commutator[x_,y_] /; MatchQ[Join[Cases[{y},Operator[__],\[Infinity]],Cases[{x},Operator[__],\[Infinity]]],
		{P,Q}/.Verbatim[Pattern][t_,z_]:>Pattern[Evaluate[ToExpression["$$$"~~ToString[t]~~"$$$"]],z]] =.  ; ]; 
	
		On[Unset::norep];
		On[Condition::condp];
	] 


QASet[{ Set[ ( Commutator[P_,Commutator[P_,Q_]] | Commutator[P_,Commutator[Q_,P_]] |
               Commutator[Commutator[P_,Q_],P_] | Commutator[Commutator[Q_,P_],P_]  ), 0 ],
		Set[ ( Commutator[Q_,Commutator[P_,Q_]] | Commutator[Q_,Commutator[Q_,P_]] |
               Commutator[Commutator[P_,Q_],Q_] | Commutator[Commutator[Q_,P_],Q_]  ), 0 ]},
		OptionsPattern[] ] :=

	Module[ {k},	

			Commutator[P,Commutator[P,Q]] = 0;
			Commutator[P,Commutator[Q,P]] = 0;
			Commutator[Commutator[P,Q],P] = 0;
			Commutator[Commutator[Q,P],P] = 0;
			Commutator[Q,Commutator[P,Q]] = 0;
			Commutator[Q,Commutator[Q,P]] = 0;
			Commutator[Commutator[P,Q],Q] = 0;
			Commutator[Commutator[Q,P],Q] = 0;
		   
			Commutator[P^n_,Q^m_] /; IntegerQ[n] && IntegerQ[m] := 
				Sum[Binomial[n,k] (m!/(m-k)!) (Commutator[P,Q]^k) ** Q^(m-k)**P^(n-k),{k,1,n}];
		   
			Commutator[Q^n_,P^m_] /; IntegerQ[n] && IntegerQ[m] :=
				Sum[Binomial[n,k] (m!/(m-k)!) ((-1)^k)(Commutator[P,Q]^k) ** P^(m-k)**Q^(n-k),{k,1,n}];
	  	 
			Commutator[P,x_] /; checkOperator[Q,x] := Commutator[P,Q]**
				(D[x,Q]/.Times->NonCommutativeMultiply) // MapAll[QAProductExpand,#]&;
		  
			Commutator[x_,P] /; checkOperator[Q,x] := -Commutator[P,Q]**
				(D[x,Q]/.Times->NonCommutativeMultiply) // MapAll[QAProductExpand,#]&;
		   
			Commutator[Q,y_] /; checkOperator[P,y] := -Commutator[P,Q]**
				(D[y,P]/.Times->NonCommutativeMultiply) // MapAll[QAProductExpand,#]&;

			Commutator[y_,Q] /; checkOperator[P,y] := Commutator[P,Q]**
				(D[y,P]/.Times->NonCommutativeMultiply) // MapAll[QAProductExpand,#]&; ]

QASet[{ SetDelayed[( Commutator[P_,Commutator[P_,Q_]] | Commutator[P_,Commutator[Q_,P_]] | 
					 Commutator[Commutator[P_,Q_],P_] | Commutator[Commutator[Q_,P_],P_]  ), 0 ],
		SetDelayed[( Commutator[Q_,Commutator[P_,Q_]] | Commutator[Q_,Commutator[Q_,P_]] |
					 Commutator[Commutator[P_,Q_],Q_] | Commutator[Commutator[Q_,P_],Q_]  ), 0 ] },
		OptionsPattern[] ] :=
	
	Module[ {k},

			Commutator[P,Commutator[P,Q]] := 0;
			Commutator[P,Commutator[Q,P]] := 0;
			Commutator[Commutator[P,Q],P] := 0;
			Commutator[Commutator[Q,P],P] := 0;
			Commutator[Q,Commutator[P,Q]] := 0;
			Commutator[Q,Commutator[Q,P]] := 0;
			Commutator[Commutator[P,Q],Q] := 0;
			Commutator[Commutator[Q,P],Q] := 0;
			
			Commutator[P^n_,Q^m_]  /; IntegerQ[n] && IntegerQ[m] := 
				Sum[Binomial[n,k] (m!/(m-k)!) (Commutator[P,Q]^k) ** Q^(m-k)**P^(n-k),{k,1,n}];

			Commutator[Q^n_,P^m_] /; IntegerQ[n] && IntegerQ[m] :=
				Sum[Binomial[n,k] (m!/(m-k)!) ((-1)^k)(Commutator[P,Q]^k) ** P^(m-k)**Q^(n-k),{k,1,n}];
	  	   
			Commutator[P,x_] /; checkOperator[Q,x] := Commutator[P,Q]**
				(D[x,Q]/.Times->NonCommutativeMultiply) // MapAll[QAProductExpand,#]&;

			Commutator[x_,P] /; checkOperator[Q,x] := -Commutator[P,Q]**
				(D[x,Q]/.Times->NonCommutativeMultiply) // MapAll[QAProductExpand,#]&;

			Commutator[Q,y_] /; checkOperator[P,y] := -Commutator[P,Q]**
				(D[y,P]/.Times->NonCommutativeMultiply) // MapAll[QAProductExpand,#]&;

			Commutator[y_,Q] /; checkOperator[P,y] := Commutator[P,Q]**
				(D[y,P]/.Times->NonCommutativeMultiply) // MapAll[QAProductExpand,#]&; ]

QAClear[{Unset[(Commutator[P_,Commutator[P_,Q_]] | Commutator[P_,Commutator[Q_,P_]] |
				Commutator[Commutator[P_,Q_],P_] | Commutator[Commutator[Q_,P_],P_] )],
		 Unset[(Commutator[Q_,Commutator[P_,Q_]] | Commutator[Q_,Commutator[Q_,P_]] | 
				Commutator[Commutator[P_,Q_],Q_] | Commutator[Commutator[Q_,P_],Q_] )]}] :=

	Module[ {},	
		   Commutator[P,Commutator[P,Q]] =.;
		   Commutator[P,Commutator[Q,P]] =.;
		   Commutator[Commutator[P,Q],P] =.;
		   Commutator[Commutator[Q,P],P] =.;
		   Commutator[Q,Commutator[P,Q]] =.;
		   Commutator[Q,Commutator[Q,P]] =.;
		   Commutator[Commutator[P,Q],Q] =.;
		   Commutator[Commutator[Q,P],Q] =.;
		   Commutator[P^n_,Q^m_] /; IntegerQ[n] && IntegerQ[m] =.;
		   Commutator[Q^n_,P^m_] /; IntegerQ[n] && IntegerQ[m] =.;
	  	 Commutator[P,x_] /; checkOperator[Q,x] =.;
		   Commutator[x_,P] /; checkOperator[Q,x] =.;
		   Commutator[Q,y_] /; checkOperator[P,y] =.;
		   Commutator[y_,Q] /; checkOperator[P,y] =.; 
	]


QASet[ lexpr_ = rexpr_, OptionsPattern[]] := Module[{temp}, 
    lexpr = temp = rexpr;
	Evaluate[ HoldPattern[lexpr] /. {Ket->Bra, Bra->Ket} ] = Hermitian[temp];  ] /; KetQ[lexpr] || BraQ[lexpr]

QASet[ lexpr_ := rexpr_, OptionsPattern[]] := Module[{},
	If[Head[rexpr]===Condition,
	Evaluate[ HoldPattern[lexpr] /. {Ket->Bra, Bra->Ket} ] := Hermitian[ rexpr[[1]] ] /; rexpr[[2]],
	Evaluate[ HoldPattern[lexpr] /. {Ket->Bra, Bra->Ket} ] := Hermitian[rexpr]];	
	lexpr := rexpr; ] /; KetQ[lexpr] || BraQ[lexpr]

QAClear[ lexpr_ =. ] := Module[{},
	lexpr =. ; 
	ReleaseHold[ Hold[lexpr =. ] /. {Ket->Bra, Bra->Ket}]] /; KetQ[lexpr] || BraQ[lexpr]


(* ::Subsection:: *)
(*Exponential Operators*)


Power[E,a_. Operator[1]] := Exp[a] /; ConstantQ[a] 

 (* Power/: Power[E, P_] ** Power[E,d_.+ Q_] := Power[E,d] Power[E,P] ** 
					Power[E,Q] /; OperatorQ[P] && OperatorQ[Q] &&
							ConstantQ[c]&& ConstantQ[d] *)

 Power/: Power[E, P_] ** Power[E, Q_] := Operator[1] /; OperatorQ[P] && OperatorQ[Q] && P == -Q
 
 Power/: Power[E,c_. P_] ** Power[E,d_. P_] := Exp[(c+d) P] /; OperatorQ[P] 

 HoldPattern[Exp[a_. Operator[0]]] := Operator[1]

 (* Power/: Power[E,P_] ** Power[E,c_] := Exp[c] ** Exp[P] /; OperatorQ[P] && ConstantQ[c] *)

 SetAttributes[BCHRelation,HoldAllComplete]
 
 BCHRelation[expr_] := 
	Module[ {sol}, 
		ClearAttributes[Plus, Orderless] ;
		sol = expr /.  Power[E, c_. (P_ + Q_)]  :> 
		Module[{ com1 = c^2 Commutator[P,Q], com2 = c^2 Commutator[Q,P]},        
		   If[Commutator[c P, com1] === Commutator[c Q, com1] === 0 ||
		      Commutator[c P, com2] === Commutator[c Q, com2] === 0 ||
		      Commutator[c P, com2] === Commutator[c Q, com1] === 0 ||
		      Commutator[c P, com1] === Commutator[c Q, com2] === 0 ||              
              Commutator[c P, com1] === Commutator[com1, c Q] === 0 ||
		      Commutator[c P, com2] === Commutator[com2, c Q] === 0 ||
              Commutator[c P, com1] === Commutator[com2, c Q] === 0 ||
		      Commutator[c P, com2] === Commutator[com1, c Q] === 0 ||
              Commutator[com1, c P] === Commutator[c Q, com1] === 0 ||
		      Commutator[com1, c P] === Commutator[com1, c Q] === 0 ||
              Commutator[com1, c P] === Commutator[c Q, com2] === 0 ||
		      Commutator[com2, c P] === Commutator[com1, c Q] === 0 ||	     
		      Commutator[com2, c P] === Commutator[c Q, com2] === 0 ||
		      Commutator[com2, c P] === Commutator[com2, c Q] === 0 ||
              Commutator[com2, c P] === Commutator[c Q, com1] === 0 ||
		      Commutator[com1, c P] === Commutator[com2, c Q] === 0 ,
			  QAProductExpand[ Exp[c P] ** Exp[c Q] ** Exp[-com1/2] ], 
			  Exp[c (P + Q)]]]  /; OperatorQ[P] && OperatorQ[Q] && ConstantQ[c] ;
        SetAttributes[Plus,Orderless];
		sol
	]


(* ::Subsection:: *)
(*Series Expansion*)


QASeries[expr_, {x_,a_, n_}] := Module[{k,t,f},
		t = (((expr /. x->a) /. Operator[1]->1) Operator[1] +
		       Sum[(( (D[expr,{x,k}] /. x->a) /. 
			    Operator[1]->1) (x-a)^k) / f[k], {k,1,n}]) /. Operator[0] -> 0;
		t = t /. f[v_] :> HoldForm[v!]]

QASeries[c_. Exp[P_] ** R_ ** Exp[Q_], n_] :=
		Module[{f, k, s},
		    s = c Sum[Nest[Commutator[P,#]&, R, k]/f[k], {k, 0, n-1}];
		    s = s /. f[v_] :> HoldForm[v!]] /; 
			OperatorQ[P] && OperatorQ[R] && P === -Q && Commutator[P,R] =!= 0




(* ::Subsection:: *)
(*Manipulate Quantum Ojects*)


QACollect[expr_,{x__}]  := Collect[expr,{x}] /. Operator[1]->1  

QAOrdering[ expr_, list_List ] := expr //. CommutationOrderRules[list]

QAPowerExpand[expr_] := expr //. rule1

QAPowerContract[expr_] := expr //. rule2

QACommutatorExpand[expr_?OperatorQ] := expr //. rule3 
 
QACommutatorExpand[expr_] := expr

QAProductExpand[expr_?(!ConstantQ[#]&)] := (expr //. rule5) //. rule4
 
QAProductExpand[expr_?ConstantQ] := expr //. NonCommutativeMultiply -> Times
 
QAProductExpand[expr_] := expr

QAExpandAll[expr_?OperatorQ] := expr // QAPowerExpand // QAProductExpand //
					QACommutatorExpand //Expand

QAExpandAll[expr_] := expr

QACommutatorContract[expr_, t : {{_, _} ...}] := 
	Fold[ReplaceRepeated, expr,                        			
		Map[First[#] ** Last[#] :>
                           Commutator[First[#], Last[#]] + Last[#] ** First[#] &, t]]

QACommutatorContract[expr_,{P_,Q_}] := expr //. P**Q :> Commutator[P,Q] + Q**P

QACommutatorContract[expr_] := expr //. P_**Q_ - Q_**P_ :> Commutator[P,Q]


(* ::Subsection:: *)
(*Auxiliar Definitions*)


CommutationOrderRules[l_List] := 
	Module[{x = {First[l]}, y = Rest[l], f, sol = {}},
	   f = Rule[#2**#1, #1**#2 + Commutator[#2, #1]] &;
	   While[Length[y] >= 1, 
		   AppendTo[sol, Outer[f, x, y]];
	    	   x = {First[y]};
		   y = Rest[y]];
		   Flatten[sol]]

expand[ NonCommutativeMultiply[P__]] :=
	Module[{temp1, temp2, l1,  temp3 = {P}},
		temp1 = Select[temp3,!ConstantQ[#]&];	
		temp2 = Select[temp3, ConstantQ[#]&];
		l1 = Length[temp1];
	        Which[l1 > 1, Times[ Apply[NonCommutativeMultiply,temp1], Apply[Times,temp2]],
	      	      l1===1, First[temp1] Apply[Times,temp2],
	      	      l1===0, Apply[Times,temp2]]]

rule1 = {(P_Plus)^n_?IntegerQ /; (OperatorQ[P]||KetQ[P]||BraQ[P]) && n > 0 :> Plus @@ Flatten[
                      Outer[NonCommutativeMultiply,Sequence @@ Table[List @@ P,{n}]]],
          Power[P_, n_?IntegerQ] /; (OperatorQ[P]||KetQ[P]||BraQ[P]) && n > 0 :>
		     Apply[NonCommutativeMultiply, Table[P,{n}]],
          Power[P_, n_?IntegerQ] /; (OperatorQ[P]||KetQ[P]||BraQ[P]) && n < -1 :>
		     Apply[NonCommutativeMultiply, Table[P^(-1),{-n}]]}

rule2 = P_^n_. ** P_^m_. :> P^(n+m) /; OperatorQ[P]

rule3 = HoldPattern[Commutator[P_, Q_]] :> P ** Q - Q ** P;

rule4 = NonCommutativeMultiply[P__] :> expand[NonCommutativeMultiply[P]];

rule5 = HoldPattern[Bra[x__] ** Ket[y__]] :> Braket[{x}, {y}];

checkOperator[Operator[P__],expr_] :=
    Complement[Cases[expr, (Operator[__]|Operator[__][__]), {0, Infinity}], 
	      {Operator[0],Operator[1]} ] ===
	      Union[Cases[expr, Operator[P], {0, Infinity}]] && expr =!= Operator[P]

checkOperator[Operator[P__][t__],expr_] :=
    Complement[Cases[expr, (Operator[__]|Operator[__][__]), {0, Infinity}],
	      {Operator[0],Operator[1]} ] ===
            Union[Cases[expr, Operator[P][t], {0, Infinity}]] && expr =!= Operator[P][t]


(* ::Subsection:: *)
(*Define New Notations*)


(* ::Subsubsection:: *)
(*Begin*)


If[!$Notebooks, 
	Print["Not a notebook front end. Notations are not loaded."];
	Goto[endNotations]; ];

If[ !ValueQ[AutoLoadNotationPalette], AutoLoadNotationPalette = True]
	
If[	AutoLoadNotationPalette === True,
	NotebookOpen[ToFileName[{DirectoryName[$InputFileName], "FrontEnd", "Palettes"},"QAPalette.nb"]];
]

If[ !ValueQ[AutoLoadNavigatorPalette], AutoLoadNavigatorPalette = True]
	
If[	AutoLoadNavigatorPalette === True,
	NotebookOpen[ToFileName[{DirectoryName[$InputFileName], "FrontEnd", "Palettes"},"QANavigator.nb"]];
]

nb1 = EvaluationNotebook[];
style = BaseStyle;


(* ::Subsubsection:: *)
(*Box Manipulations*)


inputToBoxes[x__, form_] := 
  RowBox[Drop[Flatten[Map[{",", MakeBoxes[#, form]} &, {x}]], 1]];

rowBoxToArgBoxes[x_] := Block[{Message},
   SetAttributes[Message, HoldFirst];
   If[Head[x] === RowBox,
    If[ToExpression[x] === $Failed, Sequence @@ x[[1]], x], x]];


(* ::Subsubsection:: *)
(*Modify Parenthesize and FullForm*)


SubscriptBox[RowBox[{x_, "[", RowBox[y_], "]"}], "-"] := 
  RowBox[{x, " ", SubscriptBox[RowBox[{"[", RowBox[y], "]"}], "-"]}];

RowBox[{"(", TagBox[arg_, Ket, opts___], ")"}] := TagBox[arg, Ket, opts];

RowBox[{"(", TagBox[arg_, Bra, opts___], ")"}] := TagBox[arg, Bra, opts];

RowBox[{"(", TagBox[arg_, Braket, opts___], ")"}] := TagBox[ arg, Braket, opts];

FullForm/: TagBox[StyleBox[expr_, Verbatim[ShowSpecialCharacters -> False], 
	rest__], FullForm] := TagBox[StyleBox[expr, rest], FullForm]; 


(* ::Subsubsection:: *)
(*Define Operator Notation*)


(*Power/:*) MakeBoxes[Power[Operator[x_], n_?Negative], form_] :=
  SuperscriptBox[OverscriptBox[MakeBoxes[x, form], "^"], n];

(*Power/:*) MakeBoxes[Power[Operator[x_, {y__}, {}],n_?Negative], form_] :=
  SubsuperscriptBox[OverscriptBox[MakeBoxes[x, form], "^"], inputToBoxes[y, form], n];

(*Power/:*) MakeBoxes[Power[Operator[x_,{},{z__}], n_], form_] :=
  SuperscriptBox[ RowBox[{"(", TagBox[SuperscriptBox[OverscriptBox[inputToBoxes[x, form],
	 "^"], inputToBoxes[z,form]],superscript] ,")"}], inputToBoxes[n,form]];

(*Power/:*)MakeBoxes[Power[Operator[x_, {y__}, {z__}], n_], form_] :=
  SuperscriptBox[ RowBox[{"(",TagBox[SubsuperscriptBox[OverscriptBox[MakeBoxes[x, form], "^"], 
inputToBoxes[y, form], inputToBoxes[z, form]],superscript],")"}], inputToBoxes[n,form]];


Operator /: MakeBoxes[Operator[x_], form_] :=
  OverscriptBox[MakeBoxes[x, form], "^"];

MakeExpression[
   OverscriptBox[x_, "^"], form_] :=
  MakeExpression[RowBox[{"Operator", "[", x, "]"}], form];

Operator /: MakeBoxes[Operator[x_, {}, {y__}], form_] :=
  TagBox[SuperscriptBox[OverscriptBox[MakeBoxes[x, form], "^"], 
   inputToBoxes[y, form]], superscript];

MakeExpression[
   TagBox[SuperscriptBox[OverscriptBox[x_, "^"], y_], superscript], form_] :=
  MakeExpression[
   RowBox[{"Operator", "[", 
     RowBox[{x, ",", RowBox[{"{", "}"}], ",", RowBox[{"{", y, "}"}]}],
      "]"}], form];

Operator /: MakeBoxes[Operator[x_, {y__}, {}], form_] :=
  SubscriptBox[OverscriptBox[MakeBoxes[x, form], "^"], 
   inputToBoxes[y, form]];

MakeExpression[
   SubscriptBox[OverscriptBox[x_, "^"], y_], form_] :=
  MakeExpression[
   RowBox[{"Operator", "[", 
     RowBox[{x, ",", RowBox[{"{", y, "}"}], ",", RowBox[{"{", "}"}]}],
      "]"}], form];

Operator /: MakeBoxes[Operator[x_, {y__}, {z__}], form_] :=
  TagBox[SubsuperscriptBox[OverscriptBox[MakeBoxes[x, form], "^"], 
   inputToBoxes[y, form], inputToBoxes[z, form]],superscript];

MakeExpression[
   TagBox[SubsuperscriptBox[OverscriptBox[x_, "^"], y_, z_],superscript], form_] :=
  MakeExpression[
   RowBox[{"Operator", "[", 
     RowBox[{x, ",", RowBox[{"{", y, "}"}], ",", 
       RowBox[{"{", z, "}"}]}], "]"}], form];

Operator /: MakeBoxes[Operator[x_, "H"], form_] :=
  SuperscriptBox[OverscriptBox[MakeBoxes[x, form], "^"], "\[Dagger]"];
 
MakeExpression[
   SuperscriptBox[OverscriptBox[x_, "^"], "\[Dagger]"], form_] :=
  MakeExpression[
   RowBox[{"Operator", "[", RowBox[{x, ",", "\"H\""}], "]"}], form];

Operator /: MakeBoxes[Operator[x_, "H", {}, {y__}], form_] :=
  SuperscriptBox[OverscriptBox[MakeBoxes[x, form], "^"], 
   RowBox[{inputToBoxes[y, form], ",", "\[Dagger]"}]];

MakeExpression[
   SuperscriptBox[OverscriptBox[x_, "^"], 
    RowBox[{y__, ",", "\[Dagger]"}]], form_] :=
  MakeExpression[
   RowBox[{"Operator", "[", 
     RowBox[{x, ",", "\"H\"", ",", RowBox[{"{", "}"}], ",", 
       RowBox[{"{", y, "}"}]}], "]"}], form];

Operator /: MakeBoxes[Operator[x_, "H", {y__}, {}], form_] :=
  SubsuperscriptBox[OverscriptBox[MakeBoxes[x, form], "^"], 
   inputToBoxes[y, form], "\[Dagger]"];

MakeExpression[
   SubsuperscriptBox[OverscriptBox[x_, "^"], y_, "\[Dagger]"], 
   form_] :=
  MakeExpression[
   RowBox[{"Operator", "[", 
     RowBox[{x, ",", "\"H\"", ",", RowBox[{"{", y, "}"}], ",", 
       RowBox[{"{", "}"}]}], "]"}], form];

Operator /: MakeBoxes[Operator[x_, "H", {y__}, {z__}], form_] :=
  SubsuperscriptBox[OverscriptBox[MakeBoxes[x, form], "^"], 
   inputToBoxes[y, form], 
   RowBox[{inputToBoxes[z, form], ",", "\[Dagger]"}]];

MakeExpression[
   SubsuperscriptBox[OverscriptBox[x_, "^"], y_, 
    RowBox[{z__, ",", "\[Dagger]"}]], form_] :=
  MakeExpression[
   RowBox[{"Operator", "[", 
     RowBox[{x, ",", "\"H\"", ",", RowBox[{"{", y, "}"}], ",", 
       RowBox[{"{", z, "}"}]}], "]"}], form];
       
(* ::Subsubsection:: *)
(*Define Tensor Product Notation*)     

MakeExpression[RowBox[{
      UnderoverscriptBox["\[CircleTimes]", RowBox[{a_, "=", b_}], c_], expr_}], form_] :=
  MakeExpression[ RowBox[{"DoTensorProduct", "[", 
    RowBox[{expr, ",", 
       RowBox[{"{", 
          RowBox[{a, ",", b, ",", c}], "}"}]}], "]"}], form]

(* ::Subsubsection:: *)
(*Define Bra and Ket Notation*)


Ket /: MakeBoxes[Ket[x__], form_] := 
  TagBox[RowBox[{AdjustmentBox["\[VerticalSeparator]", 
      BoxMargins -> {{-0.2, 0}, {0, 0}},
      BoxBaselineShift -> -0.1], 
     AdjustmentBox[
      TagBox[inputToBoxes[x, form], KetArgs, style -> "BraKetArg"], 
      BoxBaselineShift -> 0], "\[RightAngleBracket]"}], Ket, 
   style -> "KetWrapper"];

MakeExpression[
   TagBox[
    RowBox[{AdjustmentBox["\[VerticalSeparator]", opts1___], 
      AdjustmentBox[TagBox[x_, KetArgs, opts2___], opts3___], 
      "\[RightAngleBracket]"}], Ket, opts4___], form_] := 
  MakeExpression[RowBox[{"Ket", "[", x, "]"}], form];

Bra /: MakeBoxes[Bra[x__], form_] := 
  TagBox[RowBox[{"\[LeftAngleBracket]", 
     AdjustmentBox[
      TagBox[inputToBoxes[x, form], BraArgs, style -> "BraKetArg"], 
      BoxBaselineShift -> 0], 
     AdjustmentBox["\[VerticalSeparator]", 
      BoxMargins -> {{0, -0.2}, {0, 0}},
      BoxBaselineShift -> -0.1]}], Bra, 
   style -> "BraWrapper"];

MakeExpression[
   TagBox[
    RowBox[{"\[LeftAngleBracket]", 
      AdjustmentBox[TagBox[x_, BraArgs, opts1___], opts2___], 
      AdjustmentBox["\[VerticalSeparator]", opts3___]}], Bra, 
    opts4___], form_] :=
  MakeExpression[RowBox[{"Bra", "[", x, "]"}], form];

Braket /: MakeBoxes[Braket[{x__}, {y__}], form_] := 
  TagBox[RowBox[{"\[LeftAngleBracket]", 
     AdjustmentBox[
      RowBox[{TagBox[inputToBoxes[x, form], BraKetArgs, 
         style -> "BraKetArg"], "\[VerticalSeparator]", 
        TagBox[inputToBoxes[y, form], BraKetArgs, 
         style -> "BraKetArg"]}], BoxBaselineShift -> 0], 
     "\[RightAngleBracket]"}], Braket, style -> "BraKetWrapper"];

MakeExpression[
   TagBox[
    RowBox[{"\[LeftAngleBracket]", 
      AdjustmentBox[
       RowBox[{TagBox[x_, BraKetArgs, opts1___], 
         "\[VerticalSeparator]", TagBox[y_, BraKetArgs, opts2___]}], 
       opts3___], "\[RightAngleBracket]"}], Braket, opts4___], 
   form_] :=
  MakeExpression[
   RowBox[{"Braket", "[", 
     RowBox[{RowBox[{"{", x, "}"}], ",", RowBox[{"{", y, "}"}]}], 
     "]"}], form];

Ket /: MakeBoxes[Ket[x__, {y__}, {}], form_] :=
  TagBox[
   SubscriptBox[
    RowBox[
     {AdjustmentBox["\[VerticalSeparator]", 
       BoxMargins -> {{-0.2, 0}, {0, 0}},
       BoxBaselineShift -> -0.1], 
      AdjustmentBox[
       TagBox[inputToBoxes[x, form], KetArgs, style -> "BraKetArg"], 
       BoxBaselineShift -> 0], "\[RightAngleBracket]"}], 
    TagBox[
     AdjustmentBox[
      TagBox[inputToBoxes[y, form], KetArgs, style -> "BraKetArg"], 
      BoxBaselineShift -> 0], UDScript, style -> "RIndex"]], Ket, 
   BaseStyle -> "KetWrapper"];

MakeExpression[
   TagBox[
    SubscriptBox[
     RowBox[
      {AdjustmentBox["\[VerticalSeparator]", opts1___], 
       AdjustmentBox[TagBox[x_, KetArgs, opts2___], opts3___], 
       "\[RightAngleBracket]"}], 
     TagBox[AdjustmentBox[TagBox[y_, KetArgs, opts4___], opts5___], 
      UDScript, opts6___]], Ket, opts7___], form_] :=
  MakeExpression[
   RowBox[{"Ket", "[", 
     RowBox[{rowBoxToArgBoxes[x], ",", 
       RowBox[{"{", y, "}"}], ",", RowBox[{"{", "}"}]}], "]"}], form];
 

Ket /: MakeBoxes[Ket[x__, {}, {y__}], form_] :=
  TagBox[SuperscriptBox[
    RowBox[{AdjustmentBox["\[VerticalSeparator]", 
       BoxMargins -> {{-0.2, 0}, {0, 0}},
       BoxBaselineShift -> -0.1], 
      AdjustmentBox[
       TagBox[inputToBoxes[x, form], KetArgs, style -> "BraKetArg"], 
       BoxBaselineShift -> 0], "\[RightAngleBracket]"}], 
    TagBox[
     AdjustmentBox[
      TagBox[inputToBoxes[y, form], KetArgs, style -> "BraKetArg"], 
      BoxBaselineShift -> 0], UDScript, style -> "RIndex"]], Ket, 
   style -> "KetWrapper"];

MakeExpression[
   TagBox[
    SuperscriptBox[
     RowBox[{AdjustmentBox["\[VerticalSeparator]", opts1___], 
       AdjustmentBox[TagBox[x_, KetArgs, opts2___], opts3___], 
       "\[RightAngleBracket]"}], 
     TagBox[AdjustmentBox[TagBox[y_, KetArgs, opts4___], opts5___], 
      UDScript, opts6___]], Ket, opts7___], form_] :=
  MakeExpression[
   RowBox[{"Ket", "[", 
     RowBox[{rowBoxToArgBoxes[x], ",", 
       RowBox[{"{", "}"}], ",", RowBox[{"{", y, "}"}]}], "]"}], form];
 

Ket /: MakeBoxes[Ket[x__, {y__}, {z__}], form_] :=
  TagBox[SubsuperscriptBox[
    RowBox[{AdjustmentBox["\[VerticalSeparator]", 
       BoxMargins -> {{-0.2, 0}, {0, 0}},
       BoxBaselineShift -> -0.1], 
      AdjustmentBox[
       TagBox[inputToBoxes[x, form], KetArgs, style -> "BraKetArg"], 
       BoxBaselineShift -> 0], "\[RightAngleBracket]"}], 
    TagBox[
     AdjustmentBox[
      TagBox[inputToBoxes[y, form], KetArgs, style -> "BraKetArg"], 
      BoxBaselineShift -> 0], UDScript, style -> "RIndex"], 
    TagBox[AdjustmentBox[
      TagBox[inputToBoxes[z, form], KetArgs, style -> "BraKetArg"], 
      BoxBaselineShift -> 0], UDScript, style -> "RIndex"]], Ket, 
   style -> "KetWrapper"];

MakeExpression[
   TagBox[
    SubsuperscriptBox[
     RowBox[{AdjustmentBox["\[VerticalSeparator]", opts1___], 
       AdjustmentBox[TagBox[x_, KetArgs, opts2___], opts3___], 
       "\[RightAngleBracket]"}], 
     TagBox[AdjustmentBox[TagBox[y_, KetArgs, opts4___], opts5___], 
      UDScript, opts6___], 
     TagBox[AdjustmentBox[TagBox[z_, KetArgs, opts7___], opts8___], 
      UDScript, opts9___]], Ket, opts10___], form_] :=
  MakeExpression[
   RowBox[{"Ket", "[", 
     RowBox[{rowBoxToArgBoxes[x], ",", 
       RowBox[{"{", y, "}"}], ",", RowBox[{"{", z, "}"}]}], "]"}], 
   form];

Bra /: MakeBoxes[Bra[x__, {y__}, {}], form_] :=
  TagBox[RowBox[{SubscriptBox["\[InvisiblePrefixScriptBase]", 
      TagBox[AdjustmentBox[
        TagBox[inputToBoxes[y, form], BraArgs, style -> "BraKetArg"], 
        BoxBaselineShift -> 0], UDScript, style -> "LIndex"]], 
     RowBox[{"\[LeftAngleBracket]", 
       AdjustmentBox[
        TagBox[inputToBoxes[x, form], BraArgs, style -> "BraKetArg"], 
        BoxBaselineShift -> 0], 
       AdjustmentBox["\[VerticalSeparator]", 
        BoxMargins -> {{0, -0.2}, {0, 0}},
        BoxBaselineShift -> -0.1]}]}], Bra, 
   style -> "BraWrapper"];

MakeExpression[
   TagBox[
    RowBox[{SubscriptBox["\[InvisiblePrefixScriptBase]", 
       TagBox[AdjustmentBox[TagBox[y_, BraArgs, opts1___], opts2___], 
        UDScript, opts3___]], 
      RowBox[{"\[LeftAngleBracket]", 
        AdjustmentBox[TagBox[x_, BraArgs, opts4___], opts5___], 
        AdjustmentBox["\[VerticalSeparator]", opts6___]}]}], Bra, 
    opts7___], form_] :=
  MakeExpression[
   RowBox[{"Bra", "[", 
     RowBox[{rowBoxToArgBoxes[x], ",", 
       RowBox[{"{", y, "}"}], ",", RowBox[{"{", "}"}]}], "]"}], form];
 

Bra /: MakeBoxes[Bra[x__, {}, {y__}], form_] :=
  TagBox[RowBox[{SuperscriptBox["\[InvisiblePrefixScriptBase]", 
      TagBox[AdjustmentBox[
        TagBox[inputToBoxes[y, form], BraArgs, style -> "BraKetArg"], 
        BoxBaselineShift -> 0], UDScript, style -> "LIndex"]], 
     RowBox[{"\[LeftAngleBracket]", 
       AdjustmentBox[
        TagBox[inputToBoxes[x, form], BraArgs, style -> "BraKetArg"], 
        BoxBaselineShift -> 0], 
       AdjustmentBox["\[VerticalSeparator]", 
        BoxMargins -> {{0, -0.2}, {0, 0}},
        BoxBaselineShift -> -0.1]}]}], Bra, 
   style -> "BraWrapper"];

MakeExpression[
   TagBox[
    RowBox[{SuperscriptBox["\[InvisiblePrefixScriptBase]", 
       TagBox[AdjustmentBox[TagBox[y_, BraArgs, opts1___], opts2___], 
        UDScript, opts3___]], 
      RowBox[{"\[LeftAngleBracket]", 
        AdjustmentBox[TagBox[x_, BraArgs, opts4___], opts5___], 
        AdjustmentBox["\[VerticalSeparator]", opts6___]}]}], Bra, 
    opts7___], form_] :=
  MakeExpression[
   RowBox[{"Bra", "[", 
     RowBox[{rowBoxToArgBoxes[x], ",", 
       RowBox[{"{", "}"}], ",", RowBox[{"{", y, "}"}]}], "]"}], form];
 

Bra /: MakeBoxes[Bra[x__, {y__}, {z__}], form_] :=
  TagBox[RowBox[{SubsuperscriptBox["\[InvisiblePrefixScriptBase]", 
      TagBox[AdjustmentBox[
        TagBox[inputToBoxes[y, form], BraArgs, style -> "BraKetArg"], 
        BoxBaselineShift -> 0], UDScript, style -> "LIndex"], 
      TagBox[AdjustmentBox[
        TagBox[inputToBoxes[z, form], BraArgs, style -> "BraKetArg"], 
        BoxBaselineShift -> 0], UDScript, style -> "LIndex"]], 
     RowBox[{"\[LeftAngleBracket]", 
       AdjustmentBox[
        TagBox[inputToBoxes[x, form], BraArgs, style -> "BraKetArg"], 
        BoxBaselineShift -> 0], 
       AdjustmentBox["\[VerticalSeparator]", 
        BoxMargins -> {{0, -0.2}, {0, 0}},
        BoxBaselineShift -> -0.1]}]}], Bra, 
   style -> "BraWrapper"];

MakeExpression[
   TagBox[
    RowBox[{SubsuperscriptBox["\[InvisiblePrefixScriptBase]", 
       TagBox[AdjustmentBox[TagBox[y_, BraArgs, opts1___], opts2___], 
        UDScript, opts3___], 
       TagBox[AdjustmentBox[TagBox[z_, BraArgs, opts4___], opts5___], 
        UDScript, opts6___]], 
      RowBox[{"\[LeftAngleBracket]", 
        AdjustmentBox[TagBox[x_, BraArgs, opts7___], opts8___], 
        AdjustmentBox["\[VerticalSeparator]", opts9___]}]}], Bra, 
    opts10___], form_] :=
  MakeExpression[
   RowBox[{"Bra", "[", 
     RowBox[{rowBoxToArgBoxes[x], ",", 
       RowBox[{"{", y, "}"}], ",", RowBox[{"{", z, "}"}]}], "]"}], 
   form];

Braket /: MakeBoxes[Braket[{x__}, {y__, {z__}, {}}], form_] :=
  TagBox[SubscriptBox[
    RowBox[{"\[LeftAngleBracket]", 
      AdjustmentBox[
       RowBox[{TagBox[inputToBoxes[x, form], BraKetArgs, 
          style -> "BraKetArg"], "\[VerticalSeparator]", 
         TagBox[inputToBoxes[y, form], BraKetArgs, 
          style -> "BraKetArg"]}], BoxBaselineShift -> 0], 
      "\[RightAngleBracket]"}], 
    TagBox[AdjustmentBox[
      TagBox[inputToBoxes[z, form], KetArgs, style -> "BraKetArg"], 
      BoxBaselineShift -> 0], UDScript, style -> "RIndex"]], Braket, 
   style -> "BraKetWrapper"];

MakeExpression[
   TagBox[
    SubscriptBox[
     RowBox[{"\[LeftAngleBracket]", 
       AdjustmentBox[
        RowBox[{TagBox[x_, BraKetArgs, opts1___], 
          "\[VerticalSeparator]", TagBox[y_, BraKetArgs, opts2___]}], 
        opts3___], "\[RightAngleBracket]"}], 
     TagBox[AdjustmentBox[TagBox[z_, KetArgs, opts4___], opts5___], 
      UDScript, opts6___]], Braket, opts7___], form_] :=
  MakeExpression[
   RowBox[{"Braket", "[", 
     RowBox[{RowBox[{"{", x, "}"}], ",", 
       RowBox[{"{", 
         RowBox[{rowBoxToArgBoxes[y], ",", 
           RowBox[{"{", z, "}"}], ",", RowBox[{"{", "}"}]}], "}"}]}], 
     "]"}], form];

Braket /: MakeBoxes[Braket[{x__}, {y__, {}, {z__}}], form_] :=
  TagBox[SuperscriptBox[
    RowBox[{"\[LeftAngleBracket]", 
      AdjustmentBox[
       RowBox[{TagBox[inputToBoxes[x, form], BraKetArgs, 
          style -> "BraKetArg"], "\[VerticalSeparator]", 
         TagBox[inputToBoxes[y, form], BraKetArgs, 
          style -> "BraKetArg"]}], BoxBaselineShift -> 0], 
      "\[RightAngleBracket]"}], 
    TagBox[AdjustmentBox[
      TagBox[inputToBoxes[z, form], KetArgs, style -> "BraKetArg"], 
      BoxBaselineShift -> 0], UDScript, style -> "RIndex"]], Braket, 
   style -> "BraKetWrapper"];

MakeExpression[
   TagBox[
    SuperscriptBox[
     RowBox[{"\[LeftAngleBracket]", 
       AdjustmentBox[
        RowBox[{TagBox[x_, BraKetArgs, opts1___], 
          "\[VerticalSeparator]", TagBox[y_, BraKetArgs, opts2___]}], 
        opts3___], "\[RightAngleBracket]"}], 
     TagBox[AdjustmentBox[TagBox[z_, KetArgs, opts4___], opts5___], 
      UDScript, opts6___]], Braket, opts7___], form_] :=
  MakeExpression[
   RowBox[{"Braket", "[", 
     RowBox[{RowBox[{"{", x, "}"}], ",", 
       RowBox[{"{", 
         RowBox[{rowBoxToArgBoxes[y], ",", 
           RowBox[{"{", "}"}], ",", RowBox[{"{", z, "}"}]}], "}"}]}], 
     "]"}], form];

Braket /: MakeBoxes[Braket[{w__}, {x__, {y__}, {z__}}], form_] :=
  TagBox[SubsuperscriptBox[
    RowBox[{"\[LeftAngleBracket]", 
      AdjustmentBox[
       RowBox[{TagBox[inputToBoxes[w, form], BraKetArgs, 
          style -> "BraKetArg"], "\[VerticalSeparator]", 
         TagBox[inputToBoxes[x, form], BraKetArgs, 
          style -> "BraKetArg"]}], BoxBaselineShift -> 0], 
      "\[RightAngleBracket]"}], 
    TagBox[AdjustmentBox[
      TagBox[inputToBoxes[y, form], KetArgs, style -> "BraKetArg"], 
      BoxBaselineShift -> 0], UDScript, style -> "RIndex"], 
    TagBox[AdjustmentBox[
      TagBox[inputToBoxes[z, form], KetArgs, style -> "BraKetArg"], 
      BoxBaselineShift -> 0], UDScript, style -> "RIndex"]], Braket, 
   style -> "BraKetWrapper"];

MakeExpression[
   TagBox[
    SubsuperscriptBox[
     RowBox[{"\[LeftAngleBracket]", 
       AdjustmentBox[
        RowBox[{TagBox[w_, BraKetArgs, opts1___], 
          "\[VerticalSeparator]", TagBox[x_, BraKetArgs, opts2___]}], 
        opts3___], "\[RightAngleBracket]"}], 
     TagBox[AdjustmentBox[TagBox[y_, KetArgs, opts4___], opts5___], 
      UDScript, opts6___], 
     TagBox[AdjustmentBox[TagBox[z_, KetArgs, opts7___], opts8___], 
      UDScript, opts9___]], Braket, opts10___], form_] :=
  MakeExpression[
   RowBox[{"Braket", "[", 
     RowBox[{RowBox[{"{", w, "}"}], ",", 
       RowBox[{"{", 
         RowBox[{rowBoxToArgBoxes[x], ",", 
           RowBox[{"{", y, "}"}], ",", RowBox[{"{", z, "}"}]}], 
         "}"}]}], "]"}], form];

Braket /: MakeBoxes[Braket[{x__, {y__}, {}}, {z__}], form_] :=
  TagBox[RowBox[{SubscriptBox["\[InvisiblePrefixScriptBase]", 
      TagBox[AdjustmentBox[
        TagBox[inputToBoxes[y, form], BraArgs, style -> "BraKetArg"], 
        BoxBaselineShift -> 0], UDScript, style -> "LIndex"]], 
     RowBox[{"\[LeftAngleBracket]", 
       AdjustmentBox[
        RowBox[{TagBox[inputToBoxes[x, form], BraKetArgs, 
           style -> "BraKetArg"], "\[VerticalSeparator]", 
          TagBox[inputToBoxes[z, form], BraKetArgs, 
           style -> "BraKetArg"]}], BoxBaselineShift -> 0], 
       "\[RightAngleBracket]"}]}], Braket, style -> "BraKetWrapper"];

MakeExpression[
   TagBox[
    RowBox[{SubscriptBox["\[InvisiblePrefixScriptBase]", 
       TagBox[AdjustmentBox[TagBox[y_, BraArgs, opts1___], opts2___], 
        UDScript, opts3___]], 
      RowBox[{"\[LeftAngleBracket]", 
        AdjustmentBox[
         RowBox[{TagBox[x_, BraKetArgs, opts4___], 
           "\[VerticalSeparator]", TagBox[z_, BraKetArgs, opts5___]}],
          opts6___], "\[RightAngleBracket]"}]}], Braket, 
    style -> opts7___], form_] :=
  MakeExpression[
   RowBox[{"Braket", "[", 
     RowBox[{RowBox[{"{", 
         RowBox[{rowBoxToArgBoxes[x], ",", 
           RowBox[{"{", y, "}"}], ",", RowBox[{"{", "}"}]}], "}"}], 
       ",", RowBox[{"{", z, "}"}]}], "]"}], form];

Braket /: MakeBoxes[Braket[{x__, {}, {y__}}, {z__}], form_] :=
  TagBox[RowBox[{SuperscriptBox["\[InvisiblePrefixScriptBase]", 
      TagBox[AdjustmentBox[
        TagBox[inputToBoxes[y, form], BraArgs, style -> "BraKetArg"], 
        BoxBaselineShift -> 0], UDScript, style -> "LIndex"]], 
     RowBox[{"\[LeftAngleBracket]", 
       AdjustmentBox[
        RowBox[{TagBox[inputToBoxes[x, form], BraKetArgs, 
           style -> "BraKetArg"], "\[VerticalSeparator]", 
          TagBox[inputToBoxes[z, form], BraKetArgs, 
           style -> "BraKetArg"]}], BoxBaselineShift -> 0], 
       "\[RightAngleBracket]"}]}], Braket, style -> "BraKetWrapper"];

MakeExpression[
   TagBox[
    RowBox[{SuperscriptBox["\[InvisiblePrefixScriptBase]", 
       TagBox[AdjustmentBox[TagBox[y_, BraArgs, opts1___], opts2___], 
        UDScript, opts3___]], 
      RowBox[{"\[LeftAngleBracket]", 
        AdjustmentBox[
         RowBox[{TagBox[x_, BraKetArgs, opts4___], 
           "\[VerticalSeparator]", TagBox[z_, BraKetArgs, opts5___]}],
          opts6___], "\[RightAngleBracket]"}]}], Braket, 
    style -> opts7___], form_] :=
  MakeExpression[
   RowBox[{"Braket", "[", 
     RowBox[{RowBox[{"{", 
         RowBox[{rowBoxToArgBoxes[x], ",", 
           RowBox[{"{", "}"}], ",", RowBox[{"{", y, "}"}]}], "}"}], 
       ",", RowBox[{"{", z, "}"}]}], "]"}], form];

Braket /: MakeBoxes[Braket[{w__, {x__}, {y__}}, {z__}], form_] :=
  TagBox[RowBox[{SubsuperscriptBox["\[InvisiblePrefixScriptBase]", 
      TagBox[AdjustmentBox[
        TagBox[inputToBoxes[x, form], BraArgs, style -> "BraKetArg"], 
        BoxBaselineShift -> 0], UDScript, style -> "LIndex"], 
      TagBox[AdjustmentBox[
        TagBox[inputToBoxes[y, form], BraArgs, style -> "BraKetArg"], 
        BoxBaselineShift -> 0], UDScript, style -> "LIndex"]], 
     RowBox[{"\[LeftAngleBracket]", 
       AdjustmentBox[
        RowBox[{TagBox[inputToBoxes[w, form], BraKetArgs, 
           style -> "BraKetArg"], "\[VerticalSeparator]", 
          TagBox[inputToBoxes[z, form], BraKetArgs, 
           style -> "BraKetArg"]}], BoxBaselineShift -> 0], 
       "\[RightAngleBracket]"}]}], Braket, style -> "BraKetWrapper"];

MakeExpression[
   TagBox[
    RowBox[{SubsuperscriptBox["\[InvisiblePrefixScriptBase]", 
       TagBox[AdjustmentBox[TagBox[x_, BraArgs, opts1___], opts2___], 
        UDScript, opts3___], 
       TagBox[AdjustmentBox[TagBox[y_, BraArgs, opts4___], opts5___], 
        UDScript, opts6___]], 
      RowBox[{"\[LeftAngleBracket]", 
        AdjustmentBox[
         RowBox[{TagBox[w_, BraKetArgs, opts7___], 
           "\[VerticalSeparator]", TagBox[z_, BraKetArgs, opts8___]}],
          opts9___], "\[RightAngleBracket]"}]}], Braket, 
    style -> opts10___], form_] :=
  MakeExpression[
   RowBox[{"Braket", "[", 
     RowBox[{RowBox[{"{", 
         RowBox[{rowBoxToArgBoxes[w], ",", 
           RowBox[{"{", x, "}"}], ",", RowBox[{"{", y, "}"}]}], "}"}],
        ",", RowBox[{"{", z, "}"}]}], "]"}], form];

Braket /: 
  MakeBoxes[Braket[{w__, {x__}, {}}, {y__, {z__}, {}}], form_] :=
  TagBox[RowBox[{SubscriptBox["\[InvisiblePrefixScriptBase]", 
      TagBox[AdjustmentBox[
        TagBox[inputToBoxes[x, form], BraArgs, style -> "BraKetArg"], 
        BoxBaselineShift -> 0], UDScript, style -> "LIndex"]], 
     SubscriptBox[
      RowBox[{"\[LeftAngleBracket]", 
        AdjustmentBox[
         RowBox[{TagBox[inputToBoxes[w, form], BraKetArgs, 
            style -> "BraKetArg"], "\[VerticalSeparator]", 
           TagBox[inputToBoxes[y, form], BraKetArgs, 
            style -> "BraKetArg"]}], BoxBaselineShift -> 0], 
        "\[RightAngleBracket]"}], 
      TagBox[AdjustmentBox[
        TagBox[inputToBoxes[z, form], BraArgs, style -> "BraKetArg"], 
        BoxBaselineShift -> 0], UDScript, style -> "RIndex"]]}], 
   Braket, style -> "BraKetWrapper"];

MakeExpression[
   TagBox[
    RowBox[{SubscriptBox["\[InvisiblePrefixScriptBase]", 
       TagBox[AdjustmentBox[TagBox[x_, BraArgs, opts1___], opts2___], 
        UDScript, opts3___]], 
      SubscriptBox[
       RowBox[{"\[LeftAngleBracket]", 
         AdjustmentBox[
          RowBox[{TagBox[w_, BraKetArgs, opts4___], 
            "\[VerticalSeparator]", 
            TagBox[y_, BraKetArgs, opts5___]}], opts6___], 
         "\[RightAngleBracket]"}], 
       TagBox[AdjustmentBox[TagBox[z_, BraArgs, opts7___], opts8___], 
        UDScript, opts9___]]}], Braket, opts10___], form_] :=
  MakeExpression[
   RowBox[{"Braket", "[", 
     RowBox[{RowBox[{"{", 
         RowBox[{rowBoxToArgBoxes[w], ",", 
           RowBox[{"{", x, "}"}], ",", RowBox[{"{", "}"}]}], "}"}], 
       ",", RowBox[{"{", 
         RowBox[{rowBoxToArgBoxes[y], ",", 
           RowBox[{"{", z, "}"}], ",", RowBox[{"{", "}"}]}], "}"}]}], 
     "]"}], form];

Braket /: 
  MakeBoxes[Braket[{w__, {x__}, {}}, {y__, {}, {z__}}], form_] :=
  TagBox[RowBox[{SubscriptBox["\[InvisiblePrefixScriptBase]", 
      TagBox[AdjustmentBox[
        TagBox[inputToBoxes[x, form], BraArgs, style -> "BraKetArg"], 
        BoxBaselineShift -> 0], UDScript, style -> "LIndex"]], 
     SuperscriptBox[
      RowBox[{"\[LeftAngleBracket]", 
        AdjustmentBox[
         RowBox[{TagBox[inputToBoxes[w, form], BraKetArgs, 
            style -> "BraKetArg"], "\[VerticalSeparator]", 
           TagBox[inputToBoxes[y, form], BraKetArgs, 
            style -> "BraKetArg"]}], BoxBaselineShift -> 0], 
        "\[RightAngleBracket]"}], 
      TagBox[AdjustmentBox[
        TagBox[inputToBoxes[z, form], BraArgs, style -> "BraKetArg"], 
        BoxBaselineShift -> 0], UDScript, style -> "RIndex"]]}], 
   Braket, style -> "BraKetWrapper"];

MakeExpression[
   TagBox[
    RowBox[{SubscriptBox["\[InvisiblePrefixScriptBase]", 
       TagBox[AdjustmentBox[TagBox[x_, BraArgs, opts1___], opts2___], 
        UDScript, opts3___]], 
      SuperscriptBox[
       RowBox[{"\[LeftAngleBracket]", 
         AdjustmentBox[
          RowBox[{TagBox[w_, BraKetArgs, opts4___], 
            "\[VerticalSeparator]", 
            TagBox[y_, BraKetArgs, opts5___]}], opts6___], 
         "\[RightAngleBracket]"}], 
       TagBox[AdjustmentBox[TagBox[z_, BraArgs, opts7___], opts8___], 
        UDScript, opts9___]]}], Braket, opts10___], form_] :=
  MakeExpression[
   RowBox[{"Braket", "[", 
     RowBox[{RowBox[{"{", 
         RowBox[{rowBoxToArgBoxes[w], ",", 
           RowBox[{"{", x, "}"}], ",", RowBox[{"{", "}"}]}], "}"}], 
       ",", RowBox[{"{", 
         RowBox[{rowBoxToArgBoxes[y], ",", 
           RowBox[{"{", "}"}], ",", RowBox[{"{", z, "}"}]}], "}"}]}], 
     "]"}], form];

Braket /: 
  MakeBoxes[Braket[{w__, {}, {x__}}, {y__, {z__}, {}}], form_] :=
  TagBox[RowBox[{SuperscriptBox["\[InvisiblePrefixScriptBase]", 
      TagBox[AdjustmentBox[
        TagBox[inputToBoxes[x, form], BraArgs, style -> "BraKetArg"], 
        BoxBaselineShift -> 0], UDScript, style -> "LIndex"]], 
     SubscriptBox[
      RowBox[{"\[LeftAngleBracket]", 
        AdjustmentBox[
         RowBox[{TagBox[inputToBoxes[w, form], BraKetArgs, 
            style -> "BraKetArg"], "\[VerticalSeparator]", 
           TagBox[inputToBoxes[y, form], BraKetArgs, 
            style -> "BraKetArg"]}], BoxBaselineShift -> 0], 
        "\[RightAngleBracket]"}], 
      TagBox[AdjustmentBox[
        TagBox[inputToBoxes[z, form], BraArgs, style -> "BraKetArg"], 
        BoxBaselineShift -> 0], UDScript, style -> "RIndex"]]}], 
   Braket, style -> "BraKetWrapper"];

MakeExpression[
   TagBox[
    RowBox[{SuperscriptBox["\[InvisiblePrefixScriptBase]", 
       TagBox[AdjustmentBox[TagBox[x_, BraArgs, opts1___], opts2___], 
        UDScript, opts3___]], 
      SubscriptBox[
       RowBox[{"\[LeftAngleBracket]", 
         AdjustmentBox[
          RowBox[{TagBox[w_, BraKetArgs, opts4___], 
            "\[VerticalSeparator]", 
            TagBox[y_, BraKetArgs, opts5___]}], opts6___], 
         "\[RightAngleBracket]"}], 
       TagBox[AdjustmentBox[TagBox[z_, BraArgs, opts7___], opts8___], 
        UDScript, opts9___]]}], Braket, opts10___], form_] :=
  MakeExpression[
   RowBox[{"Braket", "[", 
     RowBox[{RowBox[{"{", 
         RowBox[{rowBoxToArgBoxes[w], ",", 
           RowBox[{"{", "}"}], ",", RowBox[{"{", x, "}"}]}], "}"}], 
       ",", RowBox[{"{", 
         RowBox[{rowBoxToArgBoxes[y], ",", 
           RowBox[{"{", z, "}"}], ",", RowBox[{"{", "}"}]}], "}"}]}], 
     "]"}], form];

Braket /: 
  MakeBoxes[Braket[{w__, {}, {x__}}, {y__, {}, {z__}}], form_] :=
  TagBox[RowBox[{SuperscriptBox["\[InvisiblePrefixScriptBase]", 
      TagBox[AdjustmentBox[
        TagBox[inputToBoxes[x, form], BraArgs, style -> "BraKetArg"], 
        BoxBaselineShift -> 0], UDScript, style -> "LIndex"]], 
     SuperscriptBox[
      RowBox[{"\[LeftAngleBracket]", 
        AdjustmentBox[
         RowBox[{TagBox[inputToBoxes[w, form], BraKetArgs, 
            style -> "BraKetArg"], "\[VerticalSeparator]", 
           TagBox[inputToBoxes[y, form], BraKetArgs, 
            style -> "BraKetArg"]}], BoxBaselineShift -> 0], 
        "\[RightAngleBracket]"}], 
      TagBox[AdjustmentBox[
        TagBox[inputToBoxes[z, form], BraArgs, style -> "BraKetArg"], 
        BoxBaselineShift -> 0], UDScript, style -> "RIndex"]]}], 
   Braket, style -> "BraKetWrapper"];

MakeExpression[
   TagBox[
    RowBox[{SuperscriptBox["\[InvisiblePrefixScriptBase]", 
       TagBox[AdjustmentBox[TagBox[x_, BraArgs, opts1___], opts2___], 
        UDScript, opts3___]], 
      SuperscriptBox[
       RowBox[{"\[LeftAngleBracket]", 
         AdjustmentBox[
          RowBox[{TagBox[w_, BraKetArgs, opts4___], 
            "\[VerticalSeparator]", 
            TagBox[y_, BraKetArgs, opts5___]}], opts6___], 
         "\[RightAngleBracket]"}], 
       TagBox[AdjustmentBox[TagBox[z_, BraArgs, opts7___], opts8___], 
        UDScript, opts9___]]}], Braket, opts10___], form_] :=
  MakeExpression[
   RowBox[{"Braket", "[", 
     RowBox[{RowBox[{"{", 
         RowBox[{rowBoxToArgBoxes[w], ",", 
           RowBox[{"{", "}"}], ",", RowBox[{"{", x, "}"}]}], "}"}], 
       ",", RowBox[{"{", 
         RowBox[{rowBoxToArgBoxes[y], ",", 
           RowBox[{"{", "}"}], ",", RowBox[{"{", z, "}"}]}], "}"}]}], 
     "]"}], form];

Braket /: 
  MakeBoxes[Braket[{v__, {w__}, {}}, {x__, {y__}, {z__}}], form_] :=
  TagBox[RowBox[{SubscriptBox["\[InvisiblePrefixScriptBase]", 
      TagBox[AdjustmentBox[
        TagBox[inputToBoxes[w, form], BraArgs, style -> "BraKetArg"], 
        BoxBaselineShift -> 0], UDScript, style -> "LIndex"]], 
     SubsuperscriptBox[
      RowBox[{"\[LeftAngleBracket]", 
        AdjustmentBox[
         RowBox[{TagBox[inputToBoxes[v, form], BraKetArgs, 
            style -> "BraKetArg"], "\[VerticalSeparator]", 
           TagBox[inputToBoxes[x, form], BraKetArgs, 
            style -> "BraKetArg"]}], BoxBaselineShift -> 0], 
        "\[RightAngleBracket]"}], 
      TagBox[AdjustmentBox[
        TagBox[inputToBoxes[y, form], BraArgs, style -> "BraKetArg"], 
        BoxBaselineShift -> 0], UDScript, style -> "RIndex"], 
      TagBox[AdjustmentBox[
        TagBox[inputToBoxes[z, form], BraArgs, style -> "BraKetArg"], 
        BoxBaselineShift -> 0], UDScript, style -> "RIndex"]]}], 
   Braket, style -> "BraKetWrapper"];

MakeExpression[
   TagBox[
    RowBox[{SubscriptBox["\[InvisiblePrefixScriptBase]", 
       TagBox[AdjustmentBox[TagBox[w_, BraArgs, opts1___], opts2___], 
        UDScript, opts3___]], 
      SubsuperscriptBox[
       RowBox[{"\[LeftAngleBracket]", 
         AdjustmentBox[
          RowBox[{TagBox[v_, BraKetArgs, opts4___], 
            "\[VerticalSeparator]", 
            TagBox[x_, BraKetArgs, opts5___]}], opts6___], 
         "\[RightAngleBracket]"}], 
       TagBox[AdjustmentBox[TagBox[y_, BraArgs, opts7___], opts8___], 
        UDScript, opts9___], 
       TagBox[AdjustmentBox[TagBox[z_, BraArgs, opts10___], 
         opts11___], UDScript, opts12___]]}], Braket, opts13___], 
   form_] :=
  MakeExpression[
   RowBox[{"Braket", "[", 
     RowBox[{RowBox[{"{", 
         RowBox[{rowBoxToArgBoxes[v], ",", 
           RowBox[{"{", w, "}"}], ",", RowBox[{"{", "}"}]}], "}"}], 
       ",", RowBox[{"{", 
         RowBox[{rowBoxToArgBoxes[x], ",", 
           RowBox[{"{", y, "}"}], ",", RowBox[{"{", z, "}"}]}], 
         "}"}]}], "]"}], form];

Braket /: 
  MakeBoxes[Braket[{v__, {}, {w__}}, {x__, {y__}, {z__}}], form_] :=
  TagBox[RowBox[{SuperscriptBox["\[InvisiblePrefixScriptBase]", 
      TagBox[
       AdjustmentBox[
        TagBox[inputToBoxes[w, form], BraArgs, style -> "BraKetArg"], 
        BoxBaselineShift -> 0], UDScript, style -> "LIndex"]], 
     SubsuperscriptBox[
      RowBox[{"\[LeftAngleBracket]", 
        AdjustmentBox[
         RowBox[{TagBox[inputToBoxes[v, form], BraKetArgs, 
            style -> "BraKetArg"], "\[VerticalSeparator]", 
           TagBox[inputToBoxes[x, form], BraKetArgs, 
            style -> "BraKetArg"]}], BoxBaselineShift -> 0], 
        "\[RightAngleBracket]"}], 
      TagBox[AdjustmentBox[
        TagBox[inputToBoxes[y, form], BraArgs, style -> "BraKetArg"], 
        BoxBaselineShift -> 0], UDScript, style -> "RIndex"], 
      TagBox[AdjustmentBox[
        TagBox[inputToBoxes[z, form], BraArgs, style -> "BraKetArg"], 
        BoxBaselineShift -> 0], UDScript, style -> "RIndex"]]}], 
   Braket, style -> "BraKetWrapper"];

MakeExpression[
   TagBox[
    RowBox[{SuperscriptBox["\[InvisiblePrefixScriptBase]", 
       TagBox[AdjustmentBox[TagBox[w_, BraArgs, opts1___], opts2___], 
        UDScript, opts3___]], 
      SubsuperscriptBox[
       RowBox[{"\[LeftAngleBracket]", 
         AdjustmentBox[
          RowBox[{TagBox[v_, BraKetArgs, opts4___], 
            "\[VerticalSeparator]", 
            TagBox[x_, BraKetArgs, opts5___]}], opts6___], 
         "\[RightAngleBracket]"}], 
       TagBox[AdjustmentBox[TagBox[y_, BraArgs, opts7___], opts8___], 
        UDScript, opts9___], 
       TagBox[AdjustmentBox[TagBox[z_, BraArgs, opts10___], 
         opts11___], UDScript, opts12___]]}], Braket, opts13___], 
   form_] :=
  MakeExpression[
   RowBox[{"Braket", "[", 
     RowBox[{RowBox[{"{", 
         RowBox[{rowBoxToArgBoxes[v], ",", 
           RowBox[{"{", "}"}], ",", RowBox[{"{", w, "}"}]}], "}"}], 
       ",", RowBox[{"{", 
         RowBox[{rowBoxToArgBoxes[x], ",", 
           RowBox[{"{", y, "}"}], ",", RowBox[{"{", z, "}"}]}], 
         "}"}]}], "]"}], form];

Braket /: 
  MakeBoxes[Braket[{v__, {w__}, {x__}}, {y__, {z__}, {}}], form_] :=
  TagBox[RowBox[{SubsuperscriptBox["\[InvisiblePrefixScriptBase]", 
      TagBox[AdjustmentBox[
        TagBox[inputToBoxes[w, form], BraArgs, style -> "BraKetArg"], 
        BoxBaselineShift -> 0], UDScript, style -> "LIndex"], 
      TagBox[AdjustmentBox[
        TagBox[inputToBoxes[x, form], BraArgs, 
         style -> "BraKetArg"], BoxBaselineShift -> 0], UDScript, 
       style -> "LIndex"]], 
     SubscriptBox[
      RowBox[{"\[LeftAngleBracket]", 
        AdjustmentBox[
         RowBox[{TagBox[inputToBoxes[v, form], BraKetArgs, 
            style -> "BraKetArg"], "\[VerticalSeparator]", 
           TagBox[inputToBoxes[y, form], BraKetArgs, 
            style -> "BraKetArg"]}], BoxBaselineShift -> 0], 
        "\[RightAngleBracket]"}], 
      TagBox[AdjustmentBox[
        TagBox[inputToBoxes[z, form], BraArgs, style -> "BraKetArg"], 
        BoxBaselineShift -> 0], UDScript, style -> "RIndex"]]}], 
   Braket, style -> "BraKetWrapper"];

MakeExpression[
   TagBox[
    RowBox[{SubsuperscriptBox["\[InvisiblePrefixScriptBase]", 
       TagBox[AdjustmentBox[TagBox[w_, BraArgs, opts1___], opts2___], 
        UDScript, opts3___], 
       TagBox[AdjustmentBox[TagBox[x_, BraArgs, opts4___], opts5___], 
        UDScript, opts6___]], 
      SubscriptBox[
       RowBox[{"\[LeftAngleBracket]", 
         AdjustmentBox[
          RowBox[{TagBox[v_, BraKetArgs, opts7___], 
            "\[VerticalSeparator]", 
            TagBox[y_, BraKetArgs, opts8___]}], opts9___], 
         "\[RightAngleBracket]"}], 
       TagBox[AdjustmentBox[TagBox[z_, BraArgs, opts10___], 
         opts11___], UDScript, opts12___]]}], Braket, opts13___], 
   form_] :=
  MakeExpression[
   RowBox[{"Braket", "[", 
     RowBox[{RowBox[{"{", 
         RowBox[{rowBoxToArgBoxes[v], ",", 
           RowBox[{"{", w, "}"}], ",", RowBox[{"{", x, "}"}]}], "}"}],
        ",", RowBox[{"{", 
         RowBox[{rowBoxToArgBoxes[y], ",", 
           RowBox[{"{", z, "}"}], ",", RowBox[{"{", "}"}]}], "}"}]}], 
     "]"}], form];

Braket /: 
  MakeBoxes[Braket[{v__, {w__}, {x__}}, {y__, {}, {z__}}], form_] :=
  TagBox[RowBox[{SubsuperscriptBox["\[InvisiblePrefixScriptBase]", 
      TagBox[AdjustmentBox[
        TagBox[inputToBoxes[w, form], BraArgs, style -> "BraKetArg"], 
        BoxBaselineShift -> 0], UDScript, style -> "LIndex"], 
      TagBox[AdjustmentBox[
        TagBox[inputToBoxes[x, form], BraArgs, 
         style -> "BraKetArg"], BoxBaselineShift -> 0], UDScript, 
       style -> "LIndex"]], 
     SuperscriptBox[
      RowBox[{"\[LeftAngleBracket]", 
        AdjustmentBox[
         RowBox[{TagBox[inputToBoxes[v, form], BraKetArgs, 
            style -> "BraKetArg"], "\[VerticalSeparator]", 
           TagBox[inputToBoxes[y, form], BraKetArgs, 
            style -> "BraKetArg"]}], BoxBaselineShift -> 0], 
        "\[RightAngleBracket]"}], 
      TagBox[AdjustmentBox[
        TagBox[inputToBoxes[z, form], BraArgs, style -> "BraKetArg"], 
        BoxBaselineShift -> 0], UDScript, style -> "RIndex"]]}], 
   Braket, style -> "BraKetWrapper"];

MakeExpression[
   TagBox[
    RowBox[{SubsuperscriptBox["\[InvisiblePrefixScriptBase]", 
       TagBox[AdjustmentBox[TagBox[w_, BraArgs, opts1___], opts2___], 
        UDScript, opts3___], 
       TagBox[AdjustmentBox[TagBox[x_, BraArgs, opts4___], opts5___], 
        UDScript, opts6___]], 
      SuperscriptBox[
       RowBox[{"\[LeftAngleBracket]", 
         AdjustmentBox[
          RowBox[{TagBox[v_, BraKetArgs, opts7___], 
            "\[VerticalSeparator]", 
            TagBox[y_, BraKetArgs, opts8___]}], opts9___], 
         "\[RightAngleBracket]"}], 
       TagBox[AdjustmentBox[TagBox[z_, BraArgs, opts10___], 
         opts11___], UDScript, opts12___]]}], Braket, opts13___], 
   form_] :=
  MakeExpression[
   RowBox[{"Braket", "[", 
     RowBox[{RowBox[{"{", 
         RowBox[{rowBoxToArgBoxes[v], ",", 
           RowBox[{"{", w, "}"}], ",", RowBox[{"{", x, "}"}]}], "}"}],
        ",", RowBox[{"{", 
         RowBox[{rowBoxToArgBoxes[y], ",", 
           RowBox[{"{", "}"}], ",", RowBox[{"{", z, "}"}]}], "}"}]}], 
     "]"}], form];

Braket /: 
  MakeBoxes[Braket[{u__, {v__}, {w__}}, {x__, {y__}, {z__}}], form_] :=
   TagBox[RowBox[{SubsuperscriptBox["\[InvisiblePrefixScriptBase]", 
      TagBox[AdjustmentBox[
        TagBox[inputToBoxes[v, form], BraArgs, style -> "BraKetArg"], 
        BoxBaselineShift -> 0.`], UDScript, style -> "LIndex"], 
      TagBox[
       AdjustmentBox[
        TagBox[inputToBoxes[w, form], BraArgs, style -> "BraKetArg"], 
        BoxBaselineShift -> 0.`], UDScript, style -> "LIndex"]], 
     SubsuperscriptBox[
      RowBox[{"\[LeftAngleBracket]", 
        AdjustmentBox[
         RowBox[{TagBox[inputToBoxes[u, form], BraKetArgs, 
            style -> "BraKetArg"], "\[VerticalSeparator]", 
           TagBox[inputToBoxes[x, form], BraKetArgs, 
            style -> "BraKetArg"]}], BoxBaselineShift -> 0.`], 
        "\[RightAngleBracket]"}], 
      TagBox[AdjustmentBox[
        TagBox[inputToBoxes[y, form], BraArgs, style -> "BraKetArg"], 
        BoxBaselineShift -> 0.`], UDScript, style -> "RIndex"], 
      TagBox[AdjustmentBox[
        TagBox[inputToBoxes[z, form], BraArgs, style -> "BraKetArg"], 
        BoxBaselineShift -> 0.`], UDScript, style -> "RIndex"]]}], 
   Braket, style -> "BraKetWrapper"];

MakeExpression[
   TagBox[
    RowBox[{SubsuperscriptBox["\[InvisiblePrefixScriptBase]", 
       TagBox[AdjustmentBox[TagBox[v_, BraArgs, opts1___], opts2___], 
        UDScript, opts3___], 
       TagBox[AdjustmentBox[TagBox[w_, BraArgs, opts4___], opts5___], 
        UDScript, opts6___]], 
      SubsuperscriptBox[
       RowBox[{"\[LeftAngleBracket]", 
         AdjustmentBox[
          RowBox[{TagBox[u_, BraKetArgs, opts7___], 
            "\[VerticalSeparator]", 
            TagBox[x_, BraKetArgs, opts8___]}], opts9___], 
         "\[RightAngleBracket]"}], 
       TagBox[AdjustmentBox[TagBox[y_, BraArgs, opts10___], 
         opts11___], UDScript, opts12___], 
       TagBox[AdjustmentBox[TagBox[z_, BraArgs, opts13___], 
         opts14___], UDScript, opts15___]]}], Braket, opts16___], 
   form_] :=
  MakeExpression[
   RowBox[{"Braket", "[", 
     RowBox[{RowBox[{"{", 
         RowBox[{rowBoxToArgBoxes[u], ",", 
           RowBox[{"{", v, "}"}], ",", RowBox[{"{", w, "}"}]}], "}"}],
        ",", RowBox[{"{", 
         RowBox[{rowBoxToArgBoxes[x], ",", 
           RowBox[{"{", y, "}"}], ",", RowBox[{"{", z, "}"}]}], 
         "}"}]}], "]"}], form];


(* ::Subsubsection:: *)
(*Define Commutator Notation*)


Commutator /: MakeBoxes[Commutator[x_, y_], form_] :=
  SubscriptBox[
   RowBox[{"[", RowBox[{MakeBoxes[x, form], ",", MakeBoxes[y, form]}],
      "]"}], "-"];

MakeExpression[
   SubscriptBox[RowBox[{"[", RowBox[{x_, ",", y_}], "]"}], "-"], 
   form_] :=
  MakeExpression[
   RowBox[{"Commutator", "[", RowBox[{x, ",", y}], "]"}], form];


(* ::Subsubsection:: *)
(*Define Conjugate Notation*)


SuperStar = Conjugate;
SuperDagger = Hermitian;

Conjugate /: MakeBoxes[Conjugate[x_?AtomQ], form_] :=
  SuperscriptBox[MakeBoxes[x, form], "*"];

Conjugate /: MakeBoxes[Conjugate[x_?(! AtomQ[#] &)], form_] :=
  SuperscriptBox[RowBox[{"(", MakeBoxes[x], ")"}], "*"];


(* ::Subsubsection:: *)
(*Define NonCommutativeMultiply Notation*)


CenterDot = NonCommutativeMultiply;


NonCommutativeMultiply /: 
 MakeBoxes[NonCommutativeMultiply[x__], form_] := 
 RowBox[Drop[
   Flatten[Map[If[Head[#] =!= Times, {"\[CenterDot]", MakeBoxes[#, form]},
	 {"\[CenterDot]", "(", MakeBoxes[#, form],")"}] &, {x}]], 1]] /; Length[{x}]>1


(* ::Subsubsection:: *)
(*Define TensorProduct Notation*)


CircleTimes = TensorProduct;


TensorProduct/:
MakeBoxes[TensorProduct[x__], form_] :=
 RowBox[Drop[
   Flatten[Map[{"\[CircleTimes]", MakeBoxes[#, form]} &, {x}]], 1]] /; Length[{x}]>1


(* ::Subsubsection:: *)
(*Style for Quantum Objects*)


QAStyles :=
  Cell[CellGroupData[{
     Cell["Bra Ket Styles", "Section"], 
     Cell[" The cells below define certain styles needed for defining Bras&Kets. These styles serve to give the correct structural properties to Bras & Kets ", "Text"],
      Cell[StyleData["BraKetArg"], SpanMinSize -> Automatic, 
      StyleMenuListing -> None,  
      AdjustmentBoxOptions -> {BoxMargins -> {{0, 0}, {0, 0}}}, 
      TagBoxOptions -> {Editable -> True, Selectable -> True}], 
     Cell[StyleData["BraWrapper"], AutoStyleOptions->{"HighlightSyntaxErrors"->False}, 
      ScriptBaselineShifts -> {0.6, 0.9},
       SpanMinSize -> 1.5, StyleMenuListing -> None, 
      AdjustmentBoxOptions -> {BoxBaselineShift->-0.01, 
							   BoxMargins -> {{0.2, -0.13}, {0, 0}}}, 
      TagBoxOptions -> {Editable -> False}], 
     Cell[StyleData["KetWrapper"], AutoStyleOptions->{"HighlightSyntaxErrors"->False}, 
      ScriptBaselineShifts -> {0.6, 0.9}, SpanMinSize -> 1.5, 
      StyleMenuListing -> None,  
      AdjustmentBoxOptions -> {BoxBaselineShift->-0.01, 
							   BoxMargins -> {{-0.2, 0.2}, {0, 0}}}, 
      TagBoxOptions -> {Editable -> False}], 
     Cell[StyleData["BraKetWrapper"], 
      ScriptBaselineShifts -> {0.6, 0.9}, SpanMinSize -> 1.5, 
      StyleMenuListing -> None, 
      AdjustmentBoxOptions -> {BoxMargins -> {{0.25, 0.25}, {0, 0}}}, 
      TagBoxOptions -> {Editable -> False}], 
     Cell[StyleData["RIndex"], StyleMenuListing -> None,  
      AdjustmentBoxOptions -> {BoxMargins -> {{-0.16, 0}, {0, 0}}}, 
      TagBoxOptions -> {Editable -> False}], 
     Cell[StyleData["LIndex"], StyleMenuListing -> None,  
      AdjustmentBoxOptions -> {BoxMargins -> {{0, -0.1}, {0, 0}}}, 
      TagBoxOptions -> {Editable -> False}]}, Closed]]; 

styleSheetInUse = 
 StyleDefinitions /. AbsoluteOptions[nb1, StyleDefinitions];

If[Head[styleSheetInUse] === Notebook,
 If[Or @@ (FreeQ[{styleSheetInUse}, #] &) /@ 
    Union @ Cases[QAStyles, StyleData[_String], {0, \[Infinity]}],
  newNBStyles = 
   styleSheetInUse /. 
    Notebook[cells_List, opts___] :> 
     Notebook[Append[cells, QAStyles], opts];
  SetOptions[nb1, StyleDefinitions -> newNBStyles]
  ]
 ];

If[Or[Head[styleSheetInUse] === FrontEnd`FileName, 
  Head[styleSheetInUse] === String],
 SetOptions[nb1, StyleDefinitions -> Notebook[{
     Cell[StyleData[StyleDefinitions -> styleSheetInUse]],
     QAStyles},
    StyleDefinitions -> "PrivateStylesheetFormatting.nb"]]
 ];



(* ::Subsubsection:: *)
(*Set New Alias For Notebook*)

qaInputAlias = 
{
	"ket" -> 
     TagBox[RowBox[{AdjustmentBox["\[VerticalSeparator]", 
         BoxMargins -> {{-0.2, 0.}, {0., 0.}},
         BoxBaselineShift -> -0.1], 
        AdjustmentBox[
         TagBox["\[SelectionPlaceholder]", KetArgs, 
          style -> "BraKetArg"], BoxBaselineShift -> 0], 
        "\[RightAngleBracket]"}], Ket, style -> "KetWrapper"],
    "bra" ->
     TagBox[
      RowBox[{"\[LeftAngleBracket]", 
        AdjustmentBox[
         TagBox["\[SelectionPlaceholder]", BraArgs, 
          style -> "BraKetArg"], BoxBaselineShift -> 0], 
        AdjustmentBox["\[VerticalSeparator]", 
         BoxMargins -> {{0, -0.2}, {0, 0}},
         BoxBaselineShift -> -0.1]}], Bra, 
      style -> "BraWrapper"],
    "braket" -> 
     TagBox[
      RowBox[{"\[LeftAngleBracket]", 
        AdjustmentBox[
         RowBox[{TagBox["\[SelectionPlaceholder]", BraKetArgs, 
            style -> "BraKetArg"], "\[VerticalSeparator]", 
           TagBox["\[SelectionPlaceholder]", BraKetArgs, 
            style -> "BraKetArg"]}], BoxBaselineShift -> 0.`], 
        "\[RightAngleBracket]"}], Braket, style -> "BraKetWrapper"],
    "ketd" ->
     TagBox[
      SubscriptBox[
       RowBox[
        {AdjustmentBox["\[VerticalSeparator]", 
          BoxMargins -> {{-0.2, 0}, {0, 0}},
          BoxBaselineShift -> -0.1], 
         AdjustmentBox[
          TagBox["\[SelectionPlaceholder]", KetArgs, 
           style -> "BraKetArg"], BoxBaselineShift -> 0], 
         "\[RightAngleBracket]"}], 
       TagBox[AdjustmentBox[
         TagBox["\[SelectionPlaceholder]", KetArgs, 
          style -> "BraKetArg"], BoxBaselineShift -> 0], UDScript, 
        style -> "RIndex"]], Ket, style -> "KetWrapper"],
    "ketu" ->
     TagBox[
      SuperscriptBox[
       RowBox[{AdjustmentBox["\[VerticalSeparator]", 
          BoxMargins -> {{-0.2, 0}, {0, 0}},
          BoxBaselineShift -> -0.1], 
         AdjustmentBox[
          TagBox["\[SelectionPlaceholder]", KetArgs, 
           style -> "BraKetArg"], BoxBaselineShift -> 0], 
         "\[RightAngleBracket]"}], 
       TagBox[AdjustmentBox[
         TagBox["\[SelectionPlaceholder]", KetArgs, 
          style -> "BraKetArg"], BoxBaselineShift -> 0], UDScript, 
        style -> "RIndex"]], Ket, style -> "KetWrapper"],
    "ketdu" ->
     TagBox[
      SubsuperscriptBox[
       RowBox[{AdjustmentBox["\[VerticalSeparator]", 
          BoxMargins -> {{-0.2, 0}, {0, 0}},
          BoxBaselineShift -> -0.1], 
         AdjustmentBox[
          TagBox["\[SelectionPlaceholder]", KetArgs, 
           style -> "BraKetArg"], BoxBaselineShift -> 0], 
         "\[RightAngleBracket]"}], 
       TagBox[AdjustmentBox[
         TagBox["\[SelectionPlaceholder]", KetArgs, 
          style -> "BraKetArg"], BoxBaselineShift -> 0], UDScript, 
        style -> "RIndex"], 
       TagBox[AdjustmentBox[
         TagBox["\[SelectionPlaceholder]", KetArgs, 
          style -> "BraKetArg"], BoxBaselineShift -> 0], UDScript, 
        style -> "RIndex"]], Ket, style -> "KetWrapper"],
    "dbra" ->
     TagBox[
      RowBox[{SubscriptBox["\[InvisiblePrefixScriptBase]", 
         TagBox[AdjustmentBox[
           TagBox["\[SelectionPlaceholder]", BraArgs, 
            style -> "BraKetArg"], BoxBaselineShift -> 0], UDScript, 
          style -> "LIndex"]], 
        RowBox[{"\[LeftAngleBracket]", 
          AdjustmentBox[
           TagBox["\[SelectionPlaceholder]", BraArgs, 
            style -> "BraKetArg"], BoxBaselineShift -> 0], 
          AdjustmentBox["\[VerticalSeparator]", 
           BoxMargins -> {{0, -0.2}, {0, 0}},
           BoxBaselineShift -> -0.1]}]}], Bra, 
      style -> "BraWrapper"],
    "ubra" ->
     TagBox[
      RowBox[{SuperscriptBox["\[InvisiblePrefixScriptBase]", 
         TagBox[AdjustmentBox[
           TagBox["\[SelectionPlaceholder]", BraArgs, 
            style -> "BraKetArg"], BoxBaselineShift -> 0], UDScript, 
          style -> "LIndex"]], 
        RowBox[{"\[LeftAngleBracket]", 
          AdjustmentBox[
           TagBox["\[SelectionPlaceholder]", BraArgs, 
            style -> "BraKetArg"], BoxBaselineShift -> 0], 
          AdjustmentBox["\[VerticalSeparator]", 
           BoxMargins -> {{0, -0.2}, {0, 0}},
           BoxBaselineShift -> -0.1]}]}], Bra, 
      style -> "BraWrapper"],
    "dubra" -> 
     TagBox[RowBox[{SubsuperscriptBox["\[InvisiblePrefixScriptBase]", 
         TagBox[AdjustmentBox[
           TagBox["\[SelectionPlaceholder]", BraArgs, 
            style -> "BraKetArg"], BoxBaselineShift -> 0], UDScript, 
          style -> "LIndex"], 
         TagBox[AdjustmentBox[
           TagBox["\[SelectionPlaceholder]", BraArgs, 
            style -> "BraKetArg"], BoxBaselineShift -> 0], UDScript, 
          style -> "LIndex"]], 
        RowBox[{"\[LeftAngleBracket]", 
          AdjustmentBox[
           TagBox["\[SelectionPlaceholder]", BraArgs, 
            style -> "BraKetArg"], BoxBaselineShift -> 0], 
          AdjustmentBox["\[VerticalSeparator]", 
           BoxMargins -> {{0, -0.2}, {0, 0}},
           BoxBaselineShift -> -0.1]}]}], Bra, 
      style -> "BraWrapper"],
    "braketd" ->
     TagBox[
      SubscriptBox[
       RowBox[{"\[LeftAngleBracket]", 
         AdjustmentBox[
          RowBox[{TagBox["\[SelectionPlaceholder]", BraKetArgs, 
             style -> "BraKetArg"], "\[VerticalSeparator]", 
            TagBox["\[SelectionPlaceholder]", BraKetArgs, 
             style -> "BraKetArg"]}], BoxBaselineShift -> 0], 
         "\[RightAngleBracket]"}], 
       TagBox[AdjustmentBox[
         TagBox["\[SelectionPlaceholder]", KetArgs, 
          style -> "BraKetArg"], BoxBaselineShift -> 0], UDScript, 
        style -> "RIndex"]], Braket, style -> "BraKetWrapper"],
    "braketu" ->
     TagBox[
      SuperscriptBox[
       RowBox[{"\[LeftAngleBracket]", 
         AdjustmentBox[
          
          RowBox[{TagBox["\[SelectionPlaceholder]", BraKetArgs, 
             style -> "BraKetArg"], "\[VerticalSeparator]", 
            TagBox["\[SelectionPlaceholder]", BraKetArgs, 
             style -> "BraKetArg"]}], BoxBaselineShift -> 0], 
         "\[RightAngleBracket]"}], 
       TagBox[AdjustmentBox[
         TagBox["\[SelectionPlaceholder]", KetArgs, 
          style -> "BraKetArg"], BoxBaselineShift -> 0], UDScript, 
        style -> "RIndex"]], Braket, style -> "BraKetWrapper"],
    "braketdu" ->
     TagBox[
      SubsuperscriptBox[
       RowBox[{"\[LeftAngleBracket]", 
         AdjustmentBox[
          RowBox[{TagBox["\[SelectionPlaceholder]", BraKetArgs, 
             style -> "BraKetArg"], "\[VerticalSeparator]", 
            TagBox["\[SelectionPlaceholder]", BraKetArgs, 
             style -> "BraKetArg"]}], BoxBaselineShift -> 0], 
         "\[RightAngleBracket]"}], 
       TagBox[AdjustmentBox[
         TagBox["\[SelectionPlaceholder]", KetArgs, 
          style -> "BraKetArg"], BoxBaselineShift -> 0], UDScript, 
        style -> "RIndex"], 
       TagBox[AdjustmentBox[
         TagBox["\[SelectionPlaceholder]", KetArgs, 
          style -> "BraKetArg"], BoxBaselineShift -> 0], UDScript, 
        style -> "RIndex"]], Braket, style -> "BraKetWrapper"],
    "dbraket" ->
     TagBox[
      RowBox[{SubscriptBox["\[InvisiblePrefixScriptBase]", 
         TagBox[AdjustmentBox[
           TagBox["\[SelectionPlaceholder]", BraArgs, 
            style -> "BraKetArg"], BoxBaselineShift -> 0], UDScript, 
          style -> "LIndex"]], 
        RowBox[{"\[LeftAngleBracket]", 
          AdjustmentBox[
           RowBox[{TagBox["\[SelectionPlaceholder]", BraKetArgs, 
              style -> "BraKetArg"], "\[VerticalSeparator]", 
             TagBox["\[SelectionPlaceholder]", BraKetArgs, 
              style -> "BraKetArg"]}], BoxBaselineShift -> 0], 
          "\[RightAngleBracket]"}]}], Braket, 
      style -> "BraKetWrapper"],
    "ubraket" ->
     TagBox[
      RowBox[{SuperscriptBox["\[InvisiblePrefixScriptBase]", 
         TagBox[AdjustmentBox[
           TagBox["\[SelectionPlaceholder]", BraArgs, 
            style -> "BraKetArg"], BoxBaselineShift -> 0], UDScript, 
          style -> "LIndex"]], 
        RowBox[{"\[LeftAngleBracket]",      
          AdjustmentBox[
           RowBox[{TagBox["\[SelectionPlaceholder]", BraKetArgs, 
              style -> "BraKetArg"], "\[VerticalSeparator]", 
             TagBox["\[SelectionPlaceholder]", BraKetArgs, 
              style -> "BraKetArg"]}], BoxBaselineShift -> 0], 
          "\[RightAngleBracket]"}]}], Braket, 
      style -> "BraKetWrapper"],
    "dubraket" -> 
     TagBox[RowBox[{SubsuperscriptBox["\[InvisiblePrefixScriptBase]", 
         TagBox[AdjustmentBox[
           TagBox["\[SelectionPlaceholder]", BraArgs, 
            style -> "BraKetArg"], BoxBaselineShift -> 0], UDScript, 
          style -> "LIndex"], 
         TagBox[AdjustmentBox[
           TagBox["\[SelectionPlaceholder]", BraArgs, 
            style -> "BraKetArg"], BoxBaselineShift -> 0], UDScript, 
          style -> "LIndex"]], 
        RowBox[{"\[LeftAngleBracket]", 
          AdjustmentBox[
           RowBox[{TagBox["\[SelectionPlaceholder]", BraKetArgs, 
              style -> "BraKetArg"], "\[VerticalSeparator]", 
             TagBox["\[SelectionPlaceholder]", BraKetArgs, 
              style -> "BraKetArg"]}], BoxBaselineShift -> 0], 
          "\[RightAngleBracket]"}]}], Braket, 
      style -> "BraKetWrapper"],
    "dbraketd" ->
     TagBox[
      RowBox[{SubscriptBox["\[InvisiblePrefixScriptBase]", 
         TagBox[AdjustmentBox[
           TagBox["\[SelectionPlaceholder]", BraArgs, 
            style -> "BraKetArg"], BoxBaselineShift -> 0], UDScript, 
          style -> "LIndex"]], 
        SubscriptBox[
         RowBox[{"\[LeftAngleBracket]", 
           AdjustmentBox[
            RowBox[{TagBox["\[SelectionPlaceholder]", BraKetArgs, 
               style -> "BraKetArg"], "\[VerticalSeparator]", 
              TagBox["\[SelectionPlaceholder]", BraKetArgs, 
               style -> "BraKetArg"]}], BoxBaselineShift -> 0], 
           "\[RightAngleBracket]"}], 
         TagBox[AdjustmentBox[
           TagBox["\[SelectionPlaceholder]", BraArgs, 
            style -> "BraKetArg"], BoxBaselineShift -> 0], UDScript, 
          style -> "RIndex"]]}], Braket, style -> "BraKetWrapper"],
    "dbraketu" ->
     TagBox[
      RowBox[{SubscriptBox["\[InvisiblePrefixScriptBase]", 
         TagBox[AdjustmentBox[
           TagBox["\[SelectionPlaceholder]", BraArgs, 
            style -> "BraKetArg"], BoxBaselineShift -> 0], UDScript, 
          style -> "LIndex"]], 
        SuperscriptBox[
         RowBox[{"\[LeftAngleBracket]", 
           AdjustmentBox[
            RowBox[{TagBox["\[SelectionPlaceholder]", BraKetArgs, 
               style -> "BraKetArg"], "\[VerticalSeparator]", 
              TagBox["\[SelectionPlaceholder]", BraKetArgs, 
               style -> "BraKetArg"]}], BoxBaselineShift -> 0], 
           "\[RightAngleBracket]"}], 
         TagBox[AdjustmentBox[
           TagBox["\[SelectionPlaceholder]", BraArgs, 
            style -> "BraKetArg"], BoxBaselineShift -> 0], UDScript, 
          style -> "RIndex"]]}], Braket, style -> "BraKetWrapper"],
    "ubraketd" ->
     TagBox[
      RowBox[{SuperscriptBox["\[InvisiblePrefixScriptBase]", 
         TagBox[AdjustmentBox[
           TagBox["\[SelectionPlaceholder]", BraArgs, 
            style -> "BraKetArg"], BoxBaselineShift -> 0], UDScript, 
          style -> "LIndex"]], 
        SubscriptBox[
         RowBox[{"\[LeftAngleBracket]", 
           AdjustmentBox[
            RowBox[{TagBox["\[SelectionPlaceholder]", BraKetArgs, 
               style -> "BraKetArg"], "\[VerticalSeparator]", 
              TagBox["\[SelectionPlaceholder]", BraKetArgs, 
               style -> "BraKetArg"]}], BoxBaselineShift -> 0], 
           "\[RightAngleBracket]"}], 
         TagBox[AdjustmentBox[
           TagBox["\[SelectionPlaceholder]", BraArgs, 
            style -> "BraKetArg"], BoxBaselineShift -> 0], UDScript, 
          style -> "RIndex"]]}], Braket, style -> "BraKetWrapper"],
    "ubraketu" ->
     TagBox[
      RowBox[{SuperscriptBox["\[InvisiblePrefixScriptBase]", 
         TagBox[AdjustmentBox[
           TagBox["\[SelectionPlaceholder]", BraArgs, 
            style -> "BraKetArg"], BoxBaselineShift -> 0], UDScript, 
          style -> "LIndex"]], 
        SuperscriptBox[
         RowBox[{"\[LeftAngleBracket]", 
           AdjustmentBox[
            RowBox[{TagBox["\[SelectionPlaceholder]", BraKetArgs, 
               style -> "BraKetArg"], "\[VerticalSeparator]", 
              TagBox["\[SelectionPlaceholder]", BraKetArgs, 
               style -> "BraKetArg"]}], BoxBaselineShift -> 0], 
           "\[RightAngleBracket]"}], 
         TagBox[AdjustmentBox[
           TagBox["\[SelectionPlaceholder]", BraArgs, 
            style -> "BraKetArg"], BoxBaselineShift -> 0], UDScript, 
          style -> "RIndex"]]}], Braket, style -> "BraKetWrapper"],
    "dbraketdu" ->
     TagBox[
      RowBox[{SubscriptBox["\[InvisiblePrefixScriptBase]", 
         TagBox[AdjustmentBox[
           TagBox["\[SelectionPlaceholder]", BraArgs, 
            style -> "BraKetArg"], BoxBaselineShift -> 0], UDScript, 
          style -> "LIndex"]], 
        SubsuperscriptBox[
         RowBox[{"\[LeftAngleBracket]", 
           AdjustmentBox[
            RowBox[{TagBox["\[SelectionPlaceholder]", BraKetArgs, 
               style -> "BraKetArg"], "\[VerticalSeparator]", 
              TagBox["\[SelectionPlaceholder]", BraKetArgs, 
               style -> "BraKetArg"]}], BoxBaselineShift -> 0], 
           "\[RightAngleBracket]"}], 
         TagBox[AdjustmentBox[
           TagBox["\[SelectionPlaceholder]", BraArgs, 
            style -> "BraKetArg"], BoxBaselineShift -> 0], UDScript, 
          style -> "RIndex"], 
         TagBox[AdjustmentBox[
           TagBox["\[SelectionPlaceholder]", BraArgs, 
            style -> "BraKetArg"], BoxBaselineShift -> 0], UDScript, 
          style -> "RIndex"]]}], Braket, style -> "BraKetWrapper"],
    "ubraketdu" ->
     TagBox[
      RowBox[{SuperscriptBox["\[InvisiblePrefixScriptBase]", 
         TagBox[
          AdjustmentBox[
           TagBox["\[SelectionPlaceholder]", BraArgs, 
            style -> "BraKetArg"], BoxBaselineShift -> 0], UDScript, 
          style -> "LIndex"]], 
        SubsuperscriptBox[
         RowBox[{"\[LeftAngleBracket]", 
           AdjustmentBox[
            RowBox[{TagBox["\[SelectionPlaceholder]", BraKetArgs, 
               style -> "BraKetArg"], "\[VerticalSeparator]", 
              TagBox["\[SelectionPlaceholder]", BraKetArgs, 
               style -> "BraKetArg"]}], BoxBaselineShift -> 0], 
           "\[RightAngleBracket]"}], 
         TagBox[AdjustmentBox[
           TagBox["\[SelectionPlaceholder]", BraArgs, 
            style -> "BraKetArg"], BoxBaselineShift -> 0], UDScript, 
          style -> "RIndex"], 
         TagBox[AdjustmentBox[
           TagBox["\[SelectionPlaceholder]", BraArgs, 
            style -> "BraKetArg"], BoxBaselineShift -> 0], UDScript, 
          style -> "RIndex"]]}], Braket, style -> "BraKetWrapper"],
    "dubraketd" ->
     TagBox[
      RowBox[{SubsuperscriptBox["\[InvisiblePrefixScriptBase]", 
         TagBox[AdjustmentBox[
           TagBox["\[SelectionPlaceholder]", BraArgs, 
            style -> "BraKetArg"], BoxBaselineShift -> 0], UDScript, 
          style -> "LIndex"], 
         TagBox[AdjustmentBox[
           TagBox["\[SelectionPlaceholder]", BraArgs, 
            style -> "BraKetArg"], BoxBaselineShift -> 0], 
          UDScript, style -> "LIndex"]], 
        SubscriptBox[
         RowBox[{"\[LeftAngleBracket]", 
           AdjustmentBox[
            RowBox[{TagBox["\[SelectionPlaceholder]", BraKetArgs, 
               style -> "BraKetArg"], "\[VerticalSeparator]", 
              TagBox["\[SelectionPlaceholder]", BraKetArgs, 
               style -> "BraKetArg"]}], BoxBaselineShift -> 0], 
           "\[RightAngleBracket]"}], 
         TagBox[AdjustmentBox[
           TagBox["\[SelectionPlaceholder]", BraArgs, 
            style -> "BraKetArg"], BoxBaselineShift -> 0], UDScript, 
          style -> "RIndex"]]}], Braket, style -> "BraKetWrapper"],
    "dubraketu" ->
     TagBox[
      RowBox[{SubsuperscriptBox["\[InvisiblePrefixScriptBase]", 
         TagBox[AdjustmentBox[
           TagBox["\[SelectionPlaceholder]", BraArgs, 
            style -> "BraKetArg"], BoxBaselineShift -> 0], UDScript, 
          style -> "LIndex"], 
         TagBox[
          AdjustmentBox[
           TagBox["\[SelectionPlaceholder]", BraArgs, 
            style -> "BraKetArg"], BoxBaselineShift -> 0], 
          UDScript, style -> "LIndex"]], 
        SuperscriptBox[
         RowBox[{"\[LeftAngleBracket]", 
           AdjustmentBox[
            RowBox[{TagBox["\[SelectionPlaceholder]", BraKetArgs, 
               style -> "BraKetArg"], "\[VerticalSeparator]", 
              TagBox["\[SelectionPlaceholder]", BraKetArgs, 
               style -> "BraKetArg"]}], BoxBaselineShift -> 0], 
           "\[RightAngleBracket]"}], 
         TagBox[AdjustmentBox[
           TagBox["\[SelectionPlaceholder]", BraArgs, 
            style -> "BraKetArg"], BoxBaselineShift -> 0], UDScript, 
          style -> "RIndex"]]}], Braket, style -> "BraKetWrapper"],
    "dubraketdu" ->
     TagBox[
      RowBox[{SubsuperscriptBox["\[InvisiblePrefixScriptBase]", 
         TagBox[AdjustmentBox[
           TagBox["\[SelectionPlaceholder]", BraArgs, 
            style -> "BraKetArg"], BoxBaselineShift -> 0.`], UDScript,
           style -> "LIndex"], 
         TagBox[AdjustmentBox[
           TagBox["\[SelectionPlaceholder]", BraArgs, 
            style -> "BraKetArg"], BoxBaselineShift -> 0.`], UDScript,
           style -> "LIndex"]], 
        SubsuperscriptBox[
         RowBox[{"\[LeftAngleBracket]", 
           AdjustmentBox[
            RowBox[{TagBox["\[SelectionPlaceholder]", BraKetArgs, 
               style -> "BraKetArg"], "\[VerticalSeparator]", 
              TagBox["\[SelectionPlaceholder]", BraKetArgs, 
               style -> "BraKetArg"]}], BoxBaselineShift -> 0.`], 
           "\[RightAngleBracket]"}], 
         TagBox[AdjustmentBox[
           TagBox["\[SelectionPlaceholder]", BraArgs, 
            style -> "BraKetArg"], BoxBaselineShift -> 0.`], UDScript,
           style -> "RIndex"], 
         TagBox[AdjustmentBox[
           TagBox["\[SelectionPlaceholder]", BraArgs, 
            style -> "BraKetArg"], BoxBaselineShift -> 0.`], UDScript,
           style -> "RIndex"]]}], Braket, style -> "BraKetWrapper"],
    "com" -> 
     SubscriptBox[
      RowBox[{"[", 
        RowBox[{"\[SelectionPlaceholder]", ",", 
          "\[SelectionPlaceholder]"}], "]"}], "-"],
    "op" -> OverscriptBox["\[SelectionPlaceholder]", "^"],
    "opd" -> SubscriptBox[OverscriptBox["\[SelectionPlaceholder]", "^"],"\[SelectionPlaceholder]"],
    "opu" -> TagBox[SuperscriptBox[OverscriptBox["\[SelectionPlaceholder]", "^"],"\[SelectionPlaceholder]"],superscript],
    "opdu" -> TagBox[SubsuperscriptBox[OverscriptBox["\[SelectionPlaceholder]", "^"],"\[SelectionPlaceholder]","\[SelectionPlaceholder]"],superscript],
    "hop" -> SuperscriptBox[OverscriptBox["\[SelectionPlaceholder]", "^"], "\[Dagger]"]

};


filteredOldInputAlias = DeleteCases[Options[nb1, InputAliases][[1, 2]], Alternatives @@ Thread[ qaInputAlias[[All,1]] -> _]];

SetOptions[nb1, InputAliases -> Join[filteredOldInputAlias, qaInputAlias]]


(* ::Subsubsection:: *)
(*End*)


Label[endNotations];


(* ::Subsection:: *)
(*End Package*)


Protect[ Evaluate[protected] ]

End[]

(*Protect[Evaluate[$Context <> "*"]]*)

SetAttributes[ Evaluate[ToExpression/@(Complement[Names @ "QuantumAlgebra`QuantumAlgebra`*", {"AutoLoadNotationPalette", "AutoLoadNavigatorPalette"}])], {ReadProtected}]

EndPackage[]
