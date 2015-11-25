(* ::Package:: *)

(* :Author: Yanbing Zhao(ustc_zzzz) *)

(* :Summary:
A Mathematica Package for Uncertainty Calculation in the Experiment of 
University Physics
*)


BeginPackage["Uncertainty`"]

If[!ValueQ[Confidence],Confidence::usage=
"Confidence gives confidence in uncentainty calculation. "];

If[!ValueQ[SetConfidence],SetConfidence::usage=
"SetConfidence[newConfidence] Set new confidence in uncertainty \
calculation. \nSetConfidence[] Reset confidence to default value(0.95). "];

If[!ValueQ[UValue],UValue::usage=
"UValue[value, uncertainty] gives uncertainty and value \
directly. \nUValue[value] gives value directly, but leave the uncertainty \
as default(0). \nUValue[value, {error, distribution}] gives value \
directly, while calculates the B type uncertainty with distribution given, \
\"Normal\", \"Uniform\", and \"Triangular\" are avaliable. For example, \
\"Normal\" means normal distribution, and \"Uniform\" means uniform \
distribution. \nUValue[value, {error}] gives value directly and calculates \
the uncertainty with the default distribution(normal distribution). \
\nUValue[{value1, value2, ...}, {error, distribution}] gives value as the \
average of a list given, calculates A type uncertainty by the list, and \
calculates B type uncertainty with the distribution given. "];

If[!ValueQ[U],U::usage=
"U[uvalue] gives the uncertainty of uncertainty value. \
For example, U[Uvalue[v, u]] returns u. "];

If[!ValueQ[V],V::usage=
"V[uvalue] gives the value of uncertainty value. \
For example, V[Uvalue[v, u]] returns v. "];

If[!ValueQ[TFactor],TFactor::usage=
"TFactor[n] gives factor of student t distribution \
with the confidence given. "];

If[!ValueQ[UValueA],UValueA::usage=
"UValueA[{value1, value2, ...}] gives A type uncertainty of a list "];

If[!ValueQ[UValueB],UValueB::usage=
"UValueB[error, distribution] gives B type uncertainty of an \
error with the distribution given, search UValue for more info. \
\nUValueB[error] gives B type uncertainty of the error \
but the default distribution(normal distribution)"];

Unprotect[Confidence,SetConfidence,UValue,U,V,TFactor,UValueA,UValueB];


Begin["`Private`"];

Confidence::illegal="The argument `1` should be between 0 and 1. ";
UValue::nonnegative="The argument `1` should be a non-negative number. "
$confidence=0.95;
Confidence:=$confidence;
SetConfidence[]:=($confidence=0.95);
SetConfidence[x_]/;If[x<1&&x>0,True,
    Message[Confidence::illegal,x],False]:=($confidence=x);
UValue[v_]:=UValue[v,0];
UValue[v_,u_?(Composition[Not,NonNegative])]:=(
    Message[UValue::nonnegative,u];UValue[v,Norm[u]]);
UValue[v_,{e_,d___}]:=UValue[v,UValueB[e,d]];
UValue[v_,u_,error__]:=UValue[v,
    Norm[Replace[{u,error},{e_,d___}->UValueB[e,d],{1}]]];
UValue[varlist_List,error___]:=UValue[Mean[varlist],
    UValueA[varlist],error];
U[UValue[v_,u_]]:=u;
V[UValue[v_,u_]]:=v;
TFactor[n_/;n>1]:=InverseCDF[StudentTDistribution[n-1],(1.+$confidence)/2];
UValueA[varlist_List/;Length[varlist]>1]:=Module[
    {\[Sigma]=StandardDeviation[varlist],n=Length[varlist]},
    TFactor[n]/Sqrt[n]\[Sigma]];
UValueB[error_,"Normal"]:=Module[{},
    error(-Sqrt[2]InverseErfc[1.+$confidence]/3)];
UValueB[error_,"Uniform"]:=Module[{},
    error (1.+$confidence)/2];
UValueB[error_,"Triangular"]:=Module[{},
    error(1-Sqrt[1.-$confidence])];
UValueB[error_:0]:=UValueB[error,"Normal"];
function_?(
    MemberQ[Attributes[#],NumericFunction]&)[f___,uv_UValue,b___]^:=Module[
	{length=Length[{f,uv,b}],
	valuelist=Replace[{f,uv,b},x_UValue->V[x],{1}],
	varlist=Flatten@Position[{f,uv,b},_UValue]},
	UValue[function@@valuelist,Norm@Table[U[Part[{f,uv,b},i]]
    (Derivative@@Array[If[i==#,1,0]&,length])[function]@@valuelist,
    {i,varlist}]]];

End[];


Protect[Confidence,SetConfidence,UValue,U,V,TFactor,UValueA,UValueB];

EndPackage[]
