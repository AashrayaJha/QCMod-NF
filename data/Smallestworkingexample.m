AttachSpec("~/GitHub/CHIMP/CHIMP.spec");
R<x> := PolynomialRing(Rationals());
F1<z1> := NumberField(x^2-x+1);
z3_1 := -z1;
z3_1^2+z3_1+1 eq 0;
R1<a1,b1,c1> := PolynomialRing(F1,3);
f1 := a1^4 + (z3_1-1)*a1^3*b1 + (3*z3_1+2)*a1^3*c1 -3*a1^2*c1^2 + (2*z3_1+2)*a1*b1^3 -
3*z3_1*a1*b1^2*c1 + 3*z3_1*a1*b1*c1^2 - 2*z3_1*a1*c1^3 - z3_1*b1^3*c1 + 3*z3_1*b1^2*c1^2 +
(-z3_1+1)*b1*c1^3 + (z3_1+1)*c1^4;
Xp := Curve(ProjectiveSpace(R1),f1);
P0:=Xp![z1,0,1];
Z1:=Transpose(Matrix(F1,3,3,[0, z1 ,z1, -z1 + 1, 0, -z1 + 1 ,-z1 + 1, z1 , 0])); 
//This is obtained from HeuristicEndomorphismRepresentation
//time Ci,Up:=ConstructCorrespondenceByCantor(Xp, P0, Z1); //In this case C is a sequence of exactly one curve.
//time _,cant:=CantorFromMatrixAmbientSplit(Xp, P0, Xp, P0, Z1); 
//time _,D:=DivisorFromMatrixAmbientSplit(Xp, P0, Xp, P0, Z1: LowerBound := 1) ;

U := Xp`U;
eqs := DefiningEquations(D);
R<y2,y1,x2,x1> := Parent(eqs[1]);
I := DefiningIdeal(D);
J := Saturation(I,Evaluate(Denominator(cant[1]),[x1,y1]));
comps := IrreducibleComponents(Scheme(Ambient(D),J));
Zi := [c : c in comps | Dimension(c) gt 0 ];
Ci := [Curve(Z) : Z in Zi];
//Below is code I used to find Z1.

// F<z> := NumberFieldExtra(x^2-x+1:prec:=300);
// z3 := -z;
// z3^2+z3+1 eq 0;
// z1^2-z1+1:=0;
// R<a,b,c> := PolynomialRing(F,3);

// f := a^4 + (z3-1)*a^3*b + (3*z3+2)*a^3*c -3*a^2*c^2 + (2*z3+2)*a*b^3 -
// 3*z3*a*b^2*c + 3*z3*a*b*c^2 - 2*z3*a*c^3 - z3*b^3*c + 3*z3*b^2*c^2 +
// (-z3+1)*b*c^3 + (z3+1)*c^4;
// X := Curve(ProjectiveSpace(R),f);
// print "The curve is", X;
// CoefficientField(X);

// //Do this in lieu  of FindEndos in pushpull.m
// time Z:= Transpose(HeuristicEndomorphismRepresentation(X : Geometric:=true)[2][1]);
// print "Endomatrix is", Z;
// P0:=X![z,0,1];
// print "A rational point on the curve is", P0;
