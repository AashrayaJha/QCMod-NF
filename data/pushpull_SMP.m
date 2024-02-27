AttachSpec("~/GitHub/CHIMP/CHIMP.spec");
Attach("pushpull.m");

R<x> := PolynomialRing(Rationals());
F<z> := NumberField(x^2-x+1);
z3 := -z;
z3^2+z3+1 eq 0;
R<a,b,c> := PolynomialRing(F,3);
f := a^4 + (z3-1)*a^3*b + (3*z3+2)*a^3*c -3*a^2*c^2 + (2*z3+2)*a*b^3 -
3*z3*a*b^2*c + 3*z3*a*b*c^2 - 2*z3*a*c^3 - z3*b^3*c + 3*z3*b^2*c^2 +
(-z3+1)*b*c^3 + (z3+1)*c^4;
X := Curve(ProjectiveSpace(R),f);
P0:=X![z,0,1];
print "A rational point on the curve is", P0;
Z:=Matrix(F,3,3,[    0    ,  -z+1 ,     -z+1, z  ,    0, z,z ,   -z+1  ,    0]);
print "Endomatrix is", Z;
//Z1 is the same matrix one obtains from ConstructCorrespondencebyCantor. NumberFieldsExtra does not play nicely
//with it so we reproduce it here.

//time _, cant := CantorFromMatrixAmbientSplit(X, P0, X, P0, Z : LowerBound := 1); //takes 58 seconds
//time _, D := DivisorFromMatrixAmbientSplit(X, P0, X, P0, Z: LowerBound := 1); //Takes 150 seconds

output_file:="data/NF-example-test-data.m";
// out_1:= Sprintf("correspondences_Cantor:=%m;",cant); //Need to rewrite constants, define parent ring
// out_2:= Sprintf("correspondences_Divisor :=%m;",D);  //Need to rewrite constants by hand, and change ring to R as below.
// Write(output_file, out_1);
// Write(output_file, out_2);

// eqs := DefiningEquations(D);
// R<y2,y1,x2,x1> := Parent(eqs[1]);
// I := DefiningIdeal(D);
// J := Saturation(I,Evaluate(Denominator(cant[1]),[x1,y1])); //This took 3300 seconds!
// out:=Sprintf("Saturation :=%m;", J);
// U := X`U;
// Write(output_file, out);
//comps := IrreducibleComponents(Scheme(Ambient(D),J));
print "Calculated correspondencebyCantor";
//Run this after loading from stored data
// Zi := [c : c in comps | Dimension(c) gt 0 ];
// Ci := [Curve(Z) : Z in Zi];
// KUp := FunctionField(Up);
// u := KUp.1;
// KXp := FunctionField(Xp);
// w := KXp.1;
// basisDR, coeffsDR, changeOfBasis := ConstructDifferentials(X, Xp, rho, KUp);
// zo, deg := EndoAction(Xp, Ci, Up, basisDR); //zo is action on allof H^1_{dR}
// //print "the action on H^1_{dR} is", Type(zo) ;

// KX<a,b> := FunctionField(X);
// zoX := [Pullback(rho, KXp!(zo[i]/Differential(u))*Differential(w)) : i in [1 .. #zo]]; //Back on X from Xp
// R<x> := PolynomialRing(Rationals());
// F<z> := NumberFieldExtra(x^2-x+1:prec:=300);
// F1<z1>:=NumberFieldExtra(x^2-x+1);
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
// time Z:= HeuristicEndomorphismRepresentation(X : Geometric:=true)[2][1];
// print "Endomatrix is", Z;
// P0:=X![z,0,1];
// print "A rational point on the curve is", P0;
// //time Ci,Up:=ConstructCorrespondence(X, P0, Z); 
//time Ci,Up:=ConstructCorrespondenceByCantor(X, P0, Z); //In this case C is a sequence of exactly one curve.
