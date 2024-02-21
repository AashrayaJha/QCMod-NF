AttachSpec("~/GitHub/CHIMP/CHIMP.spec");
Attach("pushpull.m");\

R_Q<x> := PolynomialRing(Rationals());
f :=  x^6 + 8*x^4 + 10*x^3 + 20*x^2 + 12*x + 9; p := 3; print "p = ",p;
X := HyperellipticCurve(f);
Xp, rho := ChangeCoordinatesHyp(X);
Z:= FindEndoMatrix(Xp);
print "Endomatrix is", Z;
P0:=AnyRationalPoint(Xp);
print "A rational point on the curve is", P0;
time Ci,Up:=ConstructCorrespondenceByCantor(Xp, P0, Z); //In this case C is a sequence of exactly one curve.

KUp := FunctionField(Up);
u := KUp.1;
KXp := FunctionField(Xp);
w := KXp.1;
basisDR, coeffsDR, changeOfBasis := ConstructDifferentials(X, Xp, rho, KUp);
zo, deg := EndoAction(Xp, Ci, Up, basisDR); //zo is action on allof H^1_{dR}
//print "the action on H^1_{dR} is", Type(zo) ;

KX<a,b> := FunctionField(X);
zoX := [Pullback(rho, KXp!(zo[i]/Differential(u))*Differential(w)) : i in [1 .. #zo]]; //Back on X from Xp



//load "data/quartic-test-data.m";
// Q := Polynomial([PolynomialRing(RationalField()) | [0, 0, 1, 1], [-1, 0, -1, -1], [0, 1, -1], [1]]);
// data1:=test_data_list[1];
// A2:=AffineSpace(Rationals(),2);
// y:=A2.1;x:=A2.2;
// // printf "%o,Q";
// Q:=y^3 + (-x^2+x)*y^2+(-x^3-x^2-1)*y+x^3+x^2;
// X:=Curve(A2,Q);
// FindEndoMatrix(X);

