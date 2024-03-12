AttachSpec("~/GitHub/CHIMP/CHIMP.spec");
Attach("pushpull.m");\

A2<x,y>:=AffineSpace(Rationals(),2);
Q := y^4 +(-13*x + 7)*y^3 + (22*x^2 - 4*x - 8)*y^2 + (8*x^3 - 28*x^2 + 12*x + 4)*y - 8*x^3 + 16*x^2 - 8*x;
Xp:=Curve(A2,Q);
X:=ProjectiveClosure(Xp);
print "The curve is", Xp;
Z:= FindEndoMatrix(X);
print "Endomatrix is", Z;
P0:=X![1,0,0];
print "A rational point on the curve is", P0;
//time Ci,Up:=ConstructCorrespondence(X, P0, Z); //This didn't finish in 2 hours.
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

