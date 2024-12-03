AttachSpec("~/GitHub/CHIMP/CHIMP.spec");
R<x> := PolynomialRing(Rationals());
// I picked the polredabs field to avoid confusions
F<z> := NumberFieldExtra(x^2-x+1:prec:=300);
z3 := -z;
z3^2+z3+1 eq 0;
R<a,b,c> := PolynomialRing(F,3);
f := a^4 + (z3-1)*a^3*b + (3*z3+2)*a^3*c -3*a^2*c^2 + (2*z3+2)*a*b^3 -\
3*z3*a*b^2*c + 3*z3*a*b*c^2 - 2*z3*a*c^3 - z3*b^3*c + 3*z3*b^2*c^2 +\
(-z3+1)*b*c^3 + (z3+1)*c^4;
C := Curve(ProjectiveSpace(R),f);

// behind the scenes picks the riemannsurface
_ := PeriodMatrix(C);
X := C`riesrf;

Cycles, IntersectionMat, SymplecticTransform := HomologyBasis(X);
Omega1:= HolomorphicDifferentials(X); //The output is x^iy^jdx/g where [i,j] are the three integer tuples first outputted
                                      // and g=d_y(f)
HeuristicEndomorphismRepresentation(C : Geometric:=true);

