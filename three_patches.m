
AttachSpec("QCMod.spec");
AttachSpec("coleman.spec");
import "applications.m": Qp_points;
import "singleintegrals.m": is_bad, xy_coordinates;
import "misc.m": function_field;
load "hodge_snip.m";


function transform(Q, A)
  
  // Apply the projective transformation given by the 2x2 matrix A to the 
  // curve defined by Q(x,y)=0. A acts on X and Z, but not Y.
  y := Parent(Q).1;
  x := BaseRing(Parent(Q)).1;
  K := BaseRing(BaseRing(Q));
  PK3<X,Y,Z>:=PolynomialRing(K,3);
  Q_dehom:=PK3!0;
  d := Degree(Q);
  for i:=0 to d do
    for j:=0 to Degree(Coefficient(Q,i)) do
      Q_dehom +:= Coefficient(Coefficient(Q,i),j)*Y^i*X^j;
    end for;
  end for;
  //We should probably have A^(-1), since det(A) is not always 1
  transformed_vars := [A[2,2]*X-A[1,2]*Z,Y,-A[2,1]*X+A[1,1]*Z];
  Q_hom := Homogenization(Q_dehom, Z);
  Q_trans_hom:=Evaluate(Q_hom, transformed_vars);
  Q_trans_dehom := Evaluate(Q_trans_hom, [X,Y,1]);
  Q_trans := Parent(Q)!0;
  for i:=0 to Degree(Q_trans_dehom, Y) do
    for j:=0 to Degree(Q_trans_dehom, Y)  do
      Q_trans +:= K!Coefficient(Coefficient(Q_trans_dehom, Y,i), X, j)*x^j*y^i;
    end for;
  end for;

  return Q_trans;
end function;

K<zeta3> := CyclotomicField(3);
Q := Polynomial([PolynomialRing(CyclotomicField(3)) | [[ RationalField() | 1, 1 ], [ RationalField() | 1, -1 ], [ RationalField() | 0, 3 ], [ RationalField() | 0, -1 ]], [[ RationalField() | 0, -2 ], [ RationalField() | 0, 3 ], [ RationalField() | 0, -3 ], [ RationalField() | 2, 2 ]], [[ RationalField() | -3, 0 ]], [[ RationalField() | 2, 3 ], [ RationalField() | -1, 1 ]], [[ RationalField() | 1, 0 ]]]);
v1 := ideal<Integers(CyclotomicField(3)) | \[ 13, 0 ], \[ 4, 1 ]>;
v2 := ideal<Integers(CyclotomicField(3)) | \[ 13, 0 ], \[ 10, 1 ]>;
p := 13;
N := 6;

C := CurveFromBivariate(Q);


y := Parent(Q).1;
x := BaseRing(Parent(Q)).1;
// Affine patch z=1:
Qz := Q;
// Affine patch x=1:
Qx := transform(Q, [[0,1],[1,0]]);


// Affine patch y=1
Qy := -zeta3*x^3*y + (2*zeta3 + 2)*x^3 + 3*zeta3*x^2*y^2 - 3*zeta3*x^2*y + (-zeta3 + 1)*x*y^3 + 3*zeta3*x*y^2 + (zeta3 - 1)*x + (zeta3 + 1)*y^4 - 2*zeta3*y^3 - 3*y^2 + (3*zeta3 + 2)*y + 1;
// Now make monic
Qy := transform(Qy, [[2,2],[0,1]]);

lc := LeadingCoefficient(Qy);
Qy := Parent(Qy)!(lc^3*Evaluate(Qy, y/lc));
//IsIsomorphicPlaneQuartics(C, CurveFromBivariate(Qy));
time datay1 := ColemanData(Qy, v1, N);
time datax1 := ColemanData(Qx, v1, N);
time dataz1 := ColemanData(Qz, v1, N);
//time datax2 := ColemanData(Qx, v2, N);
