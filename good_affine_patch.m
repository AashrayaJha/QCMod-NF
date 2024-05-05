AttachSpec("QCMod.spec");
AttachSpec("coleman.spec");
import "applications.m": Qp_points;
import "singleintegrals.m": is_bad, xy_coordinates;

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

function good_transform(Q, A, embeddings : N := 7)
  // check if the transformation given by the 2x2 matrix A yields an
  // affine patch without bad residue disks. A acts on X and Z, but not Y.
  v_1,v_2 := Explode(embeddings);
  Q_trans := transform(Q, A);
  new_data1 := ColemanData(Q_trans, v_1, N : useU:=false);
  new_Qppoints_1 := Qp_points(new_data1 : Nfactor := 1.5);
//pts1 := [xy_coordinates(P, new_data1):P in new_Qppoints_1 | not is_bad(P, new_data1)];
  good1 := &and[not(is_bad(P, new_data1)) : P in new_Qppoints_1];
  if good1 then 
    new_data2 := ColemanData(Q_trans, v_2, N : useU:=false);
    new_Qppoints_2 := Qp_points(new_data2 : Nfactor := 1.5);
    good2 := &and[not(is_bad(P, new_data2)) : P in new_Qppoints_2];
    if good2 then
      return true, Q_trans, new_data1, new_data2, new_Qppoints_1, new_Qppoints_2;
    end if;
  end if;
  return false, _, _, _, _;
end function;

K<u> := CyclotomicField(3);
Q := Polynomial([PolynomialRing(CyclotomicField(3)) | [[ RationalField() | 1, 1 ], [ RationalField() | 1, -1 ], [ RationalField() | 0, 3 ], [ RationalField() | 0, -1 ]], [[ RationalField() | 0, -2 ], [ RationalField() | 0, 3 ], [ RationalField() | 0, -3 ], [ RationalField() | 2, 2 ]], [[ RationalField() | -3, 0 ]], [[ RationalField() | 2, 3 ], [ RationalField() | -1, 1 ]], [[ RationalField() | 1, 0 ]]]);
v_1 := ideal<Integers(CyclotomicField(3)) | \[ 13, 0 ], \[ 4, 1 ]>;
v_2 := ideal<Integers(CyclotomicField(3)) | \[ 13, 0 ], \[ 10, 1 ]>;
p := 13;
A := [[K!1,1], [1,-11]];
// One gets this A by noticing that under both embeddings there is no 
// Q13-point with x-coordinate 11 on the curve defined by Q=0, and that 
// the only bad disks are at infinity.

bool, Q_trans, new_data1, new_data2, new_Qppoints_1, new_Qppoints_2 := good_transform(Q, A, [v_1,v_2]);
"transformed Q:";
Q_trans;
// Check that there are no bad disks
assert bool;
// Check that the curves are actually isomorphic
C1 := CurveFromBivariate(Q);
C2 := CurveFromBivariate(Q_trans);
//DO1 := DixmierOhnoInvariants(C1);
//DO2 := DixmierOhnoInvariants(C2);
//assert DixmierOhnoInvariantsEqual(DO1, DO2);
// The above works, but only checks isomorphism over \bar{K}
assert IsIsomorphicPlaneQuartics(C1,C2);
