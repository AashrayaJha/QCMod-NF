AttachSpec("QCMod.spec");
AttachSpec("coleman.spec");
import "applications.m": Qp_points;
import "singleintegrals.m": is_bad, xy_coordinates;
import "misc.m": function_field;
load "hodge_snip.m";


function inf_deg(Q_trans)
  d:=Degree(Q_trans);
  K := BaseRing(BaseRing(Q_trans));
  // find the points at infinity:
  Kx := RationalFunctionField(K);
  Kxy := PolynomialRing(Kx);
  FF := function_field(Q_trans); // function field of curve over K
  infplaces:=InfinitePlaces(FF);
  infplacesKinf := infplaces;
//     "a,b = ", a,b;
  "inf places", infplacesKinf;
  Kinf := K;
  for Kinf_new in infplacesKinf do
    if not IsOne(Degree(Kinf_new)) then
      // field generated by points at infinity
      //time norm_clos := NormalClosure(AbsoluteField(ResidueClassField(Kinf_new)));
      //time Kinf := Compositum(Kinf, norm_clos);
      time split := SplittingField(CharacteristicPolynomial(ResidueClassField(Kinf_new).1));
      time Kinf := Compositum(Kinf, split);
    end if;
  end for;
  dinf := AbsoluteDegree(Kinf);
  "degree of Kinf", dinf; 
  return dinf;
end function;


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

  return Q_trans, Q_trans_hom;
end function;

function good_transform(Q, A, embeddings : N := 6, empty := [1,2])
  // check if the transformation given by the 2x2 matrix A yields an
  // affine patch without bad residue disks with respect to the ith
  // embedding of F into Qp, where i runs through empty. 
  // Note that A acts on X and Z, but not Y.
  v1,v2 := Explode(embeddings);
  Q_trans := transform(Q, A);
  try 
    new_data1 := ColemanData(Q_trans, v1, N);
  catch e;
    e;
    return false, _, _, _, _, _;
  end try;
  new_Qppoints_1 := Qp_points(new_data1 : Nfactor := 1.5);
pts1 := [xy_coordinates(P, new_data1):P in new_Qppoints_1 | not is_bad(P, new_data1)];
  good1 := &and[not(is_bad(P, new_data1)) : P in new_Qppoints_1];
  if 1 notin empty or good1 then 
    try 
      new_data2 := ColemanData(Q_trans, v2, N : useU:=false);
    catch e;
      e;
      return false, _, _, _, _, _;
    end try;
    new_Qppoints_2 := Qp_points(new_data2 : Nfactor := 1.5);
    good2 := &and[not(is_bad(P, new_data2)) : P in new_Qppoints_2];
    if 2 notin empty or good2 then
      return true, Q_trans, new_data1, new_data2, new_Qppoints_1, new_Qppoints_2;
    end if;
  end if;
  return false, _, _, _, _, _;
end function;



function good_patches_modp(data1, data2, p : empty := [1,2], N := 6) 
  Q := data1`Q;
  v1 := data1`v;
  v2 := data2`v;
  F := BaseRing(BaseRing(Q));
  u := F.1;
  Qp1, sigma1 := Completion(F, v1);
  Qp2, sigma2 := Completion(F, v2);
  Fp := GF(p);
  Qppoints1 := Qp_points(data1);
  Qppoints2 := Qp_points(data2);
  xs1 := [xy_coordinates(pt, data1)[1] : pt in Qppoints1 | not pt`inf];
  xs2 := [xy_coordinates(pt, data2)[1] : pt in Qppoints2 | not pt`inf];
  sigmas := [sigma1, sigma2];
  xs := [xs1, xs2];
  good_coeffs_modp := [];
  for a,b in [0..p-1] do
    if &and[Fp!sigmas[i](a+b*F.1) notin xs[i] : i in empty] then
      "trying a+b*u", a+b*u;
      A := [[1,1], [1,-(a+b*u)]];
      bool, Q_trans, new_data1, new_data2, new_Qppoints_1, new_Qppoints_2 := good_transform(Q, A, [v1,v2] : empty := empty);
      "bool", bool;
      if bool then
        C2 := CurveFromBivariate(Q_trans);
        //DO1 := DixmierOhnoInvariants(C1);
        //DO2 := DixmierOhnoInvariants(C2);
        //assert DixmierOhnoInvariantsEqual(DO1, DO2);
        // The above works, but only checks isomorphism over \bar{K}
        C1 := CurveFromBivariate(Q);
        assert IsIsomorphicPlaneQuartics(C1,C2);
        infdeg := inf_deg(Q_trans);
        "works!", [a,b];
        Append(~good_coeffs_modp, [a,b,infdeg]);
      end if;
    end if;
  end for;
  return good_coeffs_modp;
end function;


function inf_poly(Q)
  Q_inf := transform(Q, [[0,1],[1,0]]);
  y := Parent(Q).1;
  infpoly := &+[Evaluate(Coefficient(Q_inf, i), 0)*y^i : i in [0..Degree(Q_inf)]];
  return infpoly;
end function;

K<zeta3> := CyclotomicField(3);
Q := Polynomial([PolynomialRing(CyclotomicField(3)) | [[ RationalField() | 1, 1 ], [ RationalField() | 1, -1 ], [ RationalField() | 0, 3 ], [ RationalField() | 0, -1 ]], [[ RationalField() | 0, -2 ], [ RationalField() | 0, 3 ], [ RationalField() | 0, -3 ], [ RationalField() | 2, 2 ]], [[ RationalField() | -3, 0 ]], [[ RationalField() | 2, 3 ], [ RationalField() | -1, 1 ]], [[ RationalField() | 1, 0 ]]]);
v1 := ideal<Integers(CyclotomicField(3)) | \[ 13, 0 ], \[ 4, 1 ]>;
v2 := ideal<Integers(CyclotomicField(3)) | \[ 13, 0 ], \[ 10, 1 ]>;
p := 13;
//data1 := ColemanData(Q, v1, 6);
//data2 := ColemanData(Q, v2, 6);
/*
good_modp := good_patches_modp(data1, data2, 13 : empty := [1], N := 8);
// One gets this A by noticing that under both embeddings there is no 
// Q13-point with x-coordinate 11 on the curve defined by Q=0, and that 
// the only bad disks are at infinity.
*/

Patches := [];
good_modp := [ [ 0, 12 ], [ 3, 11 ], [ 5, 10 ], [ 7, 4 ], [ 8, 1 ], [ 8, 9 ], [ 10, 3 ], [ 11, 0 ] ];

good_modp1 :=  
[   [ 0, 7], [ 0, 12], [ 1, 4], [ 1, 9], [ 2, 1], [ 2, 6], [ 3, 11], [ 4, 0], [ 4, 8], [ 5, 5], [ 5, 10], [ 6, 7], [ 7, 4], [ 7, 12], [ 8, 1], [ 8, 9], [ 9, 6], [ 9, 11], [ 10, 3], [ 10, 8], [ 11, 0], [ 11, 5], [ 12, 2], [ 12, 10]];


good_modp2 :=
[    [ 0, 2], [ 0, 3], [ 0, 8], [ 0, 12], [ 1, 3], [ 1, 6], [ 1, 12], [ 2, 3], [ 2, 7], [ 2, 10], [ 2, 11], [ 3, 2], [ 3, 7], [ 3, 11], [ 4, 2], [ 4, 5], [ 4, 6], [ 5, 2], [ 5, 6], [ 5, 9], [ 5, 10], [ 6, 0], [ 6, 1], [ 6, 6], [ 6, 10], [ 7, 1], [ 7, 4], [ 7, 5], [ 7, 10], [ 8, 1], [ 8, 5], [ 8, 8], [ 8, 9], [ 9, 0], [ 9, 5], [ 9, 12], [ 10, 0], [ 10, 3], [ 10, 4], [ 10, 9], [ 11, 0], [ 11, 4], [ 11, 7], [ 11, 8], [ 12, 4], [ 12, 8], [ 12, 11], [ 12, 12] ];


SetLogFile("gmp2.log");
"good_modp2";
for pair in good_modp do
  "trying a,b=", pair;
  //for i,j in [-200..200] do
  for i,j in [-0..0] do
    a := i*p+pair[1];
    b := j*p+pair[2];
    A := [[1,1], [1,-(a+b*u)]];
    Q_trans := transform(Q, A);
    //Q_trans := transform(Q_trans, [[0,1], [1,0]]);
    "Q_trans", Q_trans;
    //bool, Q_trans, new_data1, new_data2, new_Qppoints_1, new_Qppoints_2 := good_transform(Q, A, [v1,v2] : N := 6, empty := [1,2]);
//    if bool then
      d:=Degree(Q_trans);

      // find the points at infinity:

      Kx := RationalFunctionField(K);
      Kxy := PolynomialRing(Kx);

      FF := function_field(Q_trans); // function field of curve over K
      infplaces:=InfinitePlaces(FF);
      infplacesKinf := infplaces;
      //"a,b = ", a,b;
      //"inf places", infplacesKinf;

      Kinf := K;
      for Kinf_new in infplacesKinf do
        if not IsOne(Degree(Kinf_new)) then
          // field generated by points at infinity
          //time norm_clos := NormalClosure(AbsoluteField(ResidueClassField(Kinf_new)));
          //time Kinf := Compositum(Kinf, norm_clos);
          //time split := SplittingField(CharacteristicPolynomial(ResidueClassField(Kinf_new).1));
        s1:=CharacteristicPolynomial(ResidueClassField(Kinf_new).1);
s2:=Parent(s1)!inf_poly(Q_trans);
Evaluate(s1, -x) - s2;
s1;s2;
#GaloisGroup(s1);
#GaloisGroup(s2);
          galois := GaloisGroup(CharacteristicPolynomial(ResidueClassField(Kinf_new).1));
       //   time Kinf := Compositum(Kinf, split);
        end if;
      end for;
      //dinf := AbsoluteDegree(Kinf); //dinf := Degree(split);
      //"degree of Kinf", dinf; 
      if #galois lt 24 then
      //  "yay!!!"; a,b,galois,Q_trans;
      end if;
      //Append(~Patches, <pair, dinf>); //      <pair, Q_trans>; //      Append(~Patches, <pair, Q_trans>);
//    end if;
  end for;
end for;
//Qs := [t[2] : t in Patches];
//
//
/*
Qxy<x,y> := PolynomialRing(Rationals(), 2);
Qs :=
[
y^4 + ((21*u + 10)*x + (2*u + 3))*y^3 + (-3*x^2 + 6*x - 3)*y^2 + ((-3418*u - 2988)*x^3 + (711*u - 576)*x^2 + (123*u + 180)*x - 10*u - 2)*y + (-1703*u + 445)*x^4 + (1579*u - 1407)*x^3 + (264*u + 1083)*x^2 + (-146*u - 121)*x + 6*u,
y^4 + ((16*u + 12)*x + (2*u + 3))*y^3 + (-3*x^2 + 6*x - 3)*y^2 + ((-683*u - 1924)*x^3 + (738*u + 45)*x^2 + (69*u + 147)*x - 10*u - 2)*y + (-1212*u - 612)*x^4 + (1464*u + 186)*x^3 + (-153*u + 540)*x^2 + (-105*u - 114)*x + 6*u,
y^4 + ((12*u + 13)*x + (2*u + 3))*y^3 + (-3*x^2 + 6*x - 3)*y^2 + ((542*u - 720)*x^3 + (606*u + 330)*x^2 + (30*u + 120)*x - 10*u - 2)*y + (-584*u - 734)*x^4 + (970*u + 642)*x^3 + (-318*u + 198)*x^2 + (-74*u - 106)*x + 6*u,
y^4 + ((-2*u + 9)*x + (2*u + 3))*y^3 + (-3*x^2 + 6*x - 3)*y^2 + ((-170*u + 254)*x^3 + (-150*u + 114)*x^2 + (-54*u + 18)*x - 10*u - 2)*y + (162*u + 144)*x^4 + (-108*u + 48)*x^3 + (-72*u - 144)*x^2 + (12*u - 48)*x + 6*u,
y^4 + ((-9*u + 7)*x + (2*u + 3))*y^3 + (-3*x^2 + 6*x - 3)*y^2 + ((-853*u - 684)*x^3 + (-447*u - 255)*x^2 + (-96*u - 33)*x - 10*u - 2)*y + (-182*u + 133)*x^4 + (-56*u - 69)*x^3 + (177*u - 45)*x^2 + (55*u - 19)*x + 6*u,
y^4 + ((7*u + 15)*x + (2*u + 3))*y^3 + (-3*x^2 + 6*x - 3)*y^2 + ((1171*u + 812)*x^3 + (321*u + 561)*x^2 + (-24*u + 87)*x - 10*u - 2)*y + (258*u - 387)*x^4 + (192*u + 675)*x^3 + (-423*u - 189)*x^2 + (-33*u - 99)*x + 6*u,
y^4 + ((-7*u + 11)*x + (2*u + 3))*y^3 + (-3*x^2 + 6*x - 3)*y^2 + ((-1413*u - 398)*x^3 + (-591*u - 123)*x^2 + (-108*u - 15)*x - 10*u - 2)*y + (-10*u + 491)*x^4 + (-184*u - 231)*x^3 + (135*u - 219)*x^2 + (53*u - 41)*x + 6*u,
y^4 + ((-14*u + 9)*x + (2*u + 3))*y^3 + (-3*x^2 + 6*x - 3)*y^2 + ((-2330*u - 2662)*x^3 + (-966*u - 726)*x^2 + (-150*u - 66)*x - 10*u - 2)*y + (-978*u + 12)*x^4 + (336*u - 36)*x^3 + (540*u + 36)*x^2 + (96*u - 12)*x + 6*u
];


    
 
good_trans_1 := [u+2, 4*u+1, 6*u+9];
good_affines_1 := [];
for a in good_trans_1 do
  A := [[1,1], [1,-a]];
  bool, Q_trans, new_data1, new_data2, new_Qppoints_1, new_Qppoints_2 := good_transform(Q, A, [v1,v2] : N := 6, empty := [1]);
  bool;
  "a", a;
  "Q_trans", Q_trans;
  Append(~good_affines_1, Q_trans);
  infdeg(Q_trans);
end for;

Qy := -u*x^3*y + (2*u + 2)*x^3 + 3*u*x^2*y^2 - 3*u*x^2*y + (-u + 1)*x*y^3 + 3*u*x*y^2 + (u - 1)*x + (u + 1)*y^4 - 2*u*y^3 - 3*y^2 + (3*u + 2)*y + 1
*/

/*
zeropolys := [];
zeropolys1 := [];
zeropolys2 := [];

L<i,j> := FunctionField(K, 2);
Lx<x> := PolynomialRing(L);
Lxy<y> := PolynomialRing(Lx);
Q := y^4 + ((zeta3 - 1)*x + (3*zeta3 + 2))*y^3 - 3*y^2 + ((2*zeta3 + 2)*x^3 - 3*zeta3*x^2 + 3*zeta3*x - 2*zeta3)*y - zeta3*x^3 + 3*zeta3*x^2 + (-zeta3 + 1)*x + zeta3 + 1;
for pair in good_modp do
  "trying a,b=", pair;
  A := [[1,-pair[1]-i*p-(pair[2]+j*p)*zeta3], [0,1]];
  Q_trans := transform(Q, A);
  zeropoly := &+[Evaluate(Coefficient(Q_trans, i), 0)*y^i : i in [0..Degree(Q_trans)]];
  Append(~zeropolys, zeropoly);
end for;
for pair in good_modp1 do
  "trying a,b=", pair;
  A := [[1,-pair[1]-i*p-(pair[2]+j*p)*zeta3], [0,1]];
  Q_trans := transform(Q, A);
  zeropoly := &+[Evaluate(Coefficient(Q_trans, i), 0)*y^i : i in [0..Degree(Q_trans)]];
  Append(~zeropolys1, zeropoly);
end for;
for pair in good_modp2 do
  "trying a,b=", pair;
  A := [[1,-pair[1]-i*p-(pair[2]+j*p)*zeta3], [0,1]];
  Q_trans := transform(Q, A);
  zeropoly := &+[Evaluate(Coefficient(Q_trans, i), 0)*y^i : i in [0..Degree(Q_trans)]];
  Append(~zeropolys2, zeropoly);
end for;

*/