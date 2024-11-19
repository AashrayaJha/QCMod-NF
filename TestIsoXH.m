load "galrep.m";

function ThreeTorsionReps(E);
  // magma's ThreeTorsionReps from Geometry/CrvG1/3descent-testeq.m,
  // slightly adapted to work over Q(zeta3)
  a1,a2,a3,a4,a6 := Explode(aInvariants(E));
  b2,b3,b6,b8 := Explode(bInvariants(E));
  c4,c6 := Explode(cInvariants(E));
  P := PolynomialRing(BaseRing(E)); X := P.1;
  F := X^4 - 6*c4*X^2 - 8*c6*X - 3*c4^2;
  if c4 eq 0 then F := ExactQuotient(F,X); end if;
  FF := Factorization(F);
  FF := [Evaluate(f[1],X^2) : f in FF];
  FF := &cat[Factorization(f): f in FF]; 
  FF := [f[1] : f in FF];
  function torspt(L,phi)
    x := (phi^2 - b2)/12;
    y := (phi^4 - c4)/(48*phi) - (a1*x + a3)/2;
    return E(L)![x,y];
  end function;
  TT := <torspt(BaseRing(E),Roots(f)[1][1]) : f in FF | Degree(f) eq 1>;
  if c4 eq 0 then 
    flag,sqrt := IsSquare(-c6/6);
    if flag then 
      Append(~TT,E(BaseRing(E))![-b2/12,sqrt/12 + a1*b2/24 - a3/2]);
      Append(~TT,E(BaseRing(E))![-b2/12,-sqrt/12 + a1*b2/24 - a3/2]);
    else 
      L<u> := NumberField(X^2 + c6/6);
      Append(~TT,E(L)![-b2/12,u/12 + a1*b2/24 - a3/2]);
    end if;
  end if;
  F2 := [f : f in FF | Degree(f) gt 1];
  F2 := Sort(F2,func<f,g|Degree(f) - Degree(g)>);
  for f in F2 do
    L<u> := NumberField(f);
    Append(~TT,torspt(L,u));
  end for;
  return TT;
end function;


K<zeta> := CyclotomicField(3);
RR<A,B,C> := PolynomialRing(K,3);
j := (-3238903991430*zeta - 4786881515880)^3/9^27; // Aash's file
j:= 2^3 * 5^3 * 19^(-27) * (-15- 2*zeta)^3 * (452 +431*zeta)^3 * (53515427 + 41314914*zeta)^3; // RSZB

E := WeierstrassModel(EllipticCurveWithjInvariant(j));
reps := ThreeTorsionReps(E);
L<u> := Parent(reps[1,1]);
EL := ChangeRing(E, L);
f := ChangeRing(MinimalPolynomial(reps[1,1]), L);
rts := Roots(f);
//"starting Galois computation";
//G, local_rts, data := GaloisGroup(L);
eqn := DefiningEquation(EL);
y := Parent(eqn).2;
y1poly := UnivariatePolynomial(Evaluate(eqn, [rts[1,1], y,1]));
y2poly := UnivariatePolynomial(Evaluate(eqn, [rts[2,1], y,1]));
y1s := Roots(y1poly);
y2s := Roots(y2poly);
assert #y1s*#y2s gt 0;
y1 := y1s[1,1]; y2 := y2s[1,1];
P1 := EL![rts[1,1],y1]; P2 := EL![rts[2,1],y2];
B := [P1,P2];
I, GL2:= GaloisImage(E,B,3);


bases := [];
for i,j,k,l in [0..2] do  
  P1 := i*B[1]+j*B[2];
  P2 := k*B[1]+l*B[2];
  if P1 ne Zero(EL) and P2 notin [P1,-P1, Zero(EL)] and [P1,P2] notin bases then
    Append(~bases, [P1,P2]);
  end if;
end for;
Is := [GaloisImage(E,bas,3) : bas in bases];


H3 := sub<GL2 | [GL2!Matrix(GF(3), 2,2, list) :  list in [[0,26,4,6], [10,1,25,26]]]>;
&and[I eq H3 : I in Is];
assert IsNormal(H3, GL2);



 Z27 := Integers(27);
 GL227 := GL(2,Z27);
 h1 := Matrix(Z27, 2,2,[0,26,4,6]);
 h1 := GL227!h1;
 h2 := Matrix(Z27, 2,2,[10,1,25,26]);
 h2 := GL227!h2;                     
 H := sub<GL227 | [h1,h2]>;


 Z9 := Integers(9);
 GL29 := GL(2,Z9);
 h1 := Matrix(Z9, 2,2,[0,26,4,6]);
 h1 := GL29!h1;
 h2 := Matrix(Z9, 2,2,[10,1,25,26]);
 h2 := GL29!h2;                     
 H9 := sub<GL29 | [h1,h2]>;

pol9 := ExactQuotient(DivisionPolynomial(E, 9), DivisionPolynomial(E, 3));
L9<w> := NumberField(pol9);
M9<u> := SplittingField(pol9);
