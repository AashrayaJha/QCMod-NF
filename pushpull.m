declare verbose LocalHeights, 3;

intrinsic ChangeCoordinatesHyp(X::CrvHyp) -> RngUPolElt, RngUPolElt, MapSch
  {input: a monic polynomial of even degree f defining a hyperelliptic curve X:y^2=f(x)
  outut: two polynomials h, g such that the curve Xp:y^2+hy=g is isomorphic to the above.
  Also returns the isomorphism X -> Xp}
  f,hh := HyperellipticPolynomials(X);
  assert IsEven(Degree(f)) and IsZero(hh);
  require IsMonic(f) : "Non-monic polynomials not yet supported";
  //could change to square leading coeffs easily, but for now we will leave it at monic
  R<x> := Parent(f);
  _<yy> := PolynomialRing(Parent(f));
  mon := x^(Degree(f) div 2);
  adding := mon;
  yPrime := yy+mon;
  g := f - R!Evaluate(yPrime^2,0);
  while Degree(g) ge Degree(mon) do
    yPrime +:= LeadingCoefficient(g)*x^(Degree(g)-Degree(mon))/2;
    adding +:= LeadingCoefficient(g)*x^(Degree(g)-Degree(mon))/2;
    g := f - R!Evaluate(yPrime^2,0);
  end while;
  h := R!Coefficient(yPrime^2,1);
  Xp := HyperellipticCurve(g,h);

  changeOfY := X.2 - &+[Coefficient(adding,i)*X.1^i : i in [0..Degree(adding)]];
  coefs, mons := CoefficientsAndMonomials(changeOfY);
  deg := Degree(changeOfY);
  changeOfY := &+[coefs[i]*mons[i]*X.3^(deg-Degree(mons[i])) : i in [1..#mons]];
  rho := map< X -> Xp | [X.1, changeOfY, X.3]>;
  _, rho := IsIsomorphism(rho);
  return Xp,rho;
end intrinsic;

intrinsic FindEndoMatrix(X::Crv : tracezero := true) -> AlgMatElt
    {Tries to return a nontrival trace 0 endomorphism, given as Z_* on H^0.}
    g := Genus(X);
    Per := PeriodMatrix(CurveExtra(X));
    G,m := GeometricEndomorphismRepresentation(Per, RationalsExtra());
    E := EndomorphismRepresentation(G,RationalsExtra(),m);

    g := Genus(X);
    require #E gt 1 : "Neron-Severi rank is not large enough.";
    n := #E;
    for i in [1..n] do
      if E[i][1] eq IdentityMatrix(Rationals(),g) then
        continue;
      else
        Z := E[i][1];
        if IsNilpotent(Z) then //why this condition??
          continue;
        end if;
        break;
      end if;
    end for;

    if tracezero then
      tr := Integers()!Trace(Z);
      if IsZero(tr) then
        vprintf LocalHeights, 1: "Choosing %o with minimal polynomial %o.\n", Z, MinimalPolynomial(Z);
        return Transpose(Z);
      else
        m := LCM(tr,Ncols(Z));
        newZ := Transpose((m div tr)*Z- (m div Ncols(Z))*IdentityMatrix(Rationals(),Ncols(Z)));
        assert IsZero(Trace(newZ));
        vprintf LocalHeights, 1: "Choosing %o with minimal polynomial %o.\n", newZ, MinimalPolynomial(newZ);
        return newZ;
      end if;
    end if;
    return Transpose(Z);
end intrinsic;

intrinsic AnyRationalPoint(X::CrvHyp) -> Pt
    {Any rational point}
    C, pi := SimplifiedModel(X);
    f := HyperellipticPolynomials(C);
    d := LCM([Denominator(c): c in Coefficients(f)]);
    Xprime := HyperellipticCurve(f* d^2);
    b, phi := IsIsomorphic(Xprime, X);
    pts := RationalPoints(Xprime : Bound := 1000);
    finitepts := [phi(pt) : pt in pts | phi(pt)[3] ne 0];
    finiteptsNotWeierstrass := [p: p in finitepts | not IsWeierstrassPlace(Place(X!p))];
    require not #finiteptsNotWeierstrass eq 0 : "There are no rational finite non-Weierstrass points.";
    return finiteptsNotWeierstrass[1];
end intrinsic;

intrinsic ConstructCorrespondence(Xp::Crv, P0::Pt, T::AlgMatElt) -> SeqEnum, Crv
    {Construct the correspondence associated to T}
    T := Transpose(T);
    vprintf LocalHeights, 1: "Constructing Correspondence.\n";
    _, D := DivisorFromMatrixAmbientSplit(Xp, P0, Xp, P0, T : LowerBound := 1);
    U := Xp`U;
    eqs := DefiningEquations(D);
    R<y2,y1,x2,x1> := Parent(eqs[1]);
    vprintf LocalHeights, 1: "Computing irreducible components.\n";
    comps := IrreducibleComponents(D);
    Zi := [c : c in comps | Dimension(c) gt 0 ];
    Ci := [Curve(Z) : Z in Zi];
    return Ci, U;
end intrinsic;

intrinsic ConstructCorrespondenceByCantor(Xp::Crv, P0::Pt, T::AlgMatElt) -> SeqEnum, Crv
  {Construct the correspondence associated to T}
 
  vprintf LocalHeights, 1: "Constructing Correspondence using Cantor representation.\n";
  time _, cant := CantorFromMatrixAmbientSplit(Xp, P0, Xp, P0, T : LowerBound := 1);
  vprintf LocalHeights, 1: "Constructing Correspondence from divisor.\n";
  time _, D := DivisorFromMatrixAmbientSplit(Xp, P0, Xp, P0, T: LowerBound := 1);
  U := X`U;
  eqs := DefiningEquations(D);
  R<y2,y1,x2,x1> := Parent(eqs[1]);
  I := DefiningIdeal(D);
  vprintf LocalHeights, 1: "Saturating ideal.\n";
  J := Saturation(I,Evaluate(Denominator(cant[1]),[x1,y1]));
  vprintf LocalHeights, 1: "Computing irreducible components.\n";
  comps := IrreducibleComponents(Scheme(Ambient(D),J));
  Zi := [c : c in comps | Dimension(c) gt 0 ];
  Ci := [Curve(Z) : Z in Zi];
  return Ci, U;
end intrinsic;

// intrinsic TraceByGroebnerBasis(pi::MapSch, g::FldFunFracSchElt, C::Crv) -> .
//   {Take the trace of g by pi}

//   require IsProjective((Domain(pi))) and IsProjective((Codomain(pi))): "The map must have projective domain and codomain";
//   degpi := Degree(pi);
//   gnum := Numerator(g);
//   gden := Denominator(g);

//   deqns := DefiningEquations(C);

//   Rz<z, x1, y1, x2, y2> := PolynomialRing(Rationals(), 5);
//   gnum := Evaluate(gnum, [y2, y1, x2, x1]);
//   gden := Evaluate(gden, [y2, y1, x2, x1]);
//   deqns := [ Evaluate(d, [y2, y1, x2, x1]) : d in deqns];

//   I := ideal<Rz | deqns cat [gden*z -gnum]>;  // this is the correspondence and the two curves
//   elimI := EliminationIdeal(I, {x1, y1, z});

//   R1<xx1,yy1> := PolynomialRing(Rationals(),2);
//   R1z<zz> := PolynomialRing(R1);

//   polys := [ Evaluate(zpoly, [zz, xx1, yy1, 0, 0]) : zpoly in Basis(elimI)];
//   sortedbyz := Sort(polys,func<x,y | Degree(x)-Degree(y)>);
//   sortedbyz := [ s : s in sortedbyz |Degree(s) gt 0];
//   zpoly := sortedbyz[1]; //this is the smallest z-degree poly (assert that z-degree divides the degpi)

//   assert degpi mod Degree(zpoly) eq 0;

//   z2coef := Coefficient(zpoly,Degree(zpoly) - 1);
//   lczpoly := Coefficient(zpoly, Degree(zpoly));
//   scalefactor := degpi / Degree(zpoly);
//   tr := scalefactor* (-z2coef/lczpoly);
//   KU<u,v> := FunctionField(Codomain(pi));
//   tr := Evaluate(tr,[u,v]);
//   vprint LocalHeights, 2 : "The trace is", tr;

//   return tr;
// end intrinsic;

intrinsic Trace(pi::MapSch, g::FldFunFracSchElt, C::Crv) -> .
 {}
  require IsProjective((Domain(pi))) and IsProjective((Codomain(pi))): "The map must have projective domain and codomain";
  degpi := Degree(pi);
  norms := [Pushforward(pi,g-a) : a in [0..degpi]];
  //Norm(s-a) = constant term of P(x+a)
  Rx<x>:=PolynomialRing(Rationals());
  M := [[Coefficient((x+i)^j, 0): j in [0 .. degpi]] : i in [0 ..degpi]];
  Minv := Matrix(M)^(-1);
  ans := [&+[row[j]*norms[j] : j in [1..Ncols(Minv)]] : row in Rows(Minv)];
  lc := ans[degpi+1];
  tr :=  -ans[degpi]/lc;
  vprint LocalHeights, 2 : "The trace is", tr;
  return tr;

end intrinsic;


intrinsic EndoAction(Xp::Crv, Ci::SeqEnum, U::Crv, omegas::SeqEnum) -> SeqEnum, SeqEnum
    {Act by the corresopndence Ci on set of omegas}
    fp, hp := HyperellipticPolynomials(Xp);
    R<y2,y1,x2,x1> := Parent(Equations(Ci[1])[1]);
    KU := FunctionField(U);
    u := KU.1;
    v := KU.2;
    omega0 := Differential(u)/(2*v + Evaluate(hp,u));
    PU := ProjectiveClosure(U);

    Zomegas := [];
    for j->C in Ci do
      Append(~Zomegas, []);
      PC := ProjectiveClosure(C);
      vprint LocalHeights, 1 : "Computing the action of the endomorphism on differential forms";
      pi1proj := map<PC-> PU | [x1,y1,1]>;
      pi2proj := map<PC-> PU | [x2,y2,1]>;
      d1 := Degree(pi1proj);
      d2 := Degree(pi2proj);
      //Z^*(omega) = pi1_*(pi2^*(omega)) = pi1_*(endo) . omega_0
      for omega in omegas do
        endo := Pullback(pi2proj, omega)/Pullback(pi1proj, omega0);
        tr := Trace(pi1proj, endo, C);
        //tr2 := TraceByGroebnerBasis(pi1proj,endo,C);
        Zomega := tr*omega0;
        Append(~Zomegas[j], Zomega);
      end for;
    end for;
    ans := [&+[Zomegas[j][i] : j in [1..#Ci]] : i in [1..#omegas]];
    return ans, [d1,d2];
end intrinsic;

intrinsic ConstructDifferentials(X ::CrvHyp, Xp ::CrvHyp, rho::MapIsoSch, KUp::FldFunFracSch) -> SeqEnum, SeqEnum
    {Given isomorphic hypereliptic curves X, Xp (under rho:X to Xp), with X of the form y^2 + 0y= fX
     it returns a basis for the second kind differentials on the curve Xp of the form v^2 + hp v = fp by first finding them on X}
    g := Genus(X);
    R := BaseRing(X);
    fp, hp := HyperellipticPolynomials(Xp);
    u := KUp.1;
    v := KUp.2;

    P<t> := PuiseuxSeriesRing(R); //R = Rationals(), x = 1/t parameter
    fX,hX := HyperellipticPolynomials(X);
    assert IsZero(hX);
    sqrtf := Sqrt(Evaluate(fX,1/t));  //y^2 = f(1/t), sqrt(f(1/t)) = y
    //y = v +1/2h

    differentials := [-(1/t^2) * (1/t)^(i-1)/(2*sqrtf) : i in [1..2*g+1]];
    //differentials is basisDR evaluated in x= 1/t

    coeffs := [Coefficient(differentials[i],-1) : i in [1..2*g+1]];
    a := coeffs[g+1];
    vprintf LocalHeights, 1: "The residues of x^idx/(2y) at infinity are %o.\n", coeffs;
    vprintf LocalHeights, 2: "The residue of x^gdx/(2y) is  %o.\n", a;

    KX := FunctionField(X);
    x := KX.1;
    y := KX.2;

    KXp := FunctionField(Xp);
    w := KXp.1;
    z := KXp.2;

    differentialsinx := [ (x)^(i-1)*Differential(x)/(2*y) : i in [1..2*g+1]];
    holodiff := [ Pullback(rho^(-1), differentialsinx[i]) :  i in [1..g]];
    holodiffcoerced := [ KUp!(holodiff[i] / Differential(w))*Differential(u)  : i in [1..g] ]; //coerce to have the right parent
    nonholodiff:= [coeffs[i]*Pullback(rho^(-1), differentialsinx[g+1]) - a*Pullback(rho^(-1), differentialsinx[i]) : i in [g+2..2*g+1]];
    nonholodiffcoerced := [KUp!(nonholodiff[i]/Differential(w))*Differential(u) : i in [1..#nonholodiff]];

    // This is the basis of differentials in Xp coming from the canonical basis in X
    // holoDiff:=[u^(i-1)*Differential(u)/(2*v + Evaluate(hp,u)) : i in [1..g]];
    // nonholoDiff:= [coeffs[i]*u^g*Differential(u)/(2*v + Evaluate(hp,u)) - a*(u^(i-1)*Differential(u)/(2*v + Evaluate(hp,u))) : i in [g+2..2*g+1]];

    basisDR := holodiffcoerced cat nonholodiffcoerced;
    vprintf LocalHeights, 3: "The basis for the non-holomorphic differential forms that we chose is %o.\n", nonholodiff;
    matrixChange := [[i eq j select 1 else 0 :j in [1..2*g+1]]: i in [1..g]] cat [[0 : j in [1..g]] cat [coeffs[i]] cat [g+j eq i select -a else 0 : j in [2..g+1]] : i in [g+2..2*g+1]];
    return basisDR, coeffs, Matrix(matrixChange);
end intrinsic;

function GetCoefficientDerivative(poly,degree,i)
  if i eq 0 then
    return 0;
  end if;
  return i*Coefficient(poly,degree-(i-1));
end function;

function GetCoefficient(poly,degree,i)
  if degree-i lt 0 then
    return 0;
  end if;
  return Coefficient(poly,degree-i);
end function;

function GetCoefficientSmallPoly(Q,degree,i)
  // We organize the coefficients as C0 + C1x^1+...
  if degree lt i then
    return 0;
  else
    return Coefficient(Q,degree-i);
  end if;
end function;

function BoundDegree(f,D,num)
  return Max(0,Degree(num)-(Degree(f)+Degree(D))+1);
end function;

intrinsic GetMatrix(f::RngUPolElt,D::RngUPolElt,Q::RngUPolElt,num::RngUPolElt,bound::RngIntElt) -> ModMatFldElt
  {}
  // Assume the big polynomial has degree deg (this is deg(num+1) elts)
  // Assume the small polynomial has degree degf (this is degf+1 elts)
  deg := Degree(num);
  degf := Degree(f);
  R := Parent(f);
  PprimeTimes := R!(2*D*f);
  Ptimes := R!(2*Derivative(D)*f-2*(D*Derivative(Q)*f/Q)+D*Derivative(f));
  matrix := [];
  for degree in [0..deg-1] do
    vecDeg := [];
    Append(~vecDeg,GetCoefficient(Ptimes,degree,0));
    for i in [1..degree] do
      Append(~vecDeg,GetCoefficientDerivative(PprimeTimes,degree,i)+GetCoefficient(Ptimes,degree,i));
    end for;
    Append(~vecDeg,GetCoefficientDerivative(PprimeTimes,degree,degree+1));
    // We need to add one more column for P
    vecDeg := vecDeg cat [0:i in [degree+2..deg]];
    Append(~matrix, [vecDeg[i] : i in [1..bound+1]] cat [GetCoefficientSmallPoly(Q,degree,i): i in [0..degf-2]]);
  end for;
  vv := [GetCoefficientDerivative(PprimeTimes,deg,i)+GetCoefficient(Ptimes,deg,i): i in [0..bound]];
  Append(~matrix, vv cat [GetCoefficientSmallPoly(Q,deg,i): i in [0..degf-2]]);
  return Matrix(matrix);
end intrinsic;

intrinsic ReduceInCohomology(omega::DiffCrvElt, X::CrvHyp) -> DiffCrvElt, SeqEnum
  {Given a differential form omega on X, return the reduction in de Rham cohomology as a differential form,
  and as a sequence of coeffcients on the basis elements of de Rham cohomology x^i dx / (2y)}
  f, h := HyperellipticPolynomials(X);
  require h eq 0 : "Must be working on a y^2 = f model";
  Rx<x> := Parent(f);
  U := Curve(omega);
  KUp<u,v>:= FunctionField(U);
  rat := omega/Differential(u); // g(u,v) du, where this is a rational function
  KX<a,b>:= FunctionField(X);
  ratX := KX!rat;
  coefs, mons := CoefficientsAndMonomials(ratX); //A + By
  if #mons eq 2 then
    ratMin := coefs[2]*mons[2] - coefs[1]*mons[1]; // A - By (hyperellpitic involution)
  elif #mons eq 1 and mons[1] ne 1 then
    ratMin := -coefs[1]*mons[1];
  elif #mons eq 0 then
    ratMin := 0;
  else
    ratMin := coefs[1]*mons[1];
  end if;
  gOdd := (1/2)*(ratX - ratMin); // 1/2( A + By (A -By)) = By

//ratX - gOdd = d(stuff)
//y*ratX - gOdd*y = d(stuff)y
//gEven = gOdd* y*2, so
//y*ratX - gEven/2 = d(stuff)y
//ratX - gEven/(2y) = d(stuff)
//then ratX = gEven/(2y) + dstuff

  _, mons := CoefficientsAndMonomials(gOdd);
  gEven := mons[1]*gOdd*2; //this is now a function just in x, By^2

  if Index(Sprint(gEven), "/") eq 0 then
    num := eval ReplaceString(Sprint(gEven),"a","x");
    num := Rx!num;
    den := Rx!1;
  else
    num, den := Explode(Split(Sprint(gEven),"/"));
    num1 := eval num;
    den1 := eval den;
    assert num1/den1 eq gEven;
    textg2 := ReplaceString(Sprint(gEven),"a","x");
    num, den := Explode(Split(textg2,"/"));
    num := eval num;
    den := eval den;
    den := Rx!den;
    num := Rx!num;
  end if;

  if den eq Rx!1 and Degree(num) le Degree(f)-2 then
    coeffsnum := Coefficients(num);
    return omega, coeffsnum cat [0 : i in [ 1..(Degree(f)- 1) - #coeffsnum]];
  end if;

  degNum := Degree(num);
  D := Rx!(den/GCD(f*Derivative(den),den));
  bound := BoundDegree(f,D,num);
  vprintf LocalHeights, 1: "Starting with bound %o.\n", bound;
  matrix := Transpose(GetMatrix(f,D,den,num,bound));
  vec := Vector([Coefficient(num,i) : i in [0..degNum]]);
  solved := false;
  while not solved do
    try
      sol,kernel := Solution(matrix,vec);
      vprintf LocalHeights, 1 : "The dimension of the space of solutions for the reduction is %o.\n", Dimension(kernel);
      solved := true;
    catch e
      bound +:=1;
      matrix := Transpose(GetMatrix(f,D,den,num,bound));
      vprintf LocalHeights, 1: "Increasing bound to %o.\n", bound;
    end try;
  end while;
  P := &+[sol[i]*x^(i-1) : i in [1..bound+1]];
  hh := P*D/den;
  poly := &+[sol[i+bound+1]*x^(i-1) : i in [1..Degree(f)-1]];
  require num/den -( 2*f*Derivative(hh) + hh*Derivative(f) + poly) eq 0 : "Something went wrong computing the reduction.\n";
  P := &+[sol[i]*x^(i-1) : i in [1..bound+1]];
  return &+[sol[i+bound+1]*a^(i-1)/(2*b) : i in [1..Degree(f)-1]]*Differential(b), [sol[i+bound+1] : i in [1..Degree(f)-1]];
end intrinsic;

intrinsic MakeTraceZero(Z :: ModMatFldElt ) -> ModMatFldElt
  {Turn Z into a trace zero matrix.}
    tr := Integers()!Trace(Z);
    if IsZero(tr) then
      return Z;
    else
      m := LCM(tr,Ncols(Z));
      newZ := (m div tr)*Z- (m div Ncols(Z))*IdentityMatrix(Rationals(),Ncols(Z));
      assert IsZero(Trace(newZ));
      return newZ;
    end if;
end intrinsic;


intrinsic MatrixH1DR(X::CrvHyp, T::AlgMatElt) -> ModMatFldElt, SeqEnum
    {Given X hyperelliptic curve, and T an endomorphism on the H^0, return the matrix for T acting on H^1 de Rham}
    f, h := HyperellipticPolynomials(X);
    g := Genus(X);
    Rx<x> := Parent(f);
    if h ne 0 then error "h must be 0"; end if;
    Xp, rho := ChangeCoordinatesHyp(X);
    fp, hp := HyperellipticPolynomials(Xp);
    vprintf LocalHeights, 1: "Using auxiliary model %o. \n", Xp;
    KXp := FunctionField(Xp);
    P0p := AnyRationalPoint(Xp);
    Ci, Up := ConstructCorrespondenceByCantor(Xp, P0p, T);

    KUp := FunctionField(Up);
    u := KUp.1;

    KXp := FunctionField(Xp);
    w := KXp.1;

    basisDR, coeffsDR, changeOfBasis := ConstructDifferentials(X, Xp, rho, KUp);
    zo, deg := EndoAction(Xp, Ci, Up, basisDR); // This requires a y^2+h*y=f(x) model with h not 0

    KX<a,b> := FunctionField(X);
    zoX := [Pullback(rho, KXp!(zo[i]/Differential(u))*Differential(w)) : i in [1 .. #zo]];

    action := [];
    for dd in zoX do
      zomegabar, sol := ReduceInCohomology(dd, X);
      Append(~action,sol);
    end for;
    action := Transpose(Matrix(action));
    return RightInverseMatrix(Transpose(changeOfBasis))*action, deg; //deg is needed for precision bounds
end intrinsic;
