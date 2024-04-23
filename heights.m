freeze;

/////////////////////////////////////
// functions for computing heights //
/////////////////////////////////////

// July 21: JSM/JB added precision estimates

import "auxpolys.m": log;
import "singleintegrals.m": coleman_integrals_on_basis_divisors, is_bad, lie_in_same_disk,
  local_coord, tiny_integral_prec, eval_poly_Qp, eval_ff_mat_Qp;
import "misc.m": fun_field, minprec;
import "applications.m": are_equal_records;

function frob_equiv_iso(G,data,N)

  // Given a matrix G defining Frob^-1 on a mixed extension this returns s^phi.

  p:=data`p; v:=data`v; g:=data`g; F:=data`F; //N:=data`N; 
  K := NumberField(Order(v));
  Kv, loc := Completion(K, v);

  Qp:=pAdicField(p,N);  
  FQp:=ZeroMatrix(Qp,2*g,2*g);
  for i:=1 to 2*g do
    for j:=1 to 2*g do
      FQp[i,j]:=Qp!loc(F[i,j]);
    end for;
  end for;

  fx:=Matrix(Qp,2*g,1,[Qp!G[i+1,1]:i in [1..2*g]]);
  gx:=Matrix(Qp,1,2*g,[Qp!G[2*g+2,i+1]:i in [1..2*g]]);
  hx:=Matrix(Qp,1,1,[Qp!G[2*g+2,1]]);
  I := IdentityMatrix(Qp,2*g);
  
  hxprime := 1/(1-p)*(gx*((I-FQp)^-1)*fx+hx);
  
  s_phi := IdentityMatrix(Qp,2*g+2);
  for i in [1..2*g] do
    s_phi[i+1,1] := ((I-FQp)^(-1)*fx)[i,1];
    s_phi[2*g+2,i+1] := (gx*(FQp-p*I)^-1)[1,i];
  end for;
  s_phi[2*g+2,1] := hxprime[1,1];

  return s_phi, Min(N, minprec(s_phi));

end function;


intrinsic ParallelTransport(P1::Rec, P2::Rec, Z::AlgMatElt, eta::ModTupFldElt, data::Rec :
    prec := 0, N := 0)
  -> AlgMatElt[FldPad], RngIntElt
  {Computes the parallel transport map of the unipotent connection Lambda
  defined by Z, eta from P1 to P2.}

  // TODO: on a finite bad disk, make sure to take local parameters centered at the (very) bad point.

  x1:=P1`x; x2:=P2`x; b1:=P1`b; b2:=P2`b; Q:=data`Q; p:=data`p; v:=data`v; W0:=data`W0; Winf:=data`Winf;
  r:=data`r; basis:=data`basis; g:=data`g;
  d:=Degree(Q); Qp:=Parent(x1);
  K := BaseRing(BaseRing(Q));
  Kv, loc := Completion(K, v: Precision:=N);
  lc_r := Qp!loc(LeadingCoefficient(r));

  if IsZero(N) then N:=data`N; end if;

  C:=IdentityMatrix(Qp,2*g+2);
  if are_equal_records(P1, P2) then
    return C,data`N;
  end if;
    
  if not lie_in_same_disk(P1,P2,data) then
    error "the points do not lie in the same residue disk";
  end if;

  xt,bt,index:=local_coord(P1,prec,data); 

  if index eq 0 then // x-x(P1) is the local coordinate
    val:=Valuation(x2-x1);
  else // b[index]-b[index](P1) is the local coordinate
    val:=Valuation(b2[index]-b1[index]);
  end if;

  Qt<t> := LaurentSeriesRing(RationalField(),prec);
  xt := Qt!xt;
  ChangeUniverse(~bt, Qt);
  bt := Vector(bt);

  Qpt:=LaurentSeriesRing(Qp,prec);
  Zp:=RingOfIntegers(Qp);
  Zpt:=LaurentSeriesRing(Zp,prec);

  if P1`inf then
    xt:=1/xt;
    xt:=Qt!Qpt!xt; 
    Winv:=W0*Winf^(-1);
    bt := bt * Transpose(eval_ff_mat_Qp(Winv, xt, v, N));
    for i:=1 to d do
      bt[i]:=Qt!(Qpt!bt[i]);
    end for; 
  end if;

  if P1`inf or not is_bad(P1,data) then 
    denom:=Qt!Qpt!(1/eval_poly_Qp(r, xt, v, N));
  else
    Qp_N := pAdicField(p,N);
    Qpx := PolynomialRing(Qp_N);
    rQp := eval_poly_Qp(r, Qpx.1, v, N);
    zero:=HenselLift(rQp,x1);
    sQp:=rQp div (Qpx.1-zero);
    denom:=Qt!Qpt!((Qt!Zpt!(xt-Coefficient(xt,0)))^(-1)*(Qt!Qpt!(1/eval_poly_Qp(sQp, xt, v, N))));
  end if;

  // determine the number of points at infinity

  FF := fun_field(data);
  infpoints := &+[Degree(place) : place in InfinitePlaces(FF)];

  // compute the powerseries expansions of the basis elements of H^1(Y)

  derxt:=Qt!Qpt!Derivative(xt); 
  omegax:=[];
  for i:=1 to 2*g+infpoints-1 do
    basisxt := Vector([eval_poly_Qp(c, xt, v, N) : c in Eltseq(basis[i])]);
    for j:=1 to d do
      basisxt[j] := Qt!Qpt!basisxt[j];
    end for;
    omegax[i]:=InnerProduct(Vector(basisxt*derxt*lc_r*denom),bt);
    omegax[i]:=Qt!Qpt!omegax[i];
    assert IsZero(Coefficient(omegax[i],-1)); // second kind
  end for;

  // tiny single integrals
  //
  // compute p-adic precision
  e:=Degree(Parent(x2)); // e=1 in our applications
  maxpoleorder:=-(Minimum([Valuation(omegax[i]): i in [1..2*g]])); 
  maxdegree:=Maximum([Degree(omegax[i]): i in [1..2*g]]);
  mindegree:=Minimum([Degree(omegax[i]): i in [1..2*g]]);
  Nsingle := tiny_integral_prec(prec,1,maxpoleorder,maxdegree,mindegree,val,data);

  Omegax:=[];
  for i:=1 to 2*g do
    Omegax[i]:=-Integral(omegax[i]);
  end for;

  singleintegrals:=[];
  for i:=1 to 2*g do
    if index eq 0 then // x-x(P1) is the local coordinate
      singleintegrals[i]:=Evaluate(Omegax[i],x2-x1);
    else // b[index]-b[index](P1) is the local coordinate
      singleintegrals[i]:=Evaluate(Omegax[i],b2[index]-b1[index]);
    end if;
  end for;
 
  // tiny double integral

  dgx:=0;
  for i:=1 to 2*g do
    for j:=1 to 2*g do
      dgx +:= omegax[i]*(Qp!Z[i,j])*Omegax[j]; 
    end for;
  end for;
  for i:=1 to infpoints-1 do
    dgx -:= Qt!Qp!loc(eta[i])*omegax[2*g+i];
  end for;

  val_eta := Min([0] cat [Valuation(Qp!loc(c)) : c in Eltseq(eta)]);
  maxpoleorder:=-Valuation(dgx);
  maxdegree:=Degree(dgx);
  mindegree:=maxdegree;
  Ndouble := Integers()!tiny_integral_prec(prec,1,maxpoleorder,maxdegree,mindegree,val,data : N := N+Min(val_eta, -Floor(log(p,N))));

  gx:=Integral(dgx);

  if index eq 0 then // x-x(P1) is the local coordinate
    doubleintegral:=Evaluate(gx,x2-x1);
  else // b[index]-b[index](P1) is the local coordinate
    doubleintegral:=Evaluate(gx,b2[index]-b1[index]);
  end if;

  C:=IdentityMatrix(pAdicField(p,Ndouble), 2*g+2);
  // entries in first column (except last one) 
  for i:=1 to 2*g do
    C[i+1,1]:=-singleintegrals[i];
  end for;
  // entries in the last row (except first one)
  for i:=1 to 2*g do
    for j:=1 to 2*g do
      C[2*g+2,i+1]:=C[2*g+2,i+1]-singleintegrals[j]*(Qp!Z[j,i]);
    end for;
  end for;
  // bottom left entry
  C[2*g+2,1]:=-doubleintegral;

  return C, Ndouble;

end intrinsic;


intrinsic ParallelTransportToZ(P::Rec, Z::AlgMatElt, eta::ModTupFldElt, data::Rec :
    prec := 0, N := 0)
  -> AlgMatElt[FldPad], RngSerLaurElt, ModTupFldElt[RngSerLaur]
  {Computes the parallel transport map of the unipotent connection Lambda
  defined by Z, eta from P to to an arbitrary point in its residue disk as a
  power series in the local parameter there. The series expansions xt and bt
  of the coordinates on the curve in terms of this local parameter are also
  returned.}

  x0:=P`x; b:=P`b; Q:=data`Q; p:=data`p; v:=data`v; basis:=data`basis; g:=data`g; r:=data`r;
  d:=Degree(Q); W0:=data`W0; Winf:=data`Winf; W:=Winf*W0^(-1); Qp:=Parent(x0);
  K := BaseRing(BaseRing(Q));
  Kv, loc := Completion(K, v);
  lc_r := Qp!loc(LeadingCoefficient(r));

  if IsZero(N) then N:=data`N; end if;

  P1:=P; // TODO: on a finite bad disk, make sure to take local parameters centered at the (very) bad point.
  x1:=P1`x;

  xt,bt,index:=local_coord(P1,prec,data); 

  xtold:=xt;
  btold:=bt;

  Qt<t>:=LaurentSeriesRing(RationalField(),prec);
  xt:=Qt!xt;
  ChangeUniverse(~bt, Qt);
  bt:=Vector(bt);

  Qpt:=LaurentSeriesRing(Qp,prec);
  Zp:=RingOfIntegers(Qp);
  Zpt:=LaurentSeriesRing(Zp,prec);

  if P1`inf then
    xt:=1/xt;
    xt:=Qt!Qpt!xt; 
    Winv:=W0*Winf^(-1);
    bt:=bt*Transpose(eval_ff_mat_Qp(Winv, xt, v, N));
    for i:=1 to d do
      bt[i]:=Qt!(Qpt!bt[i]);
    end for; 
  end if;

  if P1`inf or not is_bad(P1,data) then 
    denom:=Qt!Qpt!(1/eval_poly_Qp(r, xt, v, N));
  else
    Qp_N := pAdicField(p,N);
    Qpx := PolynomialRing(Qp_N);
    rQp := eval_poly_Qp(r, Qpx.1, v, N);
    zero:=HenselLift(rQp,x1);
    sQp:=rQp div (Qpx.1-zero);
    denom:=Qt!Qpt!((Qt!Zpt!(xt-Coefficient(xt,0)))^(-1)*(Qt!Qpt!(1/eval_poly_Qp(sQp, xt, v, N))));
  end if;

  // determine the number of points at infinity

  FF:=fun_field(data);
  infpoints := 0;
  for place in InfinitePlaces(FF) do
    infpoints +:= Degree(place);
  end for;

  // compute the powerseries expansions of the basis elements of H^1(Y)
  // TODO: These get recomputed a lot. Maybe store them and update if we need them to
  // larger precision?

  derxt:=Qt!Qpt!Derivative(xt); 
  omegax:=[];
  for i:=1 to 2*g+infpoints-1 do
    basisxt := Vector([eval_poly_Qp(c, xt, v, N) : c in Eltseq(basis[i])]);
    for j:=1 to d do
      basisxt[j] := Qt!Qpt!basisxt[j];
    end for;
    omegax[i]:=InnerProduct(Vector(basisxt*derxt*lc_r*denom),bt);
    omegax[i]:=Qt!Qpt!omegax[i];
    assert IsZero(Coefficient(omegax[i],-1)); // second kind
  end for;

  // tiny single integrals
  Omegax:=[];
  for i:=1 to 2*g do
    Omegax[i]:=-Integral(omegax[i]);
  end for;

  // tiny double integral
  dgx:=0;
  for i:=1 to 2*g do
    for j:=1 to 2*g do
      dgx +:= omegax[i]*(Qp!Z[i,j])*Omegax[j]; 
    end for;
  end for;
  for i:=1 to infpoints-1 do
    dgx -:= Qt!Qp!loc(eta[i])*omegax[2*g+i];
  end for;
  gx:=Integral(dgx);

  C:=IdentityMatrix(Parent(gx),2*g+2);
  // entries in first column (except last one) 
  for i:=1 to 2*g do
    C[i+1,1]:=-Omegax[i];
  end for;
  // entries in the last row (except first one)
  for i:=1 to 2*g do
    for j:=1 to 2*g do
      C[2*g+2,i+1] -:= Omegax[j]*(Qp!Z[j,i]);
    end for;
  end for;
  // bottom left entry
  C[2*g+2,1]:=-gx;

  xt:=xtold;
  bt:=Vector(btold);

  return C,xt,bt;

end intrinsic;


height := function(Phi,betafil,gammafil,splitting,data)

  // This function computes the local p-adic height at p of a 
  // filtered phi-module given its Frobenius 
  // matrix Phi and splitting of the Hodge filtration determined by gamma_fil, beta_fil.


  p:=data`p; g:=data`g; v:=data`v;

  S:=BaseRing(Phi);
  splitting:=ChangeRing(splitting,S);

  betafil    := Vector(S,[0 : i in [1..g]] cat Eltseq(betafil));
  gammafil   := S!gammafil;
  alpha1g    := Vector(S,g,[Phi[i+1,1]:i in [1..g]]);
  alpha      := Vector(S,2*g,[Phi[i+1,1]:i in [1..2*g]]);
  s1alphaphi := ChangeRing(alpha1g*Transpose(splitting),S);
  s2alphaphi := alpha-alpha1g*Transpose(splitting);
  gammaphi   := Phi[2*g+2,1];
  betaphi    := Vector(S,2*g,[Phi[2*g+2,i+1]:i in [1..2*g]]);
  
  return gammaphi-gammafil-DotProduct(s1alphaphi,betaphi)-DotProduct(s2alphaphi,betafil);

end function;


E1_NF := function(Phips,changebases,Salpha)
  
  //vprintf QCMod, 2: " Parent(Phips[1])=%o,\n", Parent(Phips[1]);

  changebases := [ChangeRing(M, Salpha) : M in changebases];
  g := (Nrows(Phips[1]) div 2) - 1;
  assert g eq 3; // TODO: change
  d := #changebases;
  assert d eq 2; // only quadratic fields at this point
  E1s := [Vector(Salpha,[Phips[j][i,1] : i in [2..g+1]])*changebases[j] : j in [1..d] ]; 
  alpha := Salpha.1;
  // E1s = [E1(sigma1(P)), E1(sigma2(P))], entries are in Qp^g
  E1s_Salpha := [&+[E1s[j][i]*alpha^(i-1) : i in [1..g]] : j in [1..d] ];
  // E1s_Salpha = [E1(sigma1(P)), E1(sigma2(P))], entries are in Salpha
  return E1s_Salpha;
end function;


E2_NF := function(Phips,betafils,changebases,Salpha)
  changebases := [ChangeRing(M, Salpha) : M in changebases];
  g := (Nrows(Phips[1]) div 2) - 1;
  assert g eq 3; // TODO: change
  d := #betafils;
  assert d eq 2; // only quadratic fields at this point
  E2s := [Vector(Salpha,[Phips[i][2*g+2,g+1+j] - betafils[i][j] : j in [1..g]])*changebases[i] : i in [1..d]];
  alpha := Salpha.1;
  // E2s = [E2(sigma1(P)), E2(sigma2(P))], entries are in Qp^g
  E2s_Salpha := [&+[E2s[j][i]*alpha^(i-1) : i in [1..g]] : j in [1..d] ];
  // E2s_Salpha = [E2(sigma1(P)), E2(sigma2(P))], entries are in Salpha
  return E2s_Salpha;
end function;



E1_tensor_E2_NF := function(E1s, E2s)

/*
  changebases := [ChangeRing(M, Salpha) : M in changebases];
  g := (Nrows(Phips[1]) div 2) - 1;
  assert g eq 3; // TODO: change
  d := #betafils;
  assert d eq 2; // only quadratic fields at this point
  E1s := [Vector(Salpha,[Phips[j][i,1] : i in [2..g+1]])*changebases[j] : j in [1..d] ]; 
  E2s := [Vector(Salpha,[Phips[i][2*g+2,g+1+j] - betafils[i][j] : j in [1..g]])*changebases[i] : i in [1..d]];
  alpha := Salpha.1;
  // E1s = [E1(sigma1(P)), E1(sigma2(P))], entries are in Qp^g
  // E2s = [E2(sigma1(P)), E2(sigma2(P))], entries are in Qp^g
  E1s_Salpha := [&+[E1s[j][i]*alpha^(i-1) : i in [1..g]] : j in [1..d] ];
  E2s_Salpha := [&+[E2s[j][i]*alpha^(i-1) : i in [1..g]] : j in [1..d] ];
  // E1s_Salpha = [E1(sigma1(P)), E1(sigma2(P))], entries are in Salpha
  // E2s_Salpha = [E2(sigma1(P)), E2(sigma2(P))], entries are in Salpha
*/
  E1_E2s_Salpha := [];
  for i := 1 to #E1s do
    for j := 1 to #E2s do
      Append(~E1_E2s_Salpha, E1s[i]*E2s[j]);
    end for;
  end for;
  // E1_E2s_Salpha = [E1(sigma1(P)\otimes E2(sigma1(P)), E1(sigma1(P)\otimes E2(sigma2(P)), E1(sigma2(P)\otimes E2(sigma1(P)), E1(sigma2(P)\otimes E2(sigma2(P))], entries are elements of S_alpha
  return E1_E2s_Salpha;
end function;


E1_tensor_E2:=function(Phi,betafil,changebasis,data,Salpha)

  changebasis:=ChangeRing(changebasis,Salpha);
  g:=data`g;
  E1 := Vector(Salpha,[Phi[i,1] : i in [2..g+1]])*changebasis; 
  E2 := Vector(Salpha,[Phi[2*g+2,g+1+j] - betafil[j] : j in [1..g]])*changebasis;

  return &+[E1[i]*Salpha.1^(i-1) : i in [1..g]] * &+[E2[i]*Salpha.1^(i-1) : i in [1..g]], E1, E2;

end function;


expand_algebraic_function:=function(P,g,data,N,prec)

  // expands the algebraic function g with respect to the chosen parameter at P.
  // the parameter is the same as in ParallelTransportToZ.

  p := data`p;
  v := data`v;
  xt,bt,index:=local_coord(P,prec,data); 

  Qpt<t>:=LaurentSeriesRing(pAdicField(p,N),prec);
  xt:=Qpt!xt;
  bt:=[Qpt!bt[i]:i in [1..#bt]];
  return &+[eval_poly_Qp(g[i], xt, v, N)*bt[i]:i in [1..NumberOfColumns(g)]];
end function;



function height_coefficients(divisors, intersections, local_CG_heights_p, data)
  // TODO: arbitrary g

  local_CG_heights_away := [ ];
  p := data`p;
  g := data`g;
  if g ne 2 then
    error "Only genus 2 currently implemented";
  end if;
  Qp := pAdicField(p, data`N);

  for i := 1 to #intersections do
   local_CG_heights_away[i] := IsEmpty(intersections[i]) select 0
           else &+[Log(Qp!t[1])*t[2] : t in intersections[i] | t[1] ne p];
  end for;

  heights_basis := [local_CG_heights_p[i] + local_CG_heights_away[i] : i in [1..#local_CG_heights_p]];

  n := g*(g+1) div 2;
  M := ZeroMatrix(Parent(divisors[1,1,1,1]), n, n);

  integrals := [coleman_integrals_on_basis_divisors(t[2], t[1], data) : t in divisors];
  M[1,1] := integrals[1,1]*integrals[1,1];
  M[1,2] := 1/2*(integrals[1,1]*integrals[1,2]+integrals[1,1]*integrals[1,2]);
  M[1,3] := integrals[1,2]*integrals[1,2];
  M[2,1] := integrals[1,1]*integrals[2,1];
  M[2,2] := 1/2*(integrals[1,1]*integrals[2,2]+integrals[2,1]*integrals[1,2]);
  M[2,3] := integrals[1,2]*integrals[2,2];
  M[3,1] := integrals[2,1]*integrals[2,1];
  M[3,2] := 1/2*(integrals[2,1]*integrals[2,2]+integrals[2,1]*integrals[2,2]);
  M[3,3] := integrals[2,2]*integrals[2,2];

  assert Determinant(M) ne 0;
  heights_vector := Matrix(BaseRing(M),n, 1, heights_basis);

  alphas := Eltseq(M^-1*heights_vector);
  return [[alphas[1], alphas[2]], [0, alphas[3]]];
end function;



function new_height_coefficients(divisors, intersections, factors, local_CG_heights_p, data)
  // TODO: arbitrary g
  // Same as height_coefficients, but allow more general divisors.


  g := data`g;
  if g ne 2 then
    error "Only genus 2 currently implemented";
  end if;

  p := data`p;
  Qp := pAdicField(p, data`N);

  local_CG_heights_away := [ ];
  for i := 1 to #intersections do
   local_CG_heights_away[i] := IsEmpty(intersections[i]) select 0
           else &+[Log(Qp!t[1])*t[2] : t in intersections[i] | t[1] ne p];
  end for;

  factors1, factors2 := Explode(factors);
  heights_basis := [local_CG_heights_p[i] - local_CG_heights_away[i] : i in [1..#local_CG_heights_p]];
  h1 := heights_basis[1]/(factors1[1,1]*factors2[1,1]);
  h2 := heights_basis[2]/(factors1[1,2]*factors2[1,2]);
  h3 := heights_basis[3]/(factors1[2,2]*factors2[2,2]);
  regulator := h1*h3-h2^2;

  n := g*(g+1) div 2;
  M := ZeroMatrix(Parent(divisors[1,1,1,1]), n, n);

  
  integrals := [coleman_integrals_on_basis_divisors(t[2], t[1], data) : t in divisors];
  integrals[1,1] /:= factors1[1,1];
  integrals[1,2] /:= factors1[1,1];
  integrals[2,1] /:= factors2[2,2];
  integrals[2,2] /:= factors2[2,2];

  M[1,1] := integrals[1,1]*integrals[1,1];
  M[1,2] := 1/2*(integrals[1,1]*integrals[1,2]+integrals[1,1]*integrals[1,2]);
  M[1,3] := integrals[1,2]*integrals[1,2];
  M[2,1] := integrals[1,1]*integrals[2,1];
  M[2,2] := 1/2*(integrals[1,1]*integrals[2,2]+integrals[2,1]*integrals[1,2]);
  M[2,3] := integrals[1,2]*integrals[2,2];
  M[3,1] := integrals[2,1]*integrals[2,1];
  M[3,2] := 1/2*(integrals[2,1]*integrals[2,2]+integrals[2,1]*integrals[2,2]);
  M[3,3] := integrals[2,2]*integrals[2,2];

  assert Determinant(M) ne 0;  // Points should be independent
  heights_vector := Matrix(BaseRing(M),n, 1, [h1,h2,h3]);

  alphas := Eltseq(M^-1*heights_vector);
  return [[alphas[1], alphas[2]], [0, alphas[3]]], regulator;
end function;

