freeze;

//////////////////////////////////////////////////
// Functions for computing Frobenius structures //
//////////////////////////////////////////////////

import "auxpolys.m": log;
import "coho.m": ord_inf_mat;
import "froblift.m": frobenius, getrings, radix_reduce, reduce_mod_Q;
import "reductions.m": change_basis_b0binf, convert_to_Kxzzinvd, reduce_with_fs;
import "misc.m": compute_F, eval_R, fun_field, Kxzzinvd_to_R;


// July 21: JSM/JB added precision estimates
// 17/04/20: tried to extend this to allow for slightly worse singularities (in Tuitman's notation e.g. the denominators of the entries of W0 not being divisible by r). Not sure that it's worked.

// 16/04/20: To get this to work when the model is not smooth, we have to change how we compute basisR, phiomega, phiomegaZf and g0omega (which were originally written assuming that b^0 = [1,y,...,y^(d-1)] or equivalently that W0 is the identity. This just involves selectively multiplying by W0 or its inverse. This also means we have to do some additional reductions (i.e. quotient out by z=r/LeadingCoefficient(r). This is done using radix_reduce.

intrinsic FrobeniusStructure(data::Rec, Z::AlgMatElt, eta::ModTupFldElt, bpt::SeqEnum[FldRatElt] : N:=0)
  -> AlgMatElt[RngUPol[RngSerLaur[RngUPol[RngIntRes]]]], RngIntElt
  {Compute the matrix G of the (inverse) Frobenius structure on A_Z.}

  Q:=data`Q; v:=data`v; p:=data`p; g:=data`g; W0:=data`W0; Winf:=data`Winf; f0list:=data`f0list; finflist:=data`finflist; fendlist:=data`fendlist; 
  FH1U:=data`F; Nmax:=data`Nmax; basis:=data`basis; G0:=data`G0; Ginf:=data`Ginf; red_list_fin:=data`red_list_fin; red_list_inf:=data`red_list_inf; 
  basis:=data`basis; integrals:=data`integrals; quo_map:=data`quo_map; r:=data`r; Delta:=data`Delta; s:=data`s; frobmatb0r := data`frobmatb0r; 
  if IsZero(N) then N:=data`N; end if;

  K := BaseRing(BaseRing(Q));
  Kx := PolynomialRing(K);
  d:=Degree(Q); lc:=LeadingCoefficient(Delta);

  O, red, Ox, S, R := getrings(v, Nmax);
  // O = IntegerRing(p^Nmax), Ox = O[x], S = Ox[z,1/z], R = S[y]

  W0inv:=W0^(-1); // By assumption, W0 and W0inv are p-integral; same for LC(r)
  W:=Winf*W0^(-1);
  W0S:=ZeroMatrix(S,d,d);
  W0invS:=ZeroMatrix(S,d,d);
  for i in [1..d] do
    for j in [1..d] do
      W0n:=Numerator(W0[i,j]);
      W0d:=Denominator(W0[i,j]);
      W0invn:=Numerator(W0inv[i,j]);
      W0invd:=Denominator(W0inv[i,j]);
      if Degree(W0d) ge 1 then
        FW0:=Factorisation(W0d);
        md:=Maximum([FW0[i][2]:i in [1..#FW0]]);
        sS:=ExactQuotient((r/LeadingCoefficient(r))^md,W0d);
        W0S[i,j]:=(S.1)^(-md)*(Ox!(sS*W0n));
      else
        W0S[i,j]:=Ox!W0[i,j];
      end if;
      if Degree(W0invd) ge 1 then
        FW0inv:=Factorisation(W0invd);
        md:=Maximum([FW0inv[i][2]:i in [1..#FW0inv]]);
        sS:=ExactQuotient((r/LeadingCoefficient(r))^md,W0invd);
        W0invS[i,j]:=(S.1)^(-md)*(Ox!(sS*W0invn));
      else
        W0invS[i,j]:=ExactQuotient(Ox!W0invn,Ox!W0invd);
      end if;
    end for;
  end for;

  // Coerce Q into R:
  QR := R![Ox![red(c) : c in Eltseq(Coefficient(Q, i))] : i in [0 .. d]];

  // Coerce z = r/LeadingCoefficient(r) into R:
  zR := Ox![red(c) : c in Eltseq(r)] / red(LeadingCoefficient(r));

  // Coerce the f_i = f_{i,0}+f_{i,inf}+f_{i,end} into R:
  fRlist := [];
  finfb0list := [];
  for i:=1 to 2*g do
    fRlisti, finfb0list[i] := compute_F(Q, W0, Winf, f0list[i], finflist[i], fendlist[i]);
    fRlist[i] := Kxzzinvd_to_R(fRlisti, Q, red, r, R, W0);
    fRlist[i] := reduce_mod_Q(fRlist[i], QR, zR);
    fRlist[i] -:= eval_R(fRlist[i], bpt, r, v); // make sure that f_i(bpt) = 0 
    // No precision loss by BT17, Prop 4.5
  end for;

  // The matrix of Frobenius on H^1(X) is the 2gx2g top left corner of the matrix of Frobenius on H^1(U):
  FH1X := ChangeRing(Submatrix(FH1U, [1 .. 2*g], [1 .. 2*g]), Rationals());

  // Compute g0:
  // Z is p-integral by construction due to clearing denominators.
  //  vprintf QCMod, 2: " FH1X =  %o.\n\n", BaseRing(FH1X); 
  //  vprintf QCMod, 2: " Z =  %o.\n\n", BaseRing(Z); 
  A:=-Transpose(FH1X)*Z;
  //  vprintf QCMod, 2: " A =  %o.\n\n", BaseRing(A); 
  //  vprintf QCMod, 2: " R =  %o.\n\n", R;
  A:=ChangeRing(A,R);
  g0:=[];
  for i:=1 to 2*g do
    g0[i] := &+[A[i,j] * fRlist[j] : j in [1 .. 2*g]];
  end for;

  // Coerce s into R:
  //
  // Coerce basis elements of H^1(X) into R:
  // basis has integral coefficients by construction in basis_coho

  basisR:=[];
  for i:=1 to 2*g do
    basisR[i]:=R!0;
    for k in [1..d] do
      for j in [1..d] do
        basisR[i] +:= Ox![red(c) : c in Eltseq(basis[i][j])] * W0S[j,k]*(R.1)^(k-1);
        basisR[i] := radix_reduce(basisR[i],zR);
      end for;
    end for;
  end for;

  // At this stage, all coefficients are integral.
  // radix_reduce doesn't lead to prec loss by p.15 of [Tui17]

  // Compute g0*omega:

  g0omega := &+[g0[i] * basisR[i] : i in [1 .. 2*g]];
  g0omega:=reduce_mod_Q(g0omega,QR,zR);
  g0omega:=radix_reduce(g0omega,zR);    
  g0omega:=Vector(Coefficients(g0omega))*W0invS;
  // 16/04/20: this bit of the code is a bit odd looking: to get the right answer for the element of S^d g0omega,
  // we have to reduce using the function radix_reduce. Since this function only reduces elements of R, we
  // artificially make it an element of S^d by taking the dot product with [y^(j-1)], reducing, and taking coefficients.
  g0omega:=radix_reduce(&+[g0omega[i]*(R.1)^(i-1):i in [1..d]],zR);
  g0omega:=Vector(Coefficients(g0omega));

  // determine the number of points at infinity

  FF:=fun_field(data);
  infplaces:=InfinitePlaces(FF);
  infpoints := &+[Degree(infplace) : infplace in infplaces];

  // Compute phi^(*)(eta)-p(eta):
  kappa2 := Min([0] cat [Valuation(K!c, v) : c in Eltseq(eta)]);
  eta_new:=[];
  for i:=1 to d do
    sum:=PolynomialRing(RationalField())!0;
    for j:=1 to infpoints-1 do
      sum +:= Kx!(eta[j]*basis[2*g+j][i]); 
    end for;
    eta_new[i]:=sum;
  end for;
  eta:=eta_new;

  phi_eta:=frobenius(eta,Q,v,Nmax,r,frobmatb0r); // phi^(*)(eta)
  for i:=1 to d do
    phi_eta[i] -:= p*(S!Ox![red(c) : c in Eltseq(eta[i])]); // phi^(*)(eta)-p(eta),   as vector in S^4
  end for;

  // Compute phi^*(omega):
  // No precision loss, since basis of H^1(X) is p-integral.
  phiomega:=[];
  for i:=1 to 2*g do
    phiomegai := frobenius(basis[i],Q,v,Nmax,r,frobmatb0r);
    phiomega[i] := R!0;
    for j:=1 to d do
      for k in [1..d] do
        phiomega[i] +:= phiomegai[j]*W0S[j,k]*(R.1)^(k-1);
      end for;
    end for;
    phiomega[i]:=radix_reduce(phiomega[i],zR);
  end for;

  // Compute Z*f:
  // No precision loss, since Z is p-integral.
  Zf:=[];
  for i:=1 to 2*g do
    Zf[i] := &+[ChangeRing(Z,R)[i,j]*fRlist[j] : j in [1 .. 2*g]];
  end for;

  // Compute phi*omega*Zf:

  phiomegaZf:=R!0;
  for i:=1 to 2*g do
    phiomegaZf:=reduce_mod_Q(phiomegaZf+phiomega[i]*Zf[i],QR,zR);
  end for;
  phiomegaZf:=Vector(Coefficients(phiomegaZf))*W0invS;
  // 16/04/20: this looks a bit weird, but see remarks above on the computation of g0omega.
  phiomegaZf:=radix_reduce(&+[phiomegaZf[i]*(R.1)^(i-1):i in [1..d]],zR);
  phiomegaZf:=Vector(Coefficients(phiomegaZf));

  // Compute c and h
  chi:=convert_to_Kxzzinvd(phiomegaZf+phi_eta-g0omega,Q);
  c,h0,hinf,hend,minf:=reduce_with_fs(chi,Q,v,N,Nmax,r,W0,Winf,G0,Ginf,red_list_fin,red_list_inf,basis,integrals,quo_map); 
  hR:=Kxzzinvd_to_R(compute_F(Q,W0,Winf,h0,hinf,hend),Q,red,r,R,W0);
  hR:=reduce_mod_Q(hR,QR,zR);
  hR -:= eval_R(hR, bpt, r, v); // make sure that h(bpt) = 0 

  g1 := [g0[i] + O!c[i] : i in [1 .. 2*g]];

  // Compute matrix G:

  G:=ZeroMatrix(R,2*g+2,2*g+2);
  G[1,1]:=1;
  for i:=1 to 2*g do
    G[i+1,1]:=fRlist[i];
  end for;
  G[2*g+2,1]:=hR;
  for i:=1 to 2*g do
    for j:=1 to 2*g do
      G[i+1,j+1]:=FH1X[i,j];
    end for;
  end for;
  for i:=1 to 2*g do
    G[2*g+2,i+1]:=g1[i];
  end for;
  G[2*g+2,2*g+2]:=p;

  // Estimate precision loss 
  kappa_end := Min(&cat[&cat[[Valuation(c, v) : c in Coefficients(Eltseq(t)[i])] : i in [1..#Eltseq(t)]] : t in fendlist]);
  kappa_inf := Min(&cat[&cat[[Valuation(c, v) : c in Coefficients(Eltseq(t)[i])] : i in [1..#Eltseq(t)]] : t in finflist]);
  kappa1 := Min([kappa_inf, kappa_end]);
  kappa := Min([kappa1, kappa2, 0]);
  j := Nmax div 2;
  rhs_bound := p*Nmax - p*kappa + log(p, data`e0);
  repeat
    j +:= 1;
  until j - p*log(p, j) gt rhs_bound;
  f1 := Floor(log(p, (j-1)*data`e0)) + Floor(log(p, (-ord_inf_mat(W^-1)+1)*data`einf)); // Tui17, pf of Prop 4.10
  chiinf := change_basis_b0binf(chi,v,Nmax,r,W0,Winf);
  prec_loss_bd := Floor(log(p, -Min([Valuation(s): s in chi])*data`e0)); // Tui17, Prop. 3.7
  chiinf := change_basis_b0binf(chi,v,Nmax,r,W0,Winf);
  m_inf := -Min([Valuation(d) : d in Eltseq(chiinf)]) - Degree(data`r) + 1; 
  f2 := Floor(log(p, minf*data`einf)); // Tui17, Prop 3.8

  return G, N - Max(f1, f2);
end intrinsic;
