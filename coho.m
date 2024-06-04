freeze;

function mat_W0(Q)
  // Compute the matrix W0 using MaximalOrderFinite
  K := BaseRing(BaseRing(Q));
  Kt := RationalFunctionField(K);
  L := ext<Kt|Q>;
  b0 := Basis(MaximalOrderFinite(L));
  d := Degree(Q);
  return Matrix(Kt, d, d, [[Eltseq(L!b0[i])[j] : j in [1..d]] : i in [1..d]]);
end function;

function mat_Winf(Q)
  // Compute the matrix Winf using MaximalOrderFinite
  K := BaseRing(BaseRing(Q));
  Kt := RationalFunctionField(K);
  Kty := PolynomialRing(Kt);
  C := Coefficients(Q);
  Qnew := &+[Evaluate(C[i], 1/Kt.1)*Kty.1^(i-1) : i in [1..#C]];
  L := ext<Kt|Qnew>;
  binf := Basis(MaximalOrderFinite(L));
  d := Degree(Qnew);
  mat := Matrix(Kt, d, d, [[Eltseq(L!binf[i])[j] : j in [1..d]] : i in [1..d]]);
  return Evaluate(mat, 1/Kt.1);
end function;


function ddx(f)

  // Differentiate polynomial f(x)(y) with respect to x.

  C:=Coefficients(f);
  for i:=1 to #C do
    C[i]:=Derivative(C[i]);
  end for;
  return Parent(f)!C;
end function;


function ddx_mat(A)

  // Differentiate matrix of rational functions.

  for i:=1 to NumberOfRows(A) do
    for j:=1 to NumberOfColumns(A) do
      A[i,j]:=Derivative(A[i,j]);
    end for;
  end for;
  return A;
end function;


function ddx_vec(v)

  // Differentiate vector of rational functions.

  for i:=1 to #Eltseq(v) do
    v[i]:=Derivative(v[i]);
  end for;
  return v;
end function;


function reduce_mod_Q_exact(f,Q)

  // Eliminate powers of y >= d_x.

  while Degree(f) gt Degree(Q)-1 do
    f -:= LeadingCoefficient(f)*(Parent(f).1)^(Degree(f)-Degree(Q))*Q;
  end while;
  return f;
end function;


function polys_to_vec(polys,degx)

  // Converts a sequence of polynomials to a vector  

  dim:=#polys*(degx+1);
  v:=[];
  cnt:=1;
  for i:=1 to #polys do
    for j:=0 to degx do
      v[cnt]:=Coefficient(polys[i],j);
      cnt +:= 1;
    end for;
  end for;
  if #polys eq 0 then
    K := RationalField();
  else
    K := BaseRing(Parent(polys[1]));
  end if;
  V:=VectorSpace(K,dim);
  return V!v;
end function;


function ram(J0,Jinf)

  // Return the maximum finite and infinite ramification
  // indices, given the matrices J0, Jinf.

  d:=NumberOfRows(Jinf);
  K := BaseRing(Jinf);
  e_list:=[];
  for i:=1 to #J0 do
    for j:=1 to d do
      Append(~e_list,Denominator(RationalField()!J0[i][j,j]));
    end for;
  end for;
  if #e_list gt 0 then
    e_0:=Maximum(e_list);
  else
    e_0:=0;
  end if;
  
  e_list:=[];
  for i:=1 to d do
    Append(~e_list,Denominator(RationalField()!Jinf[i,i]));
  end for;
  e_inf:=Maximum(e_list);

  return e_0,e_inf;
end function;


function con_mat(Q,Delta,s)

  // Compute the matrix G.

  d:=Degree(Q);
  K := BaseRing(BaseRing(Q));
  Kx<x>:=RationalFunctionField(K);
  Kxy<y>:=PolynomialRing(Kx);
  
  Delta:=Kx!Delta;
  Q:=Kxy!Q;
  s:=Kxy!s;

  list:=[];
  list[1]:=Kxy!0;
  for i:=2 to d do
    list[i]:=-(i-1)*y^(i-2)*(s/Delta)*ddx(Q);
  end for;
  for i:=1 to #list do
    list[i]:=reduce_mod_Q_exact(list[i],Q);
  end for;
  G:=ZeroMatrix(Kx,d,d);
  for i:=1 to d do
    for j:=1 to d do
      G[i,j]:=Coefficient(list[i],j-1); // G acts on the right on row vectors
    end for;
  end for;
  return(G);

end function;


function jordan_inf(Ginf)

  // Compute Jordan form of residue matrix at infinity

  res_Ginf:=-Evaluate((1/Parent(Ginf[1,1]).1)*Evaluate(Ginf,1/Parent(Ginf[1,1]).1),0);
  Jinf,Tinf:=JordanForm(res_Ginf);
  return Jinf,Tinf,Tinf^(-1);  
end function;


function jordan_0(r,G0)

  // Compute Jordan forms of residue matrices at finite points
  K := BaseRing(r);
  Kx<x>:=PolynomialRing(K);
  r:=Kx!r;
  fac:=Factorization(r);
  J0:=[**];
  T0:=[**];
  T0inv:=[**];
  for i:=1 to #fac do
    if Degree(fac[i][1]) eq 1 then
      F:=K;
      s:=-Evaluate(fac[i][1],0);
    else
      F<s>:=ext<K | fac[i][1]>;
    end if;
    res_G0:=Evaluate(G0*r/Derivative(r),s);
    J,T:=JordanForm(res_G0);
    Append(~J0,J);
    Append(~T0,T);
    Append(~T0inv,T^(-1));
  end for;
  return J0,T0,T0inv;
end function;


function ord_0(f)

  // Compute ord_0(f), where f is a rational function.

  return Valuation(Numerator(f))-Valuation(Denominator(f));
end function;


function ord_0_mat(A)
  // Compute ord_0(A), where A is a matrix of rational functions.
  return Minimum([ord_0(a) : a in Eltseq(A)]);
end function;


function ord_r(f,r)

  // Compute ord_r(f), where f is a rational function

  if f eq 0 then
    return Infinity();
  end if;
  K := BaseRing(Parent(f));
  Kx<x> := PolynomialRing(K);
  r := Kx!r;
  vlist := [Valuation(Numerator(f), g[1]) - Valuation(Denominator(f), g[1]) : g in Factorization(r)];
  return Minimum(vlist);
end function;


function ord_r_mat(A,r)
  // Compute ord_r(A), where A is matrix of rational functions
  return Minimum([ord_r(a, r) : a in Eltseq(A)]);
end function;


function ord_inf(f)
  // Compute ord_inf(f), where f is a rational function.
  return -Degree(Numerator(f))+Degree(Denominator(f));
end function;


function ord_inf_mat(A)
  // Compute ord_inf(A), where A is a matrix of rational functions.
  return Minimum([ord_inf(a) : a in Eltseq(A)]);
end function;


function res_0(w,Q,r,J0,T0inv)

  // Compute res_0(\sum w_i b^0_i dx/r).

  d:=Degree(Q);
  K := BaseRing(BaseRing(Q));
  Kx<x> := PolynomialRing(K);
  r := Kx!r;
  fac:=Factorization(r);
  res_list:=[];
  for i:=1 to #fac do
    if Degree(fac[i][1]) eq 1 then
      s:=Parent(T0inv[i][1,1])!(-Coefficient(fac[i][1],0));
    else
      s:=Parent(T0inv[i][1,1]).1;
    end if;
    v:=Vector(Evaluate(w,s));
    v:=v*T0inv[i];
    for j:=1 to d do
      if J0[i][j,j] eq 0 then
        res_list:=res_list cat Eltseq(v[j]);
      end if;
    end for; 
  end for;  

  return Vector(res_list);

end function;


function val_Kttinv_d(v)
  // Compute the valuation of an element of K[t,1/t]^d.
  return Minimum([Valuation(c) : c in Eltseq(v)]);
end function;


function res_inf(w,Q,r,W0,Winf,Ginf,Jinf,Tinfinv)

  // Compute res_inf(\sum w_i b^0_i dx/r).

  d:=Degree(Q);
  K := BaseRing(BaseRing(Q));
  degr:=Degree(r);
  Kd:=RSpace(K,d);
  Kttinv<t>:=LaurentSeriesRing(K);
  Kttinvd:=RSpace(Kttinv,d);
  
  W:=Winf*W0^(-1);
  Winv:=W^(-1);
  w:=Kttinvd!Evaluate(w,t^(-1));
  w:=w*Evaluate(Winv,t^(-1));

  res_Ginf:=-Evaluate((1/Parent(Winf[1,1]).1)*Evaluate(Ginf,1/Parent(Winf[1,1]).1),0);

  // reduce to a cohomologous 1-form that is logarithmic at all points lying over x=inf:

  while val_Kttinv_d(w) lt -degr+1 do
    m:=-val_Kttinv_d(w)-degr+1;
    mat:=res_Ginf-m*IdentityMatrix(K,d);
    rhs:=Kd!0;
    for i:=1 to d do
      rhs[i]:=rhs[i]+Coefficient(-w[i],-m-degr+1)/LeadingCoefficient(r);
    end for;
    vbar:=rhs*mat^(-1);
    w:=w-ChangeRing(vbar,Kttinv)*t^(-m)*Evaluate(r*Ginf,t^(-1))-Evaluate(r,1/t)*m*t^(1-m)*ChangeRing(vbar,Kttinv);  
  end while;

  // now sum w_i b^{inf}_i dx/r is logarithmic at all points lying over x=inf

  w:=w*t^(degr-1);
  v:=Kd!0;
  for i:=1 to d do
    v[i]:=Coefficient(w[i],0);
  end for;

  // project v onto the eigenspace of res_Ginf of eigenvalue 0 

  v:=v*Tinfinv;
  res_list:=[];
  for i:=1 to d do
    if Jinf[i,i] eq 0 then
      res_list:=Append(res_list,v[i]);
    end if;
  end for;

  return Vector(res_list);
end function;


function basis_coho(Q,v,r,W0,Winf,G0,Ginf,J0,Jinf,T0inv,Tinfinv,useU,useY,basis0,basis1,basis2)

  // Compute a basis for H^1(X).
  K := BaseRing(BaseRing(Q));
  Kx<x>:=PolynomialRing(K);
  Kxy<y>:=PolynomialRing(Kx);
  d:=Degree(Q);
  Kxd:=RSpace(Kx,d);
  degr:=Degree(r);
 
  W:=Winf*W0^(-1);

  Winv:=W^(-1);
  ord0W:=ord_0_mat(W);
  ordinfW:=ord_inf_mat(W);
  ord0Winv:=ord_0_mat(Winv);
  ordinfWinv:=ord_inf_mat(Winv);

  // Compute a basis for E0

  deg_bound_E0:=degr-ord0W-ordinfW-2; 
  basisE0:=[];
  for i:=0 to d-1  do 
    for j:=0 to deg_bound_E0 do
      basisE0:=Append(basisE0,[i,j]);
    end for;
  end for;
  dimE0:=#basisE0;
  E0:=VectorSpace(K,dimE0);

  // Compute a matrix with kernel (E0 intersect Einf).

  matE0nEinf:=ZeroMatrix(K,dimE0,d*(-ordinfW-ordinfWinv));
  for i:=1 to dimE0 do
    temp:=RowSequence(x^(basisE0[i][2])*Winv)[basisE0[i][1]+1];
    for j:=0 to d-1 do
      for k:=0 to (-ordinfW-ordinfWinv-1) do
        matE0nEinf[i,j*(-ordinfW-ordinfWinv)+k+1]:=Coefficient(Numerator((Parent(W[1,1]).1)^(-ord0Winv)*temp[j+1]),k-ord0W-ord0Winv+degr-1);
      end for;
    end for;
  end for;  

  E0nEinf:=Kernel(matE0nEinf);

  // Compute a matrix with kernel the elements of E0 logarithmic at infinity.

  matlogforms:=ZeroMatrix(K,dimE0,d*(-ord0W-ordinfW-ordinfWinv-1));
  for i:=1 to dimE0 do
    temp:=RowSequence(x^(basisE0[i][2])*Winv)[basisE0[i][1]+1];
    for j:=0 to d-1 do
      for k:=0 to (-ord0W-ordinfW-ordinfWinv-2) do
        matlogforms[i,j*(-ord0W-ordinfW-ordinfWinv-1)+k+1]:=Coefficient(Numerator((Parent(W[1,1]).1)^(-ord0Winv)*temp[j+1]),k-ord0Winv+degr);
      end for;
    end for;
  end for;

  logforms:=E0nEinf meet Kernel(matlogforms);

  // Compute the finite residues.
  
  w:=Kxd!0;
  w[1]:=1;
  res0dim:=Dimension(Parent(res_0(w,Q,r,J0,T0inv)));
  matres0:=ZeroMatrix(K,dimE0,res0dim); 
  for i:=1 to dimE0 do
    w:=Kxd!0;
    w[basisE0[i][1]+1]:=x^(basisE0[i][2]);
    coefs:=res_0(w,Q,r,J0,T0inv); 
    for j:=1 to res0dim do
        matres0[i,j]:=coefs[j];
    end for;
  end for;

  // Compute the infinite residues. 

  w:=Kxd!0;
  w[1]:=1;
  resinfdim:=Dimension(Parent(res_inf(w,Q,r,W0,Winf,Ginf,Jinf,Tinfinv)));
  matresinf:=ZeroMatrix(K,dimE0,resinfdim); 
  for i:=1 to dimE0 do
    w:=Kxd!0;
    w[basisE0[i][1]+1]:=x^(basisE0[i][2]);
    coefs:=res_inf(w,Q,r,W0,Winf,Ginf,Jinf,Tinfinv); 
    for j:=1 to resinfdim do
        matresinf[i,j]:=coefs[j];
    end for;
  end for;

  forms2ndkind:=Kernel(matres0) meet Kernel(matresinf);
  cocyc:=E0nEinf meet forms2ndkind;

  "forms2ndkind", Dimension(forms2ndkind);
  "E0nEinf", Dimension(E0nEinf);
  forms1stkind:=logforms meet forms2ndkind;
  
  // Compute a matrix with kernel (B0 intersect Binf)

  deg_bound_B0:=-ord0W-ordinfW-1;
  basisB0:=[];
  for i:=0 to d-1  do 
    for j:=0 to deg_bound_B0 do
      basisB0:=Append(basisB0,[i,j]);
    end for;
  end for;
  dimB0:=#basisB0;
  B0:=VectorSpace(K,dimB0);

  matB0nBinf:=ZeroMatrix(K,dimB0,d*(-ordinfW-ordinfWinv));
  for i:=1 to dimB0 do
    power_x:=basisB0[i][2];
    power_y:=basisB0[i][1];
    temp:=RowSequence(x^(power_x)*Winv)[power_y+1];
    for j:=0 to d-1 do
      for k:=0 to (-ordinfW-ordinfWinv-1) do
        matB0nBinf[i,j*(-ordinfW-ordinfWinv)+k+1]:=Coefficient(Numerator((Parent(W[1,1]).1)^(-ord0Winv)*temp[j+1]),k-ord0W-ord0Winv);
      end for;
    end for;
  end for;

  // Compute d(B0 intersect Binf).

  B0nBinf:=Kernel(matB0nBinf);
  basisB0nBinf:=Basis(B0nBinf);
  dimB0nBinf:=#basisB0nBinf;  

  list:=[];
  for i:=1 to dimB0nBinf do
    vec:=basisB0nBinf[i];
    vecQxd:=Kxd!0;
    for j:=1 to dimB0 do
      vecQxd[basisB0[j][1]+1]:=vecQxd[basisB0[j][1]+1]+vec[j]*x^(basisB0[j][2]);
    end for;
    vecQxd:=vecQxd*ChangeRing(r*G0,Kx)+r*ddx_vec(vecQxd);
    coefs:=[];
    for j:=1 to dimE0 do
      power_x:=basisE0[j][2];
      power_y:=basisE0[j][1];
      coefs[j]:=Coefficient(vecQxd[power_y+1],power_x);  
    end for;
    list:=Append(list,E0!coefs);
  end for;
  matd:=Matrix(list);

  // Compute bases

  "cocyc", Dimension(cocyc);
  "list", #list;
  cobound:=sub<cocyc|list>;

  if basis0 eq [] then
    b0:=Basis(forms1stkind);
  else
    b0:=[];
    for i:=1 to #basis0 do
      b0[i]:=polys_to_vec(basis0[i],deg_bound_E0);
    end for;
  end if;  

  b5:=[];
  for i:=2 to dimB0nBinf do
    b5:=Append(b5,E0!list[i]);
  end for;
  
  dualspace:=Complement(cocyc,forms1stkind+cobound); // Take the dual w.r.t. cup product? Right now just any complement of forms of the 1st kind in H^1(X).
  if basis1 eq [] then
    b1:=Basis(dualspace);
  else
    b1:=[];
    for i:=1 to #basis1 do
      b1[i]:=polys_to_vec(basis1[i],deg_bound_E0);
    end for;
  end if;  

  basisH1X:=b0 cat b1;
  dimH1X:=#basisH1X;

  finiteregularlogarithmic:=logforms meet Kernel(matres0); // 1-forms that generate H^1(Y) over H^1(X), where Y=X-x^{-1}(\infty)
  H1YmodH1X:=Complement(finiteregularlogarithmic,forms1stkind);

  if basis2 eq [] then
    b2:=Basis(H1YmodH1X);
  else
    b2:=[];
    for i:=1 to #basis2 do
      b2[i]:=polys_to_vec(basis2[i],deg_bound_E0);
    end for;
  end if;

  b3:=Basis(Complement(E0nEinf,cocyc+H1YmodH1X));

  b4:=Basis(Complement(E0,E0nEinf));  
  
  b:=b0 cat b1 cat b2 cat b3 cat b4 cat b5;

  dimH1U:=#b0+#b1+#b2+#b3;
  dimH1Y:=#b0+#b1+#b2;

  if useU then
    dim:=dimH1U;
  elif useY then
    dim:=dimH1Y;
  else
    dim:=dimH1X;
  end if;

  p := Factorization(Norm(v))[1][1];
  vp := Valuation(K!p, v);
  for i:=1 to dim do
    valdenom := Minimum([Valuation(b[i][j], v)/vp : j in [1 .. dimE0]]);
    if valdenom lt 0 then
      b[i] *:= p^(Ceiling(-valdenom));
    end if;
  end for; 

  matb:=Matrix(b);
  quo_map:=matb^(-1);

  integrals:=[Kxd|];
  for i:=2 to dimB0nBinf do
    vec:=Kxd!0;
    for j:=1 to dimB0 do
      vec[basisB0[j][1]+1] +:= (basisB0nBinf[i][j])*x^(basisB0[j][2]);
    end for;
    Append(~integrals,LeadingCoefficient(r)*vec); // factor lc(r) here, since working with dx/z basis instead of dx/r
  end for;

  Kxd:=RSpace(Kx,d);
  basis:=[Kxd|];
  
  for i:=1 to dim do
    vec:=Kxd!0;
    for j:=1 to dimE0 do
      vec[basisE0[j][1]+1] +:= (K!(b[i][j]))*(Kx.1)^(basisE0[j][2]);
    end for;
    Append(~basis,vec);
  end for;

  return basis,integrals,quo_map;
    
end function;
