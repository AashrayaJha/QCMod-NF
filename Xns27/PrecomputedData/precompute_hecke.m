//This file collects Coleman and Hecke data, for one prime of the number field at a time. 
//The main intrinsic is ColemanData (in singleintegral.m), which in particular computes 
//Frobenius action on H^1_dR(X), for which we use a symplectic basis.

AttachSpec("QCMod.spec");
import "auxpolys.m": auxpolys, log;
import "singleintegrals.m": evalf0, is_bad, local_coord, set_point, tadicprec, teichmueller_pt, xy_coordinates;
import "misc.m": are_congruent, equivariant_splitting, eval_mat_R, eval_Q, FindQpointQp, fun_field, alg_approx_Qp, minprec, minval, minvalp;
import "applications.m": Q_points, Qp_points, roots_with_prec, separate;
import "heights.m": E1_tensor_E2, expand_algebraic_function, frob_equiv_iso, height;

K<u> := CyclotomicField(3);
OK := Integers(K);
Kx<x> := PolynomialRing(K);
Kxy<y> := PolynomialRing(Kx);
Q:= y^4 + ((-2*u + 9)*x + (2*u + 3))*y^3 + (-3*x^2 + 6*x - 3)*y^2 + ((-170*u + 254)*x^3 + (-150*u + 114)*x^2 + (-54*u + 18)*x - 10*u - 2)*y + (162*u + 144)*x^4 + (-108*u + 48)*x^3 + (-72*u - 144)*x^2 + (12*u - 48)*x + 6*u;
p := 13;
v := Factorization(p*OK)[1][1];
N := 400;

t1 := Cputime();
"Constructing symplectic basis of H1...";
h1basis, g, r, W0 := H1Basis(Q, v);
prec := 2*g;
// jsm: changed split to false, otherwise this takes forever
cpm := CupProductMatrix(h1basis, Q, g, r, W0 : prec := prec, split := false);
sympl := SymplecticBasisH1(cpm);
new_complementary_basis := [&+[sympl[i,j]*h1basis[j] : j in [1..2*g]] : i in [1..g]];
sympl_basis := [h1basis[i] : i in [1..g]] cat new_complementary_basis;
basis0 := [[sympl_basis[i,j] : j in [1..Degree(Q)]] : i in [1..g]]; // basis of regular differentials
basis1 := [[sympl_basis[i,j] : j in [1..Degree(Q)]] : i in [g+1..2*g]];  // basis of complementary subspace
"Symplectic basis constructed.";
t2 := Cputime();
t2-t1;

t1 := Cputime();
"Constructing Coleman data...";
time data := ColemanData(Q, v, N : useU:=false,  basis0:=basis0, basis1:=basis1);
"Coleman data constructed.";
t2 := Cputime();
t2 - t1;

d:=Degree(Q); 
q:=p; 
Qp:=pAdicField(p,N);

K:=BaseRing(BaseRing(Q));

F1 := data`F;
if q eq p then F1 := Submatrix(data`F,1,1,2*g,2*g); end if;// Necessary when q=p
F1inv := Transpose(F1)^(-1);
Aq_1 := Transpose(F1)+q*F1inv;   // Eichler-Shimura -> Hecke operator
// bound for precision loss in Eichler-Shimura
prec_loss_bd := Valuation(Norm(Determinant(F1inv)), p)+1;

AK := ZeroMatrix(K, 2*g, 2*g); 
bad_indices:=[**];

//Hecke correspondence
for j in [1..2*g] do
    for k in [1..2*g] do
        try
            AK[j,k] := alg_approx_Qp(Qp!Rationals()!Aq_1[j,k],v);    // recognition of integer in Zp via LLL
        catch e
            Append(~bad_indices,[j,k]);
        end try;                   
    end for;
end for;
//Finish Hecke correpondence. 
print "bad_indices", bad_indices;

// Start computing nice correspondences 
C:=ZeroMatrix(K,2*g,2*g);
for i:=1 to g do
C[i,g+i]:=1;
C[g+i,i]:=-1; 
end for;
Z1:=ZeroMatrix(K,2*g,2*g); Z2:=ZeroMatrix(K,2*g,2*g);
Zs:=[Z1,Z2];

for i in [1..2] do 
    A:=Aq_1^i;
    Zmx := (2*g*Aq_1^i-Trace(A)*IdentityMatrix(K,2*g))*C^(-1);
    for j in [1..2*g] do
        for k in [1..2*g] do
            try
                Zs[i][j,k] := alg_approx_Qp(Qp!Rationals()!Zmx[j,k],v);    // recognition of integer in Zp via LLL
            catch e
                Append(~bad_indices,[j,k]);
            end try;                   
        end for;
    end for;
end for;

//print results

output_file:="data/Hecke_good_patch_400.m";

out_AK:=Sprintf("AK_good_patch:=%m;",AK);
out_AK := ReplaceString(out_AK, "! ", " ");

out_Zs :=Sprintf("Zs_good_patch:=%m;", Zs) ;
out_Zs := ReplaceString(out_Zs, "! ", " ");

out_F := Sprintf("F:=%m;", data`F);
out_Q:=Sprintf("Q := %m;", data`Q);
out_N:=Sprintf("N:=%m;", N);
out_prec_loss_bd:=Sprintf("prec_loss_bd:=%m;", prec_loss_bd);

Write(output_file,out_Q);
Write(output_file,out_N);
Write(output_file,out_AK);
Write(output_file,out_Zs);
Write(output_file,out_F);
Write(output_file,out_prec_loss_bd);

"Coleman and Hecke data recorded.";



