//This file collects Coleman data, for one prime of the number field at a time. 
//The main intrinsic is ColemanData (in singleintegral.m), which in particular computes 
//Frobenius action on H^1_dR(X).

//AttachSpec("CHIMP/CHIMP.spec");
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
//Q := y^4 + ((-14*u + 9)*x - u + 1)*y^3 - 3*x^2*y^2 + ((-2330*u - 2662)*x^3 + (-663*u - 726)*x^2 + (-63*u - 66)*x - 2*u - 2)*y + (-978*u + 12)*x^4 + (-298*u + 1)*x^3 - 30*u*x^2 - u*x;
p := 13;
v := Factorization(p*OK)[2][1];
N := 201;

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

t1 := Cputime();
"Recording Coleman data...";
output_file := "data/NF-example-coleman-data-13_patch2_201.m";
fprintf output_file, "K<u> := CyclotomicField(3);\n";
fprintf output_file, "_<x> := PolynomialRing(K);\n";
fprintf output_file, "_<z> := LaurentSeriesRing(PolynomialRing(K));\n\n";

out := Sprintf("data := %m;\n\n", data);
out := ReplaceString(out, "\n ! ", "\n "); // hack to handle bug in Magma string formatting
out := ReplaceString(out, "EquationOrder(Polynomial(\\[1, 1, 1]))", "Integers(CyclotomicField(3))");
t2 := Cputime();
t2 - t1;

Write(output_file,out);
"Coleman data recorded.";


