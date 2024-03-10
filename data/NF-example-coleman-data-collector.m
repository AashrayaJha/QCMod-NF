AttachSpec("~/GitHub/CHIMP/CHIMP.spec");
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
Q := y^4 + (u - 1)*y^3*x + (2*u + 2)*y*x^3 + (3*u + 2)*y^3 - 3*u*y*x^2 - u*x^3 - 3*y^2 + 3*u*y*x + 3*u*x^2 - 2*u*y + (-u + 1)*x + (u + 1);

p := 13;
v := Factorization(p*OK)[2][1];
N := 160 ; //Going to try 40

"Constructing symplectic basis of H1...";
h1basis, g, r, W0 := H1Basis(Q, v);
prec := 2*g;
cpm := CupProductMatrix(h1basis, Q, g, r, W0 : prec := prec, split := true);
sympl := SymplecticBasisH1(cpm);
new_complementary_basis := [&+[sympl[i,j]*h1basis[j] : j in [1..2*g]] : i in [1..g]];
sympl_basis := [h1basis[i] : i in [1..g]] cat new_complementary_basis;
basis0 := [[sympl_basis[i,j] : j in [1..Degree(Q)]] : i in [1..g]]; // basis of regular differentials
basis1 := [[sympl_basis[i,j] : j in [1..Degree(Q)]] : i in [g+1..2*g]];  // basis of complementary subspace
"Symplectic basis constructed.";

"Constructing Coleman data...";
time data := ColemanData(Q, v, N : useU:=true,  basis0:=basis0, basis1:=basis1);
"Coleman data constructed.";

"Recording Coleman data...";
output_file := "data/NF-example-coleman-data-13_160.m";
fprintf output_file, "K<u> := CyclotomicField(3);\n";
fprintf output_file, "_<x> := PolynomialRing(K);\n";
fprintf output_file, "_<z> := LaurentSeriesRing(PolynomialRing(K));\n\n";

out := Sprintf("data := %m;\n\n", data);
out := ReplaceString(out, "\n ! ", "\n "); // hack to handle bug in Magma string formatting
out := ReplaceString(out, "EquationOrder(Polynomial(\\[1, 1, 1]))", "Integers(CyclotomicField(3))");

Write(output_file,out);
"Coleman data recorded.";
