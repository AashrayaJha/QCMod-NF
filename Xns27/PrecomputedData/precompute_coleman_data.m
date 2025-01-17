//This file collects Coleman data, for one prime of the number field at a time. 
//The main intrinsic is ColemanData (in singleintegral.m), which in particular computes 
//Frobenius action on H^1_dR(X), for which we use a symplectic basis.
//
//We use the primes above 13.

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

Q := Polynomial([PolynomialRing(CyclotomicField(3)) | [[ RationalField() | 0, 6 ], [ RationalField() | -106, -74 ], [ RationalField() | 198, -318 ], [ RationalField() | 642, 970 ], [ RationalField() | -734, -584 ]], [[ RationalField() | -2, -10 ], [ RationalField() | 120, 30 ], [ RationalField() | 330, 606 ], [ RationalField() | -720, 542 ]], [[ RationalField() | -3, 0 ], [ RationalField() | 6, 0 ], [ RationalField() | -3, 0 ]], [[ RationalField() | 3, 2 ], [ RationalField() | 13, 12 ]], [[ RationalField() | 1, 0 ]]]);

p := 13;
v := Factorization(p*OK)[2][1];
N := 20; 

"Constructing symplectic basis of H1...";
h1basis, g, r, W0 := H1Basis(Q, v);
prec := 2*g;
cpm := CupProductMatrix(h1basis, Q, g, r, W0 : prec := prec, split := false);
sympl := SymplecticBasisH1(cpm);
new_complementary_basis := [&+[sympl[i,j]*h1basis[j] : j in [1..2*g]] : i in [1..g]];
sympl_basis := [h1basis[i] : i in [1..g]] cat new_complementary_basis;
basis0 := [[sympl_basis[i,j] : j in [1..Degree(Q)]] : i in [1..g]]; // basis of regular differentials
basis1 := [[sympl_basis[i,j] : j in [1..Degree(Q)]] : i in [g+1..2*g]];  // basis of complementary subspace
"Symplectic basis constructed.";

"Constructing Coleman data...";
time data := ColemanData(Q, v, N : useY:=true,  basis0:=basis0, basis1:=basis1);
"Coleman data constructed.";

"Recording Coleman data...";
output_file := "data/Coleman_data_other_patches/Coleman_data_patch_3.m";
fprintf output_file, "K<u> := CyclotomicField(3);\n";
fprintf output_file, "_<x> := PolynomialRing(K);\n";
fprintf output_file, "_<z> := LaurentSeriesRing(PolynomialRing(K));\n\n";

out := Sprintf("data := %m;\n\n", data);
out := ReplaceString(out, "\n ! ", "\n "); // hack to handle bug in Magma string formatting
out := ReplaceString(out, "EquationOrder(Polynomial(\\[1, 1, 1]))", "Integers(CyclotomicField(3))");

Write(output_file,out);
"Coleman data recorded.";
