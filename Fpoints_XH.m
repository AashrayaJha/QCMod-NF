// load this file in the main directory.

AttachSpec("QCMod.spec");
import "applications.m": Qp_points;
import "singleintegrals.m": is_bad, xy_coordinates;

// load precomputed Coleman data -- this saves some time. Moreover, the bases 
// of H1 are already symplectic.
load "Xns27/PrecomputedData/coleman_data_precomputed.m";

// load precomputed Hecke data (Hecke operator T_13 and associated nice
// correspondences, encoded as matrices wrt symplectic basis of H1). We do 
// this because a lot of precision was
// necessary to recognize the entries as algebraic numbers, so computing
// them takes very long.
load "Xns27/PrecomputedData/hecke_data_precomputed.m";


t1 := Cputime();
p := data_1`p;
Q := data_1`Q;
y := Parent(Q).1;
x := BaseRing(Q).1;
assert Q eq y^4 + ((-2*u + 9)*x + (2*u + 3))*y^3 + (-3*x^2 + 6*x - 3)*y^2 + ((-170*u + 254)*x^3 + (-150*u + 114)*x^2 + (-54*u + 18)*x - 10*u - 2)*y + (162*u + 144)*x^4 + (-108*u + 48)*x^3 + (-72*u - 144)*x^2 + (12*u - 48)*x + 6*u;


// Check that condition (2) is satisfied. Condition (1) also holds, this is
// checked in ColemanData.
Qppoints_1 := Qp_points(data_1);
assert &and[not(is_bad(P, data_1)) : P in Qppoints_1];
Qppoints_2 := Qp_points(data_2);
assert &and[not(is_bad(P, data_2)) : P in Qppoints_2];

/* 
 * Alternatively, can use the following, so that the Coleman data gets
 * computed inside the main function OCModAffine.
 *
F<u> := CyclotomicField(3);
Fx<x> := PolynomialRing(F);
Fxy<y> := PolynomialRing(Fx);
p := 13;
N := 15;
Q := y^4 + ((-2*u + 9)*x + (2*u + 3))*y^3 + (-3*x^2 + 6*x - 3)*y^2 + ((-170*u + 254)*x^3 + (-150*u + 114)*x^2 + (-54*u + 18)*x - 10*u - 2)*y + (162*u + 144)*x^4 + (-108*u + 48)*x^3 + (-72*u - 144)*x^2 + (12*u - 48)*x + 6*u;
*/



// Curve as in RSZB
Q_RSZB := Polynomial([PolynomialRing(CyclotomicField(3)) | [[ RationalField() | 1, 1 ], [ RationalField() | 1, -1 ], [ RationalField() | 0, 3 ], [ RationalField() | 0, -1 ]], [[ RationalField() | 0, -2 ], [ RationalField() | 0, 3 ], [ RationalField() | 0, -3 ], [ RationalField() | 2, 2 ]], [[ RationalField() | -3, 0 ]], [[ RationalField() | 2, 3 ], [ RationalField() | -1, 1 ]], [[ RationalField() | 1, 0 ]]]);

C := CurveFromBivariate(Q);
C_RSZB := CurveFromBivariate(Q_RSZB);
bool, trans_seq := IsIsomorphicPlaneQuartics(C_RSZB,C);
trans := trans_seq[1];
a,b,c := GetVersion();
if b lt 28 then
  if c lt 10 then
    trans:=trans^(-1);  
  end if;
end if;
known_points_RSZB := [C_RSZB!P : P in [
  [1,0,0], // j = 1728, D = -4
  [1,u+1,0], // j = 287496, D = -16 
  [0,-u-1,1], // j = 1728, D = -4
  [1,-u-1,1], // j = -884736000, D = -43
  [u+1,-u-1,1], // j = 1728, D = -4
  [0,-u,1], // j = -3375, D = -7
  [u+1,0,1], // j = -884736, D = -19
  [2*u+2,u,1], // j = 16581375, D = -28
  [u,1,1], // j = 1728, D = -4
  [(1/2)*(u-3),(1/2)*(u+2),1], // j = (-3238903991430*u - 4786881515880)^3/19^27, non-CM
  [(1/3)*(-u-2), (1/3)*(u+2),1], // j = -147197952000, D = -67
  [(-1/2)*u,-1/2,1], // j = 1728, D = -4
  [(1/7)*(5*u+4),-1,1]// j=-262737412640768000, D = -163
]];
known_points_C := [C![&+[trans[i,j]*P[j] : j in [1..3]] : i in [1..3]] 
                                    : P in known_points_RSZB];

assert &and[P[3] ne 0 : P in known_points_C]; // no F-pts at infinity


correspondence_data:= [*AK_good_patch,Zs_good_patch,prec_loss_bd*];

known_affine_points := [[P[1], P[2]] : P in known_points_C];
//initialising points in new coordinates
print "starting QCModAffine";
SetVerbose("QCMod",3);

recovered_Kpts, done, sols, all_zeroes, double_zeroes, global_pts_local, bad_affine_K_pts, data1, data2 := QCModAffine(Q,p, known_affine_points, correspondence_data: data1:=data_1,data2:=data_2, N := 15, symplectic);

/* 
 * Alternatively, use the following to run QC without precomputed Coleman
 * data.
 *
recovered_Kpts, done, sols, all_zeroes, double_zeroes, global_pts_local, bad_affine_K_pts, data1, data2  := QCModAffine(Q,p, known_affine_points, correspondence_data: N := N);
 */

// Qpts contains the images of the known points under the 2 embeddings.
Qpts := [];
for i := 1 to #global_pts_local do
  pt1 := xy_coordinates(global_pts_local[i,1], data1);
  pt2 := xy_coordinates(global_pts_local[i,2], data2);
  Append(~Qpts, [pt1, pt2]);
end for;

"Check if known points are recovered";
for i in [1..#Qpts] do
  maxminval := 0;
  for j in [1..#sols] do
    minval1 := Min(Valuation(Qpts[i,1,1]-sols[j,1,1,1]), Valuation(Qpts[i,1,2]-sols[j,1,1,2]));
    minval2 := Min(Valuation(Qpts[i,2,1]-sols[j,1,2,1]), Valuation(Qpts[i,2,2]-sols[j,1,2,2]));
    //  print "i,j", i,j,minval1,minval2;
    //end if;
    maxminval := Max([maxminval, minval1, minval2]);
  end for;
  assert maxminval ge 4;
end for;
"All checks passed: images of known affine points among solutions";
"Number of known affine points", #Qpts;
"Number of solutions", #sols;
"Check that there are no roots that are multiple for both systems.";
assert not &or[s[2] : s in sols];
"No multiple roots!";
printf "Hence the set of K-rational points on the %o is \n%o,\n where K is the %o.\n",  C_RSZB, known_points_RSZB, BaseRing(C);

t2:=Cputime();
printf "This is the time taken: %o\n", t2-t1;



