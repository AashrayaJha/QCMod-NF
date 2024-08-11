AttachSpec("QCMod.spec");
load "data/Coleman_data_other_patches/Coleman_data_patch_3.m";
load "data/correspondences_all.m";

import "applications.m": Qp_points;
import "singleintegrals.m": is_bad, xy_coordinates;


t1:=Cputime();

Q:=data_1`Q;

Qppoints_1 := Qp_points(data_1);
assert &and[not(is_bad(P, data_1)) : P in Qppoints_1];
Qppoints_2 := Qp_points(data_2);
assert &and[not(is_bad(P, data_2)) : P in Qppoints_2];
// So under both embeddings, all residue disks are good (in particular,
// affine)


Q_RSZB := Polynomial([PolynomialRing(CyclotomicField(3)) | [[ RationalField() | 1, 1 ], [ RationalField() | 1, -1 ], [ RationalField() | 0, 3 ], [ RationalField() | 0, -1 ]], [[ RationalField() | 0, -2 ], [ RationalField() | 0, 3 ], [ RationalField() | 0, -3 ], [ RationalField() | 2, 2 ]], [[ RationalField() | -3, 0 ]], [[ RationalField() | 2, 3 ], [ RationalField() | -1, 1 ]], [[ RationalField() | 1, 0 ]]]);

C := CurveFromBivariate(Q);
C_RSZB := CurveFromBivariate(Q_RSZB);
bool, trans_seq := IsIsomorphicPlaneQuartics(C_RSZB,C);
trans := trans_seq[1];
a,b,c := GetVersion();
if b lt 28 then
  trans:=trans^(-1);  //AJ: Need this for Magma version V2.27-5
                      // TODO: Find out when this change was made.
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
  [(1/2)*(u-3),(1/2)*(u+2),1], // j = (-3238903991430*u - 4786881515880)^3/1/9^27, non-CM
  [(1/3)*(-u-2), (1/3)*(u+2),1], // j = -147197952000, D = -67
  [(-1/2)*u,-1/2,1], // j = 1728, D = -4
  [(1/7)*(5*u+4),-1,1]// j=-262737412640768000, D = -163
]];
known_points_C := [C![&+[trans[i,j]*P[j] : j in [1..3]] : i in [1..3]] 
                                    : P in known_points_RSZB];


p:=data_1`p;
v_1:=data_1`v; //We need both primes above p. Here v1,v2 will be those primes.
v_2:=data_2`v; //Will calculate them in the intrinsic.
//print "v_1: ", v_1;
//print "v_2: ", v_2;

correspondence_data:= correspondence_data_list[3]`correspondence_data;
// The -6 corresponds to prec_loss for correspondences. This was set as the val(det(Finv),p) over Q, so
// did the analogous thing over K.

//new polynomial: y=a,x=c,b=1.
//Q:= y^4 + ((-2*u + 9)*x + (2*u + 3))*y^3 + (-3*x^2 + 6*x - 3)*y^2 + ((-170*u + 254)*x^3 + (-150*u + 114)*x^2 + (-54*u + 18)*x - 10*u - 2)*y + (162*u + 144)*x^4 + (-108*u + 48)*x^3 + (-72*u - 144)*x^2 + (12*u - 48)*x + 6*u;

assert &and[P[3] ne 0 : P in known_points_C];
known_affine_points := [[P[1], P[2]] : P in known_points_C];
//initialising points in new coordinates
print "starting QCModAffine";
SetVerbose("QCMod",2);
sols, all_zeroes, double_zeroes, global_pts_local, F1_lists, F2_lists, Qppoints_1, Qppoints_2 := QCModAffine(Q,p: data1:=data_1,data2:=data_2,known_points:=known_affine_points, correspondence_data := correspondence_data, N := 15);

// Qpts contains the images of the known points under the 2 embeddings.
Qpts := [];
for i := 1 to #global_pts_local do
  pt1 := xy_coordinates(global_pts_local[i,1], data_1);
  pt2 := xy_coordinates(global_pts_local[i,2], data_2);
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
"Check that there are no multiple roots.";
assert not &or[s[2] : s in sols];
"No multiple roots!";
printf "Hence the set of K-rational points on the %o is \n%o,\n where K is the %o.",  C_RSZB, known_points_RSZB, BaseRing(C);


t2:=Cputime();
printf "This is the time taken: %o seconds. \n", t2-t1;


// SetVerbose("QCMod",3);
// all_zeroes, double_zeroes, E1_lists_1, E1_lists_2 := QCModAffine(Q,p: data1:=data_1,data2:=data_2, known_points:=known_affine_points, correspondence_data := correspondence_data, N := 13);

