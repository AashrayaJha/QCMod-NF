AttachSpec("QCMod.spec");
load "data/NF-example-coleman-data.m";
load "data/New_hecke.m";
SetLogFile("debugging_20.log");

import "singleintegrals.m": evalf0, is_bad, local_coord, set_point, tadicprec, teichmueller_pt, xy_coordinates;
import "misc.m": are_congruent, equivariant_splitting, eval_mat_R, eval_list, eval_Q, FindQpointQp, function_field, alg_approx_Qp, minprec, minval, minvalp, QpMatrix,QpSequence,QpPolynomial;


//data_1:=data_1;
//data_2:=data_2;

Q:=data_1`Q;
p:=data_1`p;
v_1:=data_1`v; //We need both primes above p. Here v1,v2 will be those primes.
v_2:=data_2`v; //Will calculate them in the intrinsic.
//print "v_1: ", v_1;
//print "v_2: ", v_2;

correspondence_data:= [*AK_160,Zs,-3*];
// The -6 correspinds to prec_loss for correspondences. This was set as the val(det(Finv),p) over Q, so
// did the analogous thing over K.
known_projective_points := [
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
  ];

known_affine_points:=[Prune(known_projective_points[i]): i in [3..13]];
//print "starting QCModAffine";

SetVerbose("QCMod",3);
sols, all_zeroes, double_zeroes, global_pts_local, F1_lists, F2_lists, Qppoints_1, Qppoints_2 := QCModAffine(Q,p: data1:=data_1,data2:=data_2, known_points:=known_affine_points, correspondence_data := correspondence_data, N := 15);

// Qpts contains the images of the known points under the 2 embeddings.
Qpts := [];
for i := 1 to 11 do
  pt1 := xy_coordinates(global_pts_local[i,1], data_1);
  pt2 := xy_coordinates(global_pts_local[i,2], data_2);
  Append(~Qpts, [pt1, pt2]);
end for;

"Check if known points are recovered";
for i in [1..#Qpts] do
  for j in [1..#sols] do
    minval1 := Min(Valuation(Qpts[i,1,1]-sols[j,1,1,1]), Valuation(Qpts[i,1,2]-sols[j,1,1,2]));
    minval2 := Min(Valuation(Qpts[i,2,1]-sols[j,1,2,1]), Valuation(Qpts[i,2,2]-sols[j,1,2,2]));
    assert minval1 ge 4 and minval2 ge 4 then
      //print "i,j", i,j,minval1,minval2;
    end if;
  end for;
end for;

"Any multiple roots?";
&or[s[2] : s in sols];

/*
"Second correspondence";

for i in [1..#Qpts] do
  for j in [1..#sols[2,i]] do
    if are_congruent(Qpts[i,1], sols[2,j,1]) and 
            are_congruent(Qpts[i,2], sols[2,j,2]) then
      "i,j", i,j;
    end if;
  end for;
end for;

common_zeroes := [**];
for i in [1..#all_zeroes[1]] do 
  common_zeroes[i] := [**];
  for j in [1..#all_zeroes[1]] do
    common_zeroes[i][j] := [**];
    if #all_zeroes[1,i]*#all_zeroes[2,i] gt 0 then
      if IsDefined(all_zeroes[1,i],j) and IsDefined(all_zeroes[2,i],j) then
        for pair1 in all_zeroes[1,i,j] do
          for pair2 in all_zeroes[2,i,j] do
            if Min(Valuation(pair1[1]-pair2[1]), Valuation(pair1[2]-pair2[2])) ge 2 then 
         //   if are_congruent(pair1, pair2) then
              Append(~common_zeroes[i,j], pair1);
            end if;
          end for;
        end for;
      end if;
    end if;
  end for;
end for;
      
*/
