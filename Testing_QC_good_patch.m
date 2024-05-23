AttachSpec("QCMod.spec");
load "data/NF-example-coleman-data-good_patch_2.m";
load "data/New_hecke_patch2_401.m";

data_1:=data_1;
data_2:=data_2;

Q:=data_1`Q;
p:=data_1`p;
v_1:=data_1`v; //We need both primes above p. Here v1,v2 will be those primes.
v_2:=data_2`v; //Will calculate them in the intrinsic.
//print "v_1: ", v_1;
//print "v_2: ", v_2;

correspondence_data:= [*AK_patch_2_401,Zs_patch_2,-3*];
// The -6 correspinds to prec_loss for correspondences. This was set as the val(det(Finv),p) over Q, so
// did the analogous thing over K.

//new polynomial: y=a,x=c,b=1.
//Q:= y^4 + ((-2*u + 9)*x + (2*u + 3))*y^3 + (-3*x^2 + 6*x - 3)*y^2 + ((-170*u + 254)*x^3 + (-150*u + 114)*x^2 + (-54*u + 18)*x - 10*u - 2)*y + (162*u + 144)*x^4 + (-108*u + 48)*x^3 + (-72*u - 144)*x^2 + (12*u - 48)*x + 6*u;
known_affine_points:=[*[1,-8*u - 4],[1,0]*];
//initialising points in new coordinates
print "starting QCModAffine";
SetVerbose("QCMod",2);
sols, all_zeroes, double_zeroes, global_pts_local, F1_lists, F2_lists, Qppoints_1, Qppoints_2 := QCModAffine(Q,p: data1:=data_1,data2:=data_2,known_points:=known_affine_points, correspondence_data := correspondence_data, N := 15);

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
    assert minval1 ge 4 and minval2 ge 4 ;
      //print "i,j", i,j,minval1,minval2;
    //end if;
  end for;
end for;
"All checks passed: images of known affine points among solutions";
"Number of known affine points", #Qpts;
"Number of solutions", #sols;
"Any multiple roots?";
&or[s[2] : s in sols];



// SetVerbose("QCMod",3);
// all_zeroes, double_zeroes, E1_lists_1, E1_lists_2 := QCModAffine(Q,p: data1:=data_1,data2:=data_2, known_points:=known_affine_points, correspondence_data := correspondence_data, N := 13);

