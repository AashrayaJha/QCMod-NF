AttachSpec("QCMod.spec");
// import "auxpolys.m": auxpolys, log;
// import "singleintegrals.m": evalf0, is_bad, local_coord, set_point, tadicprec, teichmueller_pt, xy_coordinates;
// import "misc.m": are_congruent, equivariant_splitting, eval_mat_R, eval_list, eval_Q, FindQpointQp, function_field, alg_approx_Qp, minprec, minval, minvalp, QpMatrix,QpSequence,QpPolynomial;
// import "applications.m": Q_points, Qp_points, roots_with_prec, separate;
// import "heights.m": E1_tensor_E2_NF, expand_algebraic_function, frob_equiv_iso, height;
// load "data/NF-example-coleman-data.m";
// load "data/New_hecke.m";

// data1:=data_1;
// data2:=data_2;

Q:=data_1`Q;
p:=data_1`p;
v1:=data_1`v; //We need both primes above p. Here v1,v2 will be those primes.
v2:=data_2`v; //Will calculate them in the intrinsic.
correspondence_data:= [*AK_160,Zs,-3*];
N:=15; corr_loss:=-3; prec:=15; global_base_point_index := 1; hodge_prec := 5; 
Kv,psi:=Completion(K,v1);

search_bound := 1000;
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
known_points:=[Prune(known_projective_points[i]): i in [3..13]];
Qpoints_1 := Q_points(data1,search_bound : known_points := known_points);
Qppoints_1 := Qp_points(data1 : Nfactor :=15);
good_Qpoints_1 := [P : P in Qpoints_1 | not is_bad(P, data1) and not P`inf];
i1:=9;  j1:=17; Qpt1:=good_Qpoints_1[i1]; //Indices of relevant good point (good_Qpoints_1,Qppoints_1)
i2:=8;  j2:=16; Qpt2:=good_Qpoints_1[i2]; // Indices of relevant bad point

good_Q_Qp_indices_1 := [FindQpointQp(P,Qppoints_1) : P in good_Qpoints_1];
good_coordinates_1 := [xy_coordinates(P,data1) : P in good_Qpoints_1];
good_affine_rat_pts_xy := [[alg_approx_Qp(P[1], v1), alg_approx_Qp(P[2], v1)] : P in good_coordinates_1]; 

teichpoints_1 := [**]; teichpoints_2 := [**];
for i in [1..20] do
    teichpoints_1[i] := is_bad(Qppoints_1[i],data1) select 0  else teichmueller_pt(Qppoints_1[i],data1); // No precision loss
end for;

bQ_1 := good_Qpoints_1[1]; b01 := teichmueller_pt(bQ_1,data1); 
b0pt1 := [Rationals()!c : c in xy_coordinates(b01, data1)]; //base points are the same
local_base_point_index_1 := FindQpointQp(bQ_1,Qppoints_1);
ks_1:=Exclude(good_Q_Qp_indices_1, local_base_point_index_1);

Z:=Zs[1];Tq:=AK_160;

h1basis, g, r, W0 := H1Basis(Q, v1); 
bQ_xy := good_affine_rat_pts_xy[1];
FF<y>  := function_field(Q);x := BaseRing(FF).1;
bpt   := CommonZeros([x-bQ_xy[1], y-bQ_xy[2]])[1];

Ncorr := N + Min(corr_loss, 0);

eta1,betafil1,gammafil1,hodge_loss1 := HodgeData(Q,g,W0,data1`basis,Z,bpt : r:=r, prec:=hodge_prec);
Nhodge := Ncorr + Min(0, hodge_loss1);
Z1 := QpMatrix(Z, Nhodge, v1);
Z1 := ChangeRing(Z1, Rationals());
PhiAZb_to_b01, Nptb01 := ParallelTransport(bQ_1,b01,Z1,eta1,data1:prec:=prec,N:=Nhodge);
for i := 1 to 2*g do
      PhiAZb_to_b01[2*g+2,i+1] := -PhiAZb_to_b01[2*g+2,i+1]; //this seems correct.
end for;

try 
    pt1, Npt1 := ParallelTransport(Qpt1,teichpoints_1[j1], Z1,eta1,data1:prec:=prec,N:=Nhodge); //This seems negative 
catch e
    printf "There is some error with indexing";
end try;   

G1, NG1 := FrobeniusStructure(data1,Z1,eta1,b0pt1 : N:=Nhodge); 

pt := [IntegerRing()!c : c in xy_coordinates(Qpt1, data1)]; 
ptbad:= [IntegerRing()!c : c in xy_coordinates(Qpt1, data1)];//This seems correct.
G_list1 := eval_mat_R(G1, pt, r, v1); //This does not seem correct. THis is the G_list1[17] in QCMod
Ncurrent := Min(Nhodge, NG1);
isoi1, Nisoi1 := frob_equiv_iso(G_list1,data1,Ncurrent);  //This seems to be negative of what it should be.
MNi1 := Npt1 lt Nisoi1 select Parent(pt1) else Parent(isoi1);
PhiAZb1 := MNi1!(pt1*PhiAZb_to_b01*isoi1); 
Phii1 := MNi1!(pt1*PhiAZb1);
Ni1 := Min([Ncurrent, Precision(BaseRing(Phii1)), minprec(Phii1)]);
x1, y1 := Explode(xy_coordinates(Qpt1, data1));
x2, y2 := Explode(xy_coordinates(Qpt2, data1));
gammafilP_1 := eval_list(Eltseq(gammafil1), x1, y1, v1, Ni1);
gammafilP_2 := eval_list(Eltseq(gammafil1), x2, y2, v1, Ni1);
eqsplit := equivariant_splitting(Tq);
eqsplit1 := QpMatrix(eqsplit,Ncorr,v1);
heightP1 := height(Phii1,QpSequence(Eltseq(betafil1),Ni1,v1),gammafilP_1,eqsplit1,data1);