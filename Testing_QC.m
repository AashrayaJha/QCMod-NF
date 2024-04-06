AttachSpec("QCMod.spec");
// load "data/NF-example-coleman-data.m";
// load "data/New_hecke.m";

// data_1:=data_1;
// data_2:=data_2;

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
pts := QCModAffine(Q,p: data1:=data_1,data2:=data_2, known_points:=known_affine_points, correspondence_data := correspondence_data);