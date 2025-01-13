//We verify Remark A.4 which says that the two definition of modular curves, one due to Shimura and the other in 
//RSZB give isomorphic curves for the subgroup H considered in our paper.

R := Integers(27);
G:= GeneralLinearGroup(2,R);
A := G![0, 26, 4, 6];
B := G![10, 1, 25, 26];
H := MatrixGroup< 2, R | A,B >; //This is H defined in the paper.

P := G![0,1,-1,0];
P*Transpose(A)*P^(-1) in H;
//returns true
P*Transpose(B)*P^(-1) in H;
//returns true