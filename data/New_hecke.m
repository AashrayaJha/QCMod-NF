AttachSpec9

Z27:=Integers(27);
G:= GeneralLinearGroup(2,Z27);
H1:=Matrix(Z27,2,2,[0,26,4,6]);
H2:=Matrix(Z27,2,2,[10,1,25,26]);
H:=sub<G|H1,H2>;
NH:=ncl<G|H>;
//Index(G,NH); //This is 6
//Index (NH,H); //This is 27

//Want to find alpha in M2(Z) reducing to H with determinant p, and sigma_p in Gamma_(NH) such that sigma_p=p.alpha mod 27

function findalpha(p):

w_p:=Matrix(2,2,Integers,[p,0,0,1]);
GammaN:=CongruenceSubgroup(27);  //Since we are working with the level 27 curve
gamma_1:=RandomElement(GammaN);
gamma_2:=RandomElement(GammaN);

return gamma1*w_p*gamma_2


