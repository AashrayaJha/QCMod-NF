Q3:=pAdicField(3,100);
R<y>:=PolynomialRing(Q3);
K:=TotallyRamifiedExtension(Q3,y^2-3);
S<x>:=PolynomialRing(K);
zeta_3:=(-1+K.1)/2;
psi:= x^9 + (9*zeta_3 - 9)*x^8 + (54*zeta_3 + 27)*x^7+ (54*zeta_3 - 27/2)*x^6\
+ (243*zeta_3 + 972)*x^5+ 729*zeta_3*x^4+ (2916*zeta_3 - 1458)*x^3\
+ (37179*zeta_3 + 41553)*x^2+ (6561*zeta_3 + 6561/8)*x - 63423*zeta_3 + 155277; 
IsIrreducible(psi);
G,S:=GaloisGroup(psi);//Don't think we will need this but good to have
