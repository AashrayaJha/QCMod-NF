//In this file, we check the values of the j-invarainats corresponding to the F-points of XH using maps between 
//modular curves provided by Jeremy Rouse.

K<zeta> := CyclotomicField(3);

RR<A,B,C> := PolynomialRing(K,3);
pol := A^4 + (zeta - 1)*A^3*B + (3*zeta + 2)*A^3*C - 3*A^2*C^2 +
    (2*zeta + 2)*A*B^3 - 3*zeta*A*B^2*C + 3*zeta*A*B*C^2 - 2*zeta*A*C^3 -
    zeta*B^3*C + 3*zeta*B^2*C^2 + (-zeta + 1)*B*C^3 + (zeta + 1)*C^4;
XH := Curve(ProjectiveSpace(K,2),pol);

points:= [XH![P[2],P[1],P[3]] : P in [
  [1,0,0], // j = 1728, D = -4
  [1,zeta+1,0], // j = 287496, D = -16 
  [0,-zeta-1,1], // j = 1728, D = -4
  [1,-zeta-1,1], // j = -884736000, D = -43
  [zeta+1,-zeta-1,1], // j = 1728, D = -4
  [0,-zeta,1], // j = -3375, D = -7
  [zeta+1,0,1], // j = -884736, D = -19
  [2*zeta+2,zeta,1], // j = 16581375, D = -28
  [zeta,1,1], // j = 1728, D = -4
    [ (1/3)*(-zeta-2), (1/3)*(zeta+2),1], // j = -147197952000, D = -67
  [(-1/2)*zeta,-1/2,1], // j = 1728, D = -4
  [(1/7)*(5*zeta+4),-1,1],// j=-262737412640768000, D = -163
  [(1/2)*(zeta-3),(1/2)*(zeta+2),1] // j = (-3238903991430*zeta - 4786881515880)^3/19^27, non-CM
]];

pt:=points[13];
P1 := ProjectiveSpace(K,1);

phi2 := map<XH -> P1 | [[(zeta-1)*XH.1+(-2*zeta-1)*XH.2+(zeta-1)*XH.3,(-zeta-2)*XH.1+(-zeta-2)*XH.2+(-zeta+1)*XH.3],[1/3*(-zeta + 1)*XH.1^3 + 1/3*(4*zeta - 1)*XH.1^2*XH.2 + 1/3*(-4*zeta + 1)*XH.1*XH.2^2
    + (zeta + 1)*XH.2^3 + 1/3*(2*zeta + 4)*XH.1^2*XH.3 + zeta*XH.1*XH.2*XH.3 + 
    1/3*(-8*zeta - 4)*XH.2^2*XH.3 + 1/3*(zeta - 1)*XH.1*XH.3^2 + 1/3*(5*zeta + 
    1)*XH.2*XH.3^2 + 1/3*(-2*zeta - 1)*XH.3^3,XH.2^3 + (zeta - 1)*XH.2^2*XH.3 + (-2*zeta - 1)*XH.2*XH.3^2]]>;

// Map from X_ns^+(9) to the j-line
phi3 := map<P1 -> P1 | [(-15*P1.1^9 - 81*P1.1^8*P1.2 - 27*P1.1^7*P1.2^2 + 117*P1.1^6*P1.2^3 
	- 81*P1.1^5*P1.2^4 + 189*P1.1^4*P1.2^5 - 288*P1.1^3*P1.2^6 + 216*P1.1*P1.2^8 - 96*P1.2^9)^3,
(P1.1^9 - 9*P1.1^7*P1.2^2 + 3*P1.1^6*P1.2^3 + 27*P1.1^5*P1.2^4 - 18*P1.1^4*P1.2^5 - 
	24*P1.1^3*P1.2^6 + 27*P1.1^2*P1.2^7 - 9*P1.1*P1.2^8 + P1.2^9)^3]>;

zeta_3:=zeta;

j1:=2^3 * 5^3 * 19^(-27) * (-15- 2*zeta)^3 * (452 +431*zeta)^3 * (53515427 + 41314914*zeta)^3; //RSZB
j2:=2^3*5^3*19^(-27)*(1- zeta_3)^6*(2-3*zeta_3)^(27)*(27+13*zeta_3)^3*(54+49*zeta_3)^3*(227+173*zeta_3)^3; //overleaf

phi:=phi2*phi3; //j map from XH

pt1:= P1![j1,1]; //RSZB
pt2:= P1![j2,1]; //overleaf

if phi(pt)[1] eq j1 then
  "The RSZB j-invariant is correct.";
elif phi(pt)[1] eq j2 then
  "The j-invariant in our paper is correct.";
end if;

if #RationalPoints(pt1@@phi3) eq 0 then 
  "The j-invariant in RSZB does not lift to Xns+(9).";
end if;

j_inv:=[phi(pts): pts in points]; //Checked other j-invaraints are correct.

