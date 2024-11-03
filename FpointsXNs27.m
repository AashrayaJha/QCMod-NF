R<a,b,c,d,e,f,g,h,i,j,k,l> := PolynomialRing(Rationals(),12);
// Canonical model of X_ns^+(27) from David Zywina's code.
psi := [
    a*b - a*d - a*k - b*e - b*f + b*h - 
        b*j - b*k - c^2 - c*d - c*k + j^2,
    a*c - a*d - 2*a*f + a*h - a*i - a*j - 
        a*k + 2*b*c - b*g + b*h - c*e - c*f
        + c*j + c*l + g*k - h*k - i*j - 
        j*k,
    a*b - 2*a*d - 2*a*f + a*h - b^2 - b*c - 
        b*f - b*g + b*i + b*j - c*g + c*h -
        c*i + g*k + h^2 - h*k + h*l - i^2 - 
        i*l + j^2 - j*k - k^2,
    -2*a^2 - a*b - a*c + 2*a*e + a*f - a*g - 
        a*h + a*i + a*j - b^2 - 2*b*c + b*e - 
        b*h + b*i + 2*b*j - 2*b*k + 2*c^2 - 
        c*d - c*e - c*g - c*i - c*j + d*g +
        j*l,
    -a^2 - a*b - a*f - a*h + a*i + a*j - 
        2*a*k + a*l - 3*b^2 - b*c + b*f - 
        b*k - c*e - c*g + c*h - 2*c*i - 
        c*j + d*g + g*k - h*k - i*k + j^2 
        - j*k + j*l - k^2 + k*l,
    a*b + a*e - 2*a*g + 2*a*h - a*j - a*k 
        - a*l - b^2 + 2*b*c + b*e + b*g - b*h 
        - 3*b*k + c^2 - c*f - c*g - c*h + c*i 
        + c*j + d*g - h^2 + i*l + k^2 - k*l,
    -2*a*h + a*i - a*k - 2*b^2 - b*d + c^2 - 
        2*c*e + c*f + c*g - c*h + c*i + c*j
        + d*g - h^2 - h*i - h*l + i^2 + i*k + 
        j*k + k*l,
    -a*b - a*d - 2*a*e - 2*a*f - a*i - a*j 
        - a*k - 2*b^2 + 2*b*d - b*e + b*f - 
        b*g - b*h - b*i + b*k + c*d + 2*c*f
        + h*j + i*k - i*l + j*k - k^2 + 
        k*l,
    -a^2 - a*b + a*e - a*g - a*h - a*k - 
        a*l + b*c - b*e + b*f - b*g + b*h -
        b*i - b*k + c*d - c*g + c*h + c*k 
        + c*l + d*g + g*k - h^2 - h*j - i*k
        + i*l - j^2 + k^2 - k*l + l^2,
    -a^2 + 2*a*b - a*e - a*g + a*h + a*i + 
        2*a*j - a*l + b^2 + b*c + b*d - b*e -
        b*h + b*i + c*d + 2*c*e - 2*c*g - h^2 +
        h*i + h*k + i*j - i*k + i*l + k^2 
        - 2*k*l,
    -a^2 - a*c + a*e - a*h + a*j - a*k + 
        a*l + b*d + 2*b*f - b*g + b*h - b*i
        + b*k + 2*c^2 + c*d - 2*c*e - c*h + 
        c*i + c*j + c*k + g*k - h^2 - 2*h*l
        + i^2 + i*l - j^2 + k^2 - k*l,
    2*a^2 + 2*a*b + 2*a*h - a*i + b^2 - b*c + 
        2*b*e + b*f + b*g - b*h + 2*b*i - c^2 +
        2*c*e - c*f - c*g + c*k - c*l - 
        g*k + g*l + h^2 + h*l - i^2 - j*k - 
        k*l,
    -a*b + a*e + a*f + a*i - a*j + a*l - 
        b^2 - b*f + b*h - b*i - b*k + c^2 - 
        c*d - 2*c*e - c*f + 2*c*g - c*h - 
        c*k + g*i + g*j - h*i - h*j - 
        h*l - i*j + k*l,
    -2*a^2 - a*b - a*d - a*f + a*g - a*h - 
        a*j - a*k - b*c + b*d - 2*b*e - 
        2*b*g + b*j + b*k + c^2 - c*d + c*k 
        + c*l + g*j + g*k - h*k - i*j - 
        j^2 - k*l + l^2,
    -a^2 - 2*a*e - a*f + a*g - a*h + a*j + 
        a*k - b*e - b*f - b*h + 2*b*k + 
        2*c*f - c*k - d*g + g*i - g*k + 
        g*l + 2*h*k - h*l + i*k - i*l + 
        j*k,
    -a^2 - a*b + a*c + a*d + a*e + a*f + 
        a*g + a*j + a*k + a*l + b^2 + b*e + 
        b*h + b*i + b*j + b*k - c*d - c*e 
        - c*f - c*i - d*g + g*h + g*i + g*l
        - h*l - i*j - j*k,
    a^2 - a*b + a*d + a*f - 2*a*h + a*i + 
        2*a*j - a*k + b^2 - b*d + b*f + b*g +
        b*i - b*j - c^2 + c*e - c*i - c*j + 
        d*g - g*i - g*k - h^2 + i^2 - i*k + 
        i*l + j*l,
    a^2 + 2*a*b + a*c + a*d - a*e + a*f - 
        a*h + a*i + a*j + 2*a*k - a*l + b^2 
        - b*d - b*f + b*i + b*k - c^2 + c*d + 
        c*f + c*g - c*h + c*i - c*k - g*h +
        g*i - g*k - h*i + h*k + j*k + 
        k*l,
    a^2 + a*e + a*f + a*h + a*k + a*l - 
        2*b*c + b*d + b*e + b*f + b*h + b*j
        + c*d + c*h - c*i + f*k + g*i + 
        g*k + h^2 - i^2 - i*k - j*k - k^2 + 
        k*l,
    -2*a^2 - a*b + a*d + a*f - a*g - a*h + 
        a*i - a*l - b^2 + b*d + b*f - b*h - 
        b*i - b*k + 2*c*e + c*f - c*g - c*h
        - c*j - f*k - g*i - g*l - 2*h^2 + 
        h*i + h*j + i^2 + i*j - i*k + i*l +
        j*k + k^2 - k*l,
    a*c - 2*a*d + a*e - a*g - 2*a*i - a*j +
        b^2 - b*c + b*d + b*f - 2*b*g + 2*b*h +
        b*j + c*f + c*h + c*i + c*j + c*k
        + g^2 - g*i - g*j + i*k - j^2 + j*k -
        j*l,
    -2*a*b + a*c + a*d + 2*a*g - a*h + 2*a*i
        + a*k + a*l - 2*b*c - b*d + b*g - 
        b*h + b*k - c*d - c*f + c*g + c*h -
        2*c*i - c*j + f*l + h^2 - h*j - i*k
        - j*k - k^2 + k*l,
    -a^2 - a*b - a*h - 2*a*i - 2*a*k + 2*b*c +
        b*d - b*e + b*f - b*g + 2*b*h - b*i 
        - c^2 + c*d - c*e + c*h - c*i - c*j + 
        c*l + g*k + h*i + j*l + k*l,
    a^2 + a*c - a*e + a*f + a*g - a*h + 
        a*i + a*j - a*k + 2*b^2 + b*c - 2*b*e
        + b*h + b*i - b*j + b*k - 2*c^2 + c*e
        + c*g + e*j + g*h - h^2 + i^2 - j*l,
    -a*b - a*c + a*e + a*f + a*g - a*h + 
        a*j - a*k + a*l - b*c + b*h + b*i
        - c*d - c*i - c*j + e*j + f*j + 
        g*i + g*j + g*l + h*k - h*l - 
        i*j - j*k + j*l,
    a^2 + 2*a*b + a*c - a*d - a*e - 2*a*f + 
        2*a*h - a*i - a*j + a*k + b*c - 
        b*d + b*e - b*f - b*h - b*i + c*i + 
        c*j - c*l - d*g + f*i - g*j + h^2 + 
        h*l - i^2 + i*k - i*l - j*l,
    a^2 - a*d - a*f + 2*a*h - a*i - a*j + 
        a*l - b^2 - b*c + b*d + b*e - b*i + 
        c*k - c*l + f*g + h^2 - i^2 + i*k - 
        i*l - j*l,
    a^2 + a*c - a*h - a*j + a*k + 2*b^2 + 
        b*c - b*d - b*f + b*h - c*e + 2*c*g 
        - c*h + c*i + c*j - c*k + c*l + 
        e*f - g*h - g*l - h*i + i^2 + k*l,
    -a*b - a*c - a*d + a*e - a*f - a*g - 
        a*i - a*j - a*k - 2*b^2 + b*c - b*g +
        b*h - 2*b*i + b*j - 2*b*k + c^2 - c*e
        + c*f + d*g + f*h - g*h + g*k - 
        g*l + j*k + k*l,
    a*c + a*d - a*e - a*i - a*j - a*k - 
        b^2 + b*e + 2*b*f + b*g - b*h - b*j + 
        c*e - c*h + c*i + c*j + c*k + e*f 
        + f^2 - h^2 + h*j - h*l + i^2 + i*l - 
        j^2,
    a^2 - a*b - a*c - a*g + a*h - a*j + 
        a*k - b^2 + b*e - b*h - b*i + 2*c*d + 
        c*e + c*f + e*h - h*j + h*k - i^2 + 
        i*j - i*k,
    a*b + a*c + a*d - a*e - a*f + a*h + 
        a*i - 2*b^2 + 2*b*c + b*g - b*h - b*j 
        - c^2 + c*d - c*f - c*g - c*i - c*k + 
        e*h + e*l + f*l + g*i + g*k + h*i
        - i^2 - i*k + j^2 - j*k,
    -2*a*b - a*e - a*f + a*g - a*i - 2*a*j 
        + a*k + b*c + b*d - b*g - b*h - b*i
        + b*k + c^2 + c*g - c*h - c*k + c*l 
        - d*g + e*h + e*k + g*j - h*l - 
        i*j - j^2,
    2*a*b + a*d + a*h - a*i - a*l + b^2 + 
        b*c - b*d + b*e - b*h + 2*b*i + b*j
        + c*d + c*e - 2*c*g + c*i + c*j + 
        c*k + e*g + g*l + i*j + i*l - 
        k*l,
    -a^2 - a*c + a*e + a*f + a*k + b^2 - 
        b*d - b*f + b*h - 2*c*d + c*e - c*k
        - d*g + d*l + h*i + h*k + k^2 - 
        k*l,
    a^2 + a*b + a*e + a*f - a*h + 2*a*k + 
        b^2 - b*c - b*d + b*e + b*h + b*k + 
        c*g - c*k - d*g + d*k,
    -a^2 - a*b + a*e + a*f - a*g + a*i + 
        a*j + a*k + 2*b^2 - b*c + b*d - b*h +
        b*j - b*k + c^2 - c*f - c*j + d*j -
        g*k - h*j + h*k - j^2 + k^2 - k*l,
    a^2 + a*b + a*e + a*f - a*h + a*j + 
        a*k + b^2 - b*c - b*d + b*e + b*h + 
        b*i + b*j + b*l,
    a*e + a*h - 2*a*k + a*l - b^2 + b*d + 
        b*e + 2*b*f + b*h + b*i - b*k - c*f
        - c*g - c*i + c*k + d*g + d*i + 
        g*k - h*k - i*k + i*l - j*k,
    a*b - a*d + a*e + a*f - a*g + a*i + 
        b^2 - b*c + b*d - b*f - b*h + b*i - 
        b*k + c^2 - c*e + c*g - 2*c*h + c*i - 
        c*l + d*h - g*k - h*i + i*k + 
        j*k,
    a^2 - a*b + a*g - a*j - a*k + a*l - 
        2*b^2 - b*c + b*e + b*f + b*k - c*d + 
        c*g - c*h + d*e + g*j + h*j - h*k
        - i*j - k^2 + k*l,
    -a^2 + a*b + a*c - a*d - a*e - a*f - 
        a*h + a*i - a*l + b*c - b*d - b*e -
        b*f - b*g + b*k + c*e + c*f + c*j 
        + c*l + e*i + f*i - g*i + g*k - 
        g*l - h^2 - h*k + h*l + i^2 + i*j + 
        j*k - j*l,
    2*a^2 - a*e - a*f + a*h - 2*a*i - 2*a*j + 
        b^2 + 2*b*d + b*f - b*i - c^2 + c*d + 
        c*e + c*h - c*i - c*j + d*f + h^2 + 
        h*l - i^2 - i*l - k^2 + k*l,
    -a^2 + a*c - a*d - a*f - a*g - a*j - 
        a*k - a*l - 2*b^2 + b*c + b*e + b*f -
        b*h + b*i + b*j - b*k + c^2 - c*g + 
        c*j + c*l + d^2 + d*f + d*g + g*k - 
        h^2,
    a^2 - a*b + 2*a*h - a*i - 2*a*k - b^2 + 
        b*c + 2*b*d + b*e + 2*b*f + b*g - 
        b*i - 2*b*j - b*k - 2*c^2 + c*e - c*h
        - c*i - c*j - c*k - c*l + d*e - e^2 
        - e*f + f*l + g*i - g*k + g*l + 
        h*k - h*l + i*l + j^2 + j*l
];
Xns27 := Curve(ProjectiveSpace(Rationals(),11),psi : Nonsingular := true);

K<zeta> := CyclotomicField(3);
Xns27K := BaseChange(Xns27,K);
RR<A,B,C> := PolynomialRing(K,3);
pol := A^4 + (zeta - 1)*A^3*B + (3*zeta + 2)*A^3*C - 3*A^2*C^2 +
    (2*zeta + 2)*A*B^3 - 3*zeta*A*B^2*C + 3*zeta*A*B*C^2 - 2*zeta*A*C^3 -
    zeta*B^3*C + 3*zeta*B^2*C^2 + (-zeta + 1)*B*C^3 + (zeta + 1)*C^4;
XH := Curve(ProjectiveSpace(K,2),pol);

// Staring at the Fourier expansions of modular forms that David Zywina's code produces 
// and matching those with the modular forms for the subgroup H whose determinant isn't
// surjective yields the following map.

phi := map<Xns27K -> XH | [(1/3)*((-2-zeta)*Xns27K.11 + (1-zeta)*Xns27K.12), 
(-1/3)*((1+2*zeta)*Xns27K.9 + (-1+zeta)*Xns27K.10),
(-1/3)*((1+2*zeta)*Xns27K.7 + (1-zeta)*Xns27K.8)]>;
u:=zeta;

//Changing the coordinates back to RSZB format, which we had initially changed to use Tuitman's algorithm.
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
  [(1/2)*(zeta-3),(1/2)*(zeta+2),1] // j = (-3238903991430*zeta - 4786881515880)^3/1/9^27, non-CM
]];

CM_points:=[points[i]: i in [1..12]];
non_CMpoint:=points[13];
non_CMns27:=RationalPoints(non_CMpoint@@phi);
//CMj := [0,1728,-3375,8000,-32768,54000,287496,-884736,-12288000,16581375,-884736000,-147197952000,-262537412640768000];
Rationalpts:=[**];
Kpoints:=[**];
KpointsCM:=[**];

x:=0;
for pt in points do
    pt_preimage_list := RationalPoints(pt@@phi);
    for ptinv in pt_preimage_list do
        Append(~Kpoints,ptinv);
        try
            rtlptinv:=Xns27![Rationals()!ptinv[i]: i in [1..12]];
            Append(~Rationalpts,rtlptinv);
            catch e
                x:=x+1;
        end try;            
    end for;    
end for;   

for pt in CM_points do
    pt_preimage_list := RationalPoints(pt@@phi);
    for ptinv in pt_preimage_list do
        Append(~KpointsCM,ptinv);      
    end for;    
end for;   

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




printf "The number of Q(zeta_3) points of X_ns^+(27) is %o, and of these %o are Q-points.\n", #Kpoints,#Rationalpts;
printf "X_{ns}^+(27) has %o Q(zeta_3)-points which are CM. \n", #KpointsCM;
printf "There are %o non-CM point (s), which are \n %o.\n" , #non_CMns27, non_CMns27;
// jinvlist:=[]; 
// x:=0 ; //CMcount
 
// pt:=points[1];
// jinvar := phi3(phi2(pt))[1];
// chk, disc := HasComplexMultiplication(EllipticCurveWithjInvariant(jinvar));
//      if chk  then
//         x:=x+1;
//         Append(~jinvlist, jinvar);
//      else
//         Append(~jinvlist, jinvar);
//      end if;    

// RF:=recformat<pt:Pt, jinv: FldCycElt>;
// jinvlist:=rec<RF|>;
    

// for pt in points do
//     pt_preimage := RationalPoints(pt@@phi);
// //   end for;
// //   printf "The points in XH(K) above j = %o are %o.\n",curj,X3list;
// //   printf "The points on Xns+(27)(K) above j = %o are %o.\n",curj,X4list;
// end for;  