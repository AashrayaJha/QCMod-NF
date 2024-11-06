P2<x,y,z> := ProjectiveSpace(Rationals(),2);
 C := Curve(P2, x^4+y^4+3*y^3*z-2*z^4);
 D := Curve(P2, (x+z)^4+y^4+3*y^3*z-2*z^4);
 b, phi := IsIsomorphicPlaneQuartics(C,D);
 phi;

 K<u> := CyclotomicField(3);  
 A2<x,y> := AffineSpace(K, 2);
 Q1 := y^4 + ((-2*u + 9)*x + (2*u + 3))*y^3 + (-3*x^2 + 6*x - 3)*y^2 + ((-170*u + 254)*x^3 + (-150*u + 114)*x^2 + (-54*u + 18)*x - 10*u - 2)*y + (162*u + 144)*x^4 + (-108*u + 48)*x^3 + (-72*u - 144)*x^2 + (12*u - 48)*x + 6*u;
 Q2 := y^4 + ((u - 1)*x + (3*u + 2))*y^3 - 3*y^2 + ((2*u + 2)*x^3 - 3*u*x^2 + 3*u*x - 2*u)*y - u*x^3 + 3*u*x^2 + (-u + 1)*x + u + 1;
 C1 := ProjectiveClosure(Curve(A2, Q1));
 C2 := ProjectiveClosure(Curve(A2, Q2));
 b, phi := IsIsomorphicPlaneQuartics(C1,C2);
 phi;

IsIsomorphicPlaneQuartics:Maximal;
