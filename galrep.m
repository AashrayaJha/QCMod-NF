/*
    Magma code related to the paper :
        "Computing images of Galois representations attached to elliptic curves", Forum of Mathematics, Sigma 4 (2016) e4 (79 pages), http://dx.doi.org/10.1017/fms.2015.33
    Copyright (c) 2015 by Andrew V. Sutherland
    You are welcome to use and/or modify this code under GPL version 2 or any later version (but please be sure to cite the paper)
   
    This module naively implements the function GaloisImage(E,m), which, given an elliptic curve E/Q and an integer m > 1 computes
    a subgroup G_E(m) of GL_2(Z/mZ) conjugate to the image of the mod-m representation.  It does this by explicitly computing the Galois action of the m-torsion field.
    
    This code is painfully slow, but if you are patient it will do the job (for small m).
*/

// Returns poly who's roots are precisely the x-coords of the points of order m on E (this means removing factors that are n-division polys for n dividing m)
PrimitiveDivisionPolynomial := function(E,m)
    local f;
    f:=DivisionPolynomial(E,m);
    for d in Divisors(m) do if d gt 1 and d lt m then f := ExactQuotient(f,$$(E,d)); end if; end for;
    return f;
end function;

// Magma really wants number fields to be defined by integral monic polynomials, so we make sure this happens
IsIntegrallyDefined := function(K)
    local f;
    if K eq Rationals() then return true; end if;
    if not IsAbsoluteField(K) then return false; end if;
    f := DefiningPolynomial(K);
    return IsMonic(f) and &and[c in Integers():c in Coefficients(f)];
end function;

// Redefines a number field so that it is defined in terms of the absolute minimal polynomial of a generator that is an algebraic integer
MakeIntegrallyDefined := function(K)
    local g;
    while not IsIntegrallyDefined(K) do
        g := SimpleExtension(K).1;
        f := MinimalPolynomial(g);
        g *:= &*PrimeDivisors(LCM([Denominator(c/LeadingCoefficient(f)):c in Coefficients(f)]));
        K:=NumberField(MinimalPolynomial(g));
    end while;
    return K;
end function;


// Returns a pair [P,Q] of independent generators for E[m] (the points P and Q will necessarily have order m), where E is an elliptic curve y^2=x^3+Ax+B with A,B in Q
// Be warned that this is painfully slow: unless the m-division field Q(E[m]) has very small degree you will need to be patient.
TorsionField := function(E,m)
    local C, K, L, EL, x1, y1, y1s, x2, y2, y2s, Q, P, S, phi, f, b, g;

    C := Coefficients(E);
    assert C[1] eq 0 and C[2] eq 0 and C[3] eq 0 and C[4] in Rationals() and C[5] in Rationals(); // To simplify matters, we require E to be in the form y^2=x^3+Ax+B with A,B in Q
    phi:=PrimitiveDivisionPolynomial(E,m);
    roots := Roots(phi);
    if #roots ne Degree(phi) then
        K:=SplittingField(phi);
        return $$(ChangeRing(E,MakeIntegrallyDefined(K)),m);
    end if;
    K:=BaseRing(E);
    L:=K;
    R<x>:=PolynomialRing(K);
    // Our first basis point P (of order m) will have x-coord equal to the first root of phi
    x1:=roots[1][1];
    f:=x^3+C[4]*x+C[5];
    y1s:=Evaluate(f,x1);
    b,y1:=IsSquare(y1s);  // this step is time-consuming
    // if y1 is not in L, extend L so that it is
    if not b then L := NumberField(x^2-y1s); end if;
    if L ne Rationals() and not IsAbsoluteField(L) then L:=AbsoluteField(L); end if;
    return MakeIntegrallyDefined(L);
end function;

    
// Returns a pair [P,Q] of independent generators for E[m] (the points P and Q will necessarily have order m), where E is an elliptic curve y^2=x^3+Ax+B with A,B in Q
// Be warned that this is painfully slow: unless the m-division field Q(E[m]) has very small degree you will need to be patient.
TorsionBasis := function(E,m)
    local C, K, L, EL, x1, y1, y1s, x2, y2, y2s, Q, P, S, phi, f, b, g;

    C := Coefficients(E);
    assert C[1] eq 0 and C[2] eq 0 and C[3] eq 0 and C[4] in Rationals() and C[5] in Rationals(); // To simplify matters, we require E to be in the form y^2=x^3+Ax+B with A,B in Q
    phi:=PrimitiveDivisionPolynomial(E,m);
    roots := Roots(phi);
    if #roots ne Degree(phi) then
        K:=SplittingField(phi);
        return $$(ChangeRing(E,MakeIntegrallyDefined(K)),m);
    end if;
    K:=BaseRing(E);
    L:=K;
    R<x>:=PolynomialRing(K);
    // Our first basis point P (of order m) will have x-coord equal to the first root of phi
    x1:=roots[1][1];
    f:=x^3+C[4]*x+C[5];
    y1s:=Evaluate(f,x1);
    b,y1:=IsSquare(y1s);  // this step is time-consuming
    // if y1 is not in L, extend L so that it is
    if not b then L := NumberField(x^2-y1s);  y1:=L.1; x1 := L!x1; f:= ChangeRing(f,L); E:=ChangeRing(E,L); end if;
    // make a list S of the x-coords of the points in <P> .  Note that we only need to compute multiples of P up to m/2 since -P and P have the same x-coord.
    S:=[x1];
    P:=E![x1,y1];
    Q := P+P;
    for i:=2 to Floor(m/2) do
        S:=Append(S,Q[1]);
        Q +:= P;
    end for;
    // find a root of phi that is not the x-coord of a point in S
    for r in roots do
        if r[1] in S then continue; end if;
        // Construct P2 not in <P1>. 
        x2 := r[1];
        y2s:=Evaluate(f,x2);
        b,y2:=IsSquare(y2s); // this step is *really* time consuming
        assert b; // We are guaranteed that P2 is in E(L), since its x-coord is, and so is the y-coord of P1 (see lemma in galrep paper)
        if not IsPrime(m) then // if m is not prime then we also need to verify that no multiple of Q=[x2,y2] lies in <P>
            EL:=ChangeRing(E,L);
            Q:=EL![x2,y2];
            R:=EL!0;
            fail := false;
            for i:= 1 to Floor(m/2) do
                R+:=Q;
                if R[1] in Parent(x1) and Parent(x1)!R[1] in S then fail := true; end if;
            end for;
            if fail then continue; end if;
        end if;
        break;
    end for;
    if L ne Rationals() and not IsAbsoluteField(L) then L:=AbsoluteField(L); end if;
    if not IsIntegrallyDefined(L) then return $$(ChangeRing(E,MakeIntegrallyDefined(L)),m); end if;
    EL:=ChangeRing(E,L);
    return [EL![x1,y1],EL![x2,y2]];     // note that this line will fail if x2 and y2 did not get set above
end function;

// Given a basis B for E[m], returns a list A of the points in E[m] ordered so that A[i*m+j+1] = i*B[1]+j*B[2] for i,j in [0,m-1]
TorsionPoints := function(B,m)
    local A;
    A:=[Parent(B[1])!0:i in [1..m^2]];
    for j:= 1 to m-1 do A[j+1] := A[j]+B[2]; end for;
    for i:= 1 to m-1 do
        A[i*m+1] := A[(i-1)*m+1]+B[1];
        for j:= 1 to m-1 do A[i*m+j+1] := A[i*m+j]+B[2]; end for;
    end for;
    return A;
end function;

// Given a list A of E[m] ordered so that A[i*m+j+1] = i*B[1]+j*B[2] for some basis B, and a point P, returns i and j such that P=i*B[1]+j*B[2]
TorsionPointIndex := function(A,P,m)
    k := Index(A,P);
    assert k ne 0;
    j := (k-1) mod m;
    i := ExactQuotient(k-1-j,m);
    return [i,j];   
end function;

// Given sigma in Gal(Q(E[m]/Q) and P in E[m], returns sigma(P)
Action := function(sigma,P);
    return Parent(P)![sigma(P[1]),sigma(P[2]),sigma(P[3])];
end function;

// Given sigma in Gal(Q(E[m]/Q), a basis B for E[m], a list A of E[m] ordered as above, returns a matrix M in GL(2,Z/mZ) such that M*B = sigma(B), where B is viewed as a column vector
SigmaMatrix := function(sigma,B,A,m)
    return Transpose(GL(2,Integers(m))![TorsionPointIndex(A,Action(sigma,B[1]),m),TorsionPointIndex(A,Action(sigma,B[2]),m)]);
end function;


// Computes the mod m Galois image of E as a subgroup of GL(2,Z/mZ) with respect to the basis B
GaloisImage:= function(E,B,m)
    local t0, t1, G,S,phi,A;
    //assert BaseRing(E) eq Rationals();  // we only handle curves over Q (already hard enough)
    //t0:=Cputime();
    //printf "Computing basis for E[%o]...", m;
    //t1 := Cputime();
    //B:=TorsionBasis(WeierstrassModel(E),m);
    //printf "basis B constructed over degree-%o extension in %os\n", Degree(Parent(B[1][1])), Cputime()-t1;
    //printf "Computing Gal(Q(E[%o])/Q)...", m;
    t0 := Cputime();
    G,S,phi:=AutomorphismGroup(Parent(B[1][1]));
    printf "Galois group G of order %o computed in %os\n", #G, Cputime()-t0;
    printf "Enumerating E[%o] and computing Galois action...", m;
    t1:=Cputime();
    A:=TorsionPoints(B,m);
    GL2 := GL(2,Integers(m));
    H := sub<GL2|[SigmaMatrix(phi(g),B,A,m): g in Generators(G)]>;
    printf "%os\n", Cputime()-t1;
    printf "Galois image computes as an index %o subgroup of GL(2,Z.%oZ) in total time %os\n", Index(GL(2,Integers(m)), H), m, Cputime()-t0;
    return H, GL2;
end function;


// Below are a bunch of small examples to demonstrate the functionality.
// None of them should individiually take more than a few minutes except for the last one
/*
print "Use GalrepExamples() to see some examples.";
GalrepExamples := procedure()

D:=CremonaDatabase();

m:=2;
E:=EllipticCurve(D,15,1,1);
G,B:=GaloisImage(E,m);
printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;
E:=EllipticCurve(D,14,1,1);
G,B:=GaloisImage(E,m);
printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;
E:=EllipticCurve(D,196,1,1);
G,B:=GaloisImage(E,m);
printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;

m:=3;
E:=EllipticCurve(D,14,1,1);
G,B:=GaloisImage(E,m);
printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;
E:=EllipticCurve(D,98,1,3);
G,B:=GaloisImage(E,m);
printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;
E:=EllipticCurve(D,14,1,4);
G,B:=GaloisImage(E,m);
printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;
E:=EllipticCurve(D,14,1,3);
G,B:=GaloisImage(E,m);
printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;
E:=EllipticCurve(D,338,4,1);
G,B:=GaloisImage(E,m);
printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;
E:=EllipticCurve(D,50,2,1);
G,B:=GaloisImage(E,m);
printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;
E:=EllipticCurve(D,245,1,1);
G,B:=GaloisImage(E,m);
printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;

m:=4;
E:=EllipticCurve([1,1,1,-10,-10]);
G,B:=GaloisImage(E,m);
printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;
E:=EllipticCurve([0,0,0,13,-34]);
G,B:=GaloisImage(E,m);
printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;
E:=EllipticCurve([1,1,1,-135,-660]);
G,B:=GaloisImage(E,m);
printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;
E:=EllipticCurve([1,-1,1,4,-1]);
G,B:=GaloisImage(E,m);
printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;
E:=EllipticCurve([1,1,1,-80,242]);
G,B:=GaloisImage(E,m);
printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;
E:=EllipticCurve([1,-1,0,-6,8]);
G,B:=GaloisImage(E,m);
printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;
E:=EllipticCurve([0,1,0,4,4]);
G,B:=GaloisImage(E,m);
printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;
E:=EllipticCurve([0,1,0,-114,-127]);
G,B:=GaloisImage(E,m);
printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;
E:=EllipticCurve([1,0,1,4,-6]);
G,B:=GaloisImage(E,m);
printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;

m:=5;
E:=EllipticCurve(D,11,1,1);
G,B:=GaloisImage(E,m);
printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;
E:=EllipticCurve(D,275,2,2);
G,B:=GaloisImage(E,m);
printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;
E:=EllipticCurve(D,99,4,2);
G,B:=GaloisImage(E,m);
printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;
E:=EllipticCurve(D,6975,1,1);
G,B:=GaloisImage(E,m);
printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;
E:=EllipticCurve(D,11,1,3);
G,B:=GaloisImage(E,m);
printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;
// This curve hits an internal magma error
// E:=EllipticCurve(D,11,1,2);
// time G,B:=GaloisImage(E,m);
// printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;
E:=EllipticCurve(D,50,1,1);
G,B:=GaloisImage(E,m);
printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;
E:=EllipticCurve(D,50,1,3);
G,B:=GaloisImage(E,m);
printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;
E:=EllipticCurve(D,608,2,1);
G,B:=GaloisImage(E,m);
printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;

m:=6;
E:=EllipticCurve([1,0,1,4,-6]);
G,B:=GaloisImage(E,m);
printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;
E:=EllipticCurve([1,0,1,-19,26]);
G,B:=GaloisImage(E,m);
printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;
E:=EllipticCurve([1,1,0,-1740,22184]);
G,B:=GaloisImage(E,m);
printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;
E:=EllipticCurve([1,0,1,-1,0]);
G,B:=GaloisImage(E,m);
printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;
E:=EllipticCurve([1,1,0,-475,-4187]);
G,B:=GaloisImage(E,m);
printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;
E:=EllipticCurve([0,1,0,-114,-127]);
G,B:=GaloisImage(E,m);
printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;
E:=EllipticCurve([0,-1,0,-1,0]);
G,B:=GaloisImage(E,m);
printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;

m:=7;
E:=EllipticCurve(D,2450,27,1);
G,B:=GaloisImage(E,m);
printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;

m:=8;
E:=EllipticCurve([1,1,1,-10,-10]);
G,B:=GaloisImage(E,m);
printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;

print "The last example typically takes 5 minutes, press ctrl-C if you don't want to wait...";

// This example causes magma to run out of memory
//m:=9;
//E:=EllipticCurve(D,17100,7,1);
//G,B:=GaloisImage(E,m);
//printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;

m:=10;
E:=EllipticCurve([0,-1,1,-10,-20]);
G,B:=GaloisImage(E,m);
printf "Elliptic curve E=%o has Galois image G_E(%o) = %o\n\n",CremonaReference(E), m, G;

end procedure;
*/
