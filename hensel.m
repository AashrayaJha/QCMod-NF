function hensel_lift_n(flist,p,prec)

// porting Francesca Bianchi's code
//https://github.com/bianchifrancesca/QC_elliptic_imaginary_quadratic_rank_2/blob/master/auxiliary_functions.sage
// which says
// Multivariable Hensel lifter for roots that are simple modulo `p`.
// This is essentially the code from [S15] with some minor modifications.
// To do: was lazy in one part and made it for n = 2 (fix this)
// [S15]: \B. Schmidt, "Solutions to Systems of Multivariate p-adic Power Series". Oxford MSc Thesis, 2015.
/*

**work in progress**

Example

R<s, t> := PolynomialRing(pAdicField(5, 10),2);
f1 := s + t - 2*s*t;
f2 := s - t;
a, b := hensel_lift_n([f1, f2], 5, 10);
*/

precvec := [];
k := 1;
for F in flist do
    R1 := Parent(flist[k]);
    F1 := BaseRing(R1);
    if IsExactpAdic(F1) then
        precision1 := prec;
    else
        precision1 := Precision(F1);
        if prec gt Precision(F1) then
            print "Cannot get %o  digits of precision due to precision of inputs of f1; raise precision of inputs", prec;
        elif prec lt Precision(F1) then
            precision1 := prec;
        end if;
    end if;
    Append(~precvec, precision1);
    k := k + 1;
end for;
precision := Min(precvec);

R := PolynomialRing(pAdicField(p,precision), #flist);
flistnew:=[**];
for F in flist do
    Append(~flistnew,R!F);
end for;
Jlist:=[];
for F in flistnew do
    for i in [1..#flistnew] do
        Append(~Jlist, Derivative(F,i));
    end for;
end for;
J := Matrix(#flistnew, #flistnew, Jlist);
M := Determinant(J);
if #flistnew eq 2 then
    coords:=[[i,j] : i in [0..p-1], j in [0..p-1]];
else 
    return "To do!";
end if;
roots := [**];
roots_info := [**];
nonroots := 0;
for i in [1..#coords] do
    valuesval := [Valuation(Evaluate(F,coords[i])): F in flistnew];
    min_valuesval := Minimum(valuesval);
    ord_det_J := Valuation(Evaluate(M,coords[i]));
    if Min(valuesval) gt 0 and ord_det_J eq 0 then
        roots := Append(roots, coords[i]);
        Append(~roots_info, [min_valuesval - 2*ord_det_J, ord_det_J]);
    elif min_valuesval gt 0 then
        nonroots:=nonroots+1;
    end if;
end for;

actual_roots := [**];

for r in roots do
    ind_roots := Index(roots, r);
    rt_info := roots_info[ind_roots];
    print rt_info[1];
    if Type(rt_info[1]) eq Intrinsic then
        Append(~actual_roots, Matrix(#flist, 1, r));
    else
        variables := [];
        k := 0;
        i_l := Matrix(#flist, 1, r);
        Jeval := Matrix(#flistnew, #flistnew, [Evaluate(f,r): f in Jlist]);
        B:= Transpose(Jeval)*Jeval;
        rt_info[2];
        rt_info[1];
        const1:=Ceiling(Log( ((prec-rt_info[2])/rt_info[1]))/Log(2.)) + 1 ;
        while k lt const1 and Determinant(B) ne 0 do
            A := Matrix(#flistnew, 1, [-Evaluate(f,r): f in flistnew]);
            i_l := i_l + B^(-1)*Transpose(Jeval)*A;
            for i in [1..#flist] do
                Append(~variables, i_l[i, 1]);
            end for;
            Jeval := Matrix(#flistnew, #flistnew, [Evaluate(f,variables) : f in Jlist]);
            variables := [**];
            k := k+1;
            BB:= Transpose(Jeval)*Jeval;
        end while;
        Append(~actual_roots, i_l);
      end if;
end for;

return actual_roots, nonroots;      
end function;
