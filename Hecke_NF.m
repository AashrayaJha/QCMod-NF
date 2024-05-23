AttachSpec("QCMod.spec");

//load coleman data for sufficiently high precision. 
load "data/Correspondence_31_300.m";

import "misc.m": algdepQp,lindepQp, alg_approx_Qp;
data_1:=data_1;
//data_2:=data_2;

//initially did it for both embeddings, in hopes that one could recover the numbers with trace and norms.
//This is a good check, but first just find the correspondence.

Q:=data_1`Q; g:=data_1`g; d:=Degree(Q); p:=data_1`p; v:=data_1`v; 
q:=p; N:=data_1`N;
Qp:=pAdicField(p,N);

K:=BaseRing(BaseRing(Q));

F1 := data_1`F;
if q eq p then F1 := Submatrix(data_1`F,1,1,2*g,2*g); end if;// Necessary when q=p
F1inv := Transpose(F1)^(-1);
Aq_1 := Transpose(F1)+q*F1inv;   // Eichler-Shimura -> Hecke operator
prec_loss_bd1 := Valuation(Norm(Determinant(F1inv)), p);

AK := ZeroMatrix(K, 2*g, 2*g); 
bad_indices:=[**];

//Hecke correspondence
for j in [1..2*g] do
    for k in [1..2*g] do
        try
            AK[j,k] := alg_approx_Qp(Qp!Rationals()!Aq_1[j,k],v);    // recognition of integer in Zp via LLL
        catch e
            Append(~bad_indices,[j,k]);
        end try;                   
    end for;
end for;
//Finish Hecke correpondence. 

// Start computing nice correspondences 
C:=ZeroMatrix(K,2*g,2*g);
for i:=1 to g do
C[i,g+i]:=1;
C[g+i,i]:=-1; 
end for;
Z1:=ZeroMatrix(K,2*g,2*g); Z2:=ZeroMatrix(K,2*g,2*g);
Zs:=[Z1,Z2];

for i in [1..2] do 
    A:=Aq_1^i;
    Zmx := (2*g*Aq_1^i-Trace(A)*IdentityMatrix(K,2*g))*C^(-1);
    for j in [1..2*g] do
        for k in [1..2*g] do
            try
                Zs[i][j,k] := alg_approx_Qp(Qp!Rationals()!Zmx[j,k],v);    // recognition of integer in Zp via LLL
            catch e
                Append(~bad_indices,[j,k]);
            end try;                   
        end for;
    end for;
end for;

//print results
out:=Sprintf("AK_patch_3_400:=%m;",AK);
output_file:="data/New_hecke_good_patch.m";
out_Zs :=Sprintf("Zs_patch_3_400:=%m;", Zs) ;
Write(output_file,out);
Write(output_file,out_Zs);

//below are functions for verifying results via trace and norm using both embeddings. 

function make_poly_quadratic(trace,norm)
    R<x>:=PolynomialRing(Rationals());
return x^2-trace*x+norm;
end function;

function Matrix_from_trace_norm(trace,norm,A,v)
// If the trace and norm matrices have been recognised as elements over the rationals, reconstruct
// the matrix of the original A with the minimal polynomial using the trace and norm matrices.
K := NumberField(Order(v));
Kv, loc := Completion(K, v);
n:= #Rows(trace);
M:=ZeroMatrix(K,n,n);
bad_indices_1:=[**];
for i in [1..n] do
    for j in [1..n] do
        Tr:=trace[i][j];
        Nm:=norm[i][j];
        a_Qp:=A[i][j];
        
        poly:=make_poly_quadratic(Tr,Nm);
        poly:=ChangeRing(poly,K);
        alist := [-Coefficient(fac[1], 0)/Coefficient(fac[1], 1) : fac in Factorization(poly) | Degree(fac[1]) eq 1];
        try 
            // Sort roots by how close they are in Qp to the original value
            Sort(~alist, func< a, b | Valuation(loc(a_Qp) - Qp!loc(b)) - Valuation(loc(a_Qp) - Qp!loc(a)) >);
            M[i][j]:=alist[1]; // Sort roots by how close they are in Qp to the original value
        catch e
            Append(~bad_indices_1,[i,j]);
        end try;    
    end for ;
end for;

return M;
end function;

// F2 := data_2`F;
// if q eq p then F2 := Submatrix(data_2`F,1,1,2*g,2*g); end if;// Necessary when q=p
// F2inv := Transpose(F2)^(-1);
// Aq_2 := Transpose(F2)+q*F2inv;   // Eichler-Shimura -> Hecke operator
// prec_loss_bd2 := Valuation(Norm(Determinant(F2inv)), p);
  
// A_sum := Aq_1+Aq_2;
// A_prod:=ZeroMatrix(Rationals(),6,6); 
// for i in [1..6] do
//     for j in [1..6] do
//         A_prod[i][j]:=Aq_1[i][j]*Aq_2[i][j];
//     end for;
// end for;     

// [A_prod1[i][j]:=Aq_1[i][j]*Aq_2[i][j]: i in [1..6] and j in [1..6]];
// Alist:=[*A_sum,A_prod*];
// A_Qsum:=ZeroMatrix(Rationals(),6,6);
// A_Qprod:=ZeroMatrix(Rationals(),6,6);
// A_Qlist:=[A_Qsum,A_Qprod];
// bad_indices_sp:=[[**]:i in [1..2]];


// for j in [1..2*g] do
//     for k in [1..2*g] do
//         try 
//             A_Qprod[j,k] := lindepQp(pAdicField(q, Floor(N-1))!Rationals()!A_prod[j,k]);    // recognition of integer in Zp via LLL
//         catch e  
//             Append(~bad_indices_sp[1],[j,k]);
//         end try;    
//     end for;
// end for;

// for j in [1..2*g] do
//     for k in [1..2*g] do
//         try 
//             A_Qsum[j,k] := lindepQp(pAdicField(q, Floor(N-1))!Rationals()!A_sum[j,k]);    // recognition of integer in Zp via LLL
//         catch e  
//             Append(~bad_indices_sp[2],[j,k]);
//         end try;    
//     end for;
// end for;

// Check:=[[**]: i in [1..2]];   

// for i in [1..6] do
//     for j in [1..6] do
//         sum_1:=Trace(AK[i][j]);
//         prod_1:=Norm(AK[i][j]);
//         sum_2:=Rationals()!A_Qsum[i,j];
//         prod_2:=Rationals()!A_Qprod[i,j];
//         if Rationals()!sum_1 eq Rationals()!sum_2 then 
//             continue ;
//         else     
//             Append(~Check[1],[i,j,sum_1/sum_2]);
//         end if;
//         if Rationals()!prod_1 eq Rationals()!prod_2 then 
//             continue ;
//         else     
//             Append(~Check[2],[i,j,prod_1/prod_2]);
//         end if;
//     end for;
// end for;

// out_prod:=Sprintf("AQ_prod_40:=%m;",A_Qprod);
// out_sum:=Sprintf("AQ_sum_40:=%m;",A_Qsum);

//Write(output_file,out_sum);
//Write(output_file,out_prod);

