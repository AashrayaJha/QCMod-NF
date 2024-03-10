//This file will calculate the Trace zero correspondences which we need to run the QC algorithm.

AttachSpec("~/GitHub/CHIMP/CHIMp.spec");
AttachSpec("QCMod.spec");
//load "data/NF-example-coleman-data-13_160.m"; 

//10/3/2024 AJ- I think the correspondences obtained here are correct. The number of digits seem
// to have stabilised. Also they satisfy the right minpoly, though this is a abd check
// Just going to try and work directly with the K-numbers now. 

//import "misc.m": algdepQp,lindepQp, alg_approx_Qp;

data_1:=data_1;

Q:=data_1`Q; g:=data_1`g; d:=Degree(Q); p:=data_1`p; v:=data_1`v; 
q:=p; N:=data_1`N;
Qp:=pAdicField(p,N);

K:=BaseRing(BaseRing(Q));
C:=ZeroMatrix(K,2*g,2*g);
//C is the sympectic matrix.
for i:=1 to g do 
    C[i,g+i]:=1;
    C[g+i,i]:=-1; 
end for;

F1 := data_1`F;
if q eq p then F1 := Submatrix(data_1`F,1,1,2*g,2*g); end if;// Necessary when q=p
F1inv := Transpose(F1)^(-1);
Aq_1 := Transpose(F1)+q*F1inv;   // Eichler-Shimura -> Hecke operator
prec_loss_bd1 := Valuation(Norm(Determinant(F1inv)), p);

AK := ZeroMatrix(K, 2*g, 2*g); ZK := AK;
bad_indices:=[[*[**],[**]*]:i in [1..g-1]];
Zs:=[]; As:=[];

for i in [1..g-1] do
      A := Aq_1^i; // ith power of hecke operator
      Zmx := (2*g*A-Trace(A)*IdentityMatrix(K,2*g))*C^(-1);      // Zmx is a q-adic approximation of a nice correspondence Z, ie trace 0
      //The denominator created problems earlier, let's see if this works without this. 

    for j in [1..2*g] do
        for k in [1..2*g] do
    
            try
                AK[j,k] := alg_approx_Qp(Qp!Rationals()!Aq_1[j,k],v);    // recognition of integer in Zp via LLL
            catch e
                Append(~bad_indices[i][1],[j,k]);
            end try;     
            
            try 
            ZK[j,k] := alg_approx_Qp(Qp!Rationals()!Zmx[j,k],v);  // dito
            catch e
                Append(~bad_indices[i][2],[j,k]);
            end try;   

        end for;
    end for;
    if Trace(ZK*C) ne 0 then // approximation issue. Perturbe ZQ slightly.
        if q ne p then 
            error "q-adic approximation of nice correspondence not exact.";  
        end if;  
        W:=ZeroMatrix(K,2*g, 2*g);
        W[1,g+1]:=Trace(ZK*C);
        W[g+1,1]:=-Trace(ZK*C);
        ZK:=2*ZK+W;
    end if;
    Append(~Zs,ZK);
    Append(~As,AK);
end for;
