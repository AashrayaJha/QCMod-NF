import "misc.m": algdepQp,lindepQp, alg_approx_Qp;

intrinsic HeckeCorrespondenceNF(data::Rec, N::RngIntElt) 
    -> AlgMatElt, SeqEnum, RngIntElt, SeqEnum
 
  {For i=1,...,g-1, construct a nice correspondence Zi from the ith power of
   the Hecke operator Ap using Eichler-Shimura.  N is the optional precision for the p-adic computation. 

   Input: Record which is the output of ColemanData in singleintegrals.m. 
   Warning: The p-adic precision to which one needs the Frobenius action was very high (N=344) for the 
   for the curve X in Examples/qc_Xns27quotient.m. This took about 4 hours to compute on an Apple M2 with 16 GB RAM. 

   TODO: Change Hodge and Frobenius intrinsics so that it just takes p-adic matrcies as input, and hence won't need 
   algdep and high p-adic precision. 

   Ouput: Tp, Zs= [*Z1,..Z_g-1*], prec_loss_bd, bad_indices ,where Tp is the Hecke operator for the prime above p 
   for which
   the data was cacluated, Zi=2*g*Tp^i-Trace(Tp^i), prec_loss_bd is the loss in precision expected, bad_indices are the indices
   where LLL failed to recognise the p-adic approximations of the numbers.

   Tp and Zi are encoded as matrices representing their action on H^1_dR(X).}

    Q:=data`Q; g:=data`g; d:=Degree(Q); p:=data`p; 

    Q:=data`Q; g:=data`g; d:=Degree(Q); p:=data`p; v:=data`v; 
    q:=p; N:=data`N;
    Qp:=pAdicField(p,N);
    
    K:=BaseRing(BaseRing(Q));
    C:=ZeroMatrix(K,2*g,2*g);

    //C is the sympectic matrix.
    for i:=1 to g do 
        C[i,g+i]:=1;
        C[g+i,i]:=-1; 
    end for;

    F := data`F; //Frobenius on the H1dR of the curve.
    F := Submatrix(data`F,1,1,2*g,2*g); 
    Finv := Transpose(F)^(-1);
    Aq := Transpose(F)+q*Finv;   // Eichler-Shimura -> Hecke operator
    prec_loss_bd := Valuation(Norm(Determinant(Finv)), p);

    AK := ZeroMatrix(K, 2*g, 2*g); ZK := AK; //initialising Hecke operators and nice correspondences. 
    bad_indices:=[[*[**],[**]*]:i in [1..g-1]]; //To record indices of the matrix where LLL does not work.
    Zs:=[]; As:=[];

    for i in [1..g-1] do
        A := Aq^i; // ith power of hecke operator
        Zmx := (2*g*A-Trace(A)*IdentityMatrix(K,2*g))*C^(-1);      // Zmx is a q-adic approximation of a nice correspondence Z, ie trace 0
        //The denominator created problems earlier, let's see if this works without this. 

        for j in [1..2*g] do
            for k in [1..2*g] do
        
                try
                    AK[j,k] := alg_approx_Qp(Qp!Rationals()!Aq[j,k],v);    // recognition of integer in Zp via LLL
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
                error "p-adic approximation of nice correspondence not exact.";  
            end if;  
            W:=ZeroMatrix(K,2*g, 2*g);
            W[1,g+1]:=Trace(ZK*C);
            W[g+1,1]:=-Trace(ZK*C);
            ZK:=2*ZK+W;
        end if;
        Append(~Zs,ZK);
        Append(~As,AK);
    end for;
    A:=As[1];

    return AK ,Zs, prec_loss_bd, bad_indices;
end intrinsic;


intrinsic HeckeCorrespondenceNF(Q::RngUPolElt[RngUPol], v::RngOrdIdl, N::RngIntElt :
                      useU:=false, useY:=false, basis0:=[], basis1:=[], basis2:=[], heights:=false)

                    ->AlgMatElt, SeqEnum, RngIntElt, SeqEnum
{Computes the Hecke Correspondence action given data for ColemanData}

    data:=ColemanData(Q, v, N);
    return HeckeCorrespondenceNF(data, 15);

end intrinsic;

