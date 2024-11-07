freeze;


import "auxpolys.m": auxpolys, log;
import "singleintegrals.m": evalf0, is_bad, local_coord, set_point, tadicprec, teichmueller_pt, xy_coordinates;
import "misc.m": are_congruent, equivariant_splitting, eval_mat_R, eval_list, eval_Q, FindQpointQp, function_field, alg_approx_Qp, minprec, minval, minvalp, QpMatrix,QpSequence,QpPolynomial;
import "applications.m": Q_points, Qp_points, roots_with_prec, separate;
import "heights.m": E1_NF, E2_NF,  E1_tensor_E2_NF, expand_algebraic_function, frob_equiv_iso, height;
import "hensel.m": hensel_lift_n, two_variable_padic_system_solver;


// verbose flag determines how much information is printed during the computation.
declare verbose QCMod, 4;

intrinsic QCModAffine(Q::RngUPolElt[RngUPol], p::RngIntElt, known_points::SeqEnum, 
  correspondence_data::List : N := 15, prec := 2*N, basis0 := [], basis1 := [], 
  basis2 := [], data1:=0, data2:=0, base_point := 0, rho := 0)
                //      hecke_prime := 0, unit_root_splitting := false, eqsplit := 0,
                //      height_coeffs := [], rho := 0, use_log_basis := false, use_polys:=[])
  -> SeqEnum, BoolElt, SeqEnum, SeqEnum, List, SeqEnum, Rec, Rec
  {Run quadratic Chabauty given an equation of a plane affine curve satisfying Tuitman's conditions over a quadratic imaginary number field K,
   a prime p that is split in K and of good reduction (in the Balakrishnan-Tuitman sense), a sequence of known K-points and a list of nice correspondences
    constructed using powers of the Hecke operator at p. 
    
  Outputs a list of K-points mapping under all embedding to disks where Tuitman's Frobenius lift is defined.
  Also outputs additional information.}

// INPUT
//  * Q is a polynomial in (K[x])[y] with O_K-integral coefficients, where K is a imaginary
//    quadratic number field, defining a affine plane curve Y/K which is birational to an affine open
//    of a smooth projective model X and J = Jac(X). We requre the follwoing conditions to be checked:
//    * rk(J(K)) = 2(g(X))
//    * rank of Neron--Severi of J (over K) is greater thank 1. This is the same as the rank of the ring of
//      Rosati-fixed endomorphisms of J. 
//    These conditions are not checked!
//    
//  * p is a prime that splits in K such that for both primes above p X has good 
//    reduction in the sense of Balakrishnan-Tuitman and log is an isomorphism on 
//    the (J)(K\otimesQp)\to T\otimes Qp where T=H^0(X,\Omega^1) is the tangent space at the origin.
//  * known_points is a list of known points in Y(K).
//  * correspondence_data contains three entries. The first is the action of some correspondence on H1dR, which 
//    is currently assumed to be the Hecke correspondence T_p. The second entry is a list of
//    linearly independent correspondences which are nice in the sense of BDMTV, obtained from powers of Tp.
//    The third entry is a p-adic precision loss bound from the correspondence. This can be calculated by calling
//    the intrinsic HeckeCorrespondenceNF.

//  OPTIONAL PARAMETERS
//  * N is the p-adic precision used in the computations
//  * prec is the t-adic precision used for power series computations
//  * basis0 is a basis of the holomorphic differentials
//  * basis1 is a set of g independent meromorphic differentials s.t. basis0 and 
//    basis1 generate H^1_dR(X).
//  * Together with basis0, basis1, the sequence basis2 forms a basis of H^1_dR(U), 
//    where U is the affine patch we work on (bad disks removed).
//  * base_point is a pair [x,y] specifying a rational point on X to be used as a 
//    base point for Abel Jacobi. If base_point = 0, the first good affine rational 
//    point found is used.
//
//  OUTPUT:
//  ** good_affine_K_pts_xy, complete, sol_list, zero_list, double_zero_list, 
//     global_pts_local, bad_affine_K_pts_xy, data1, data2
//  **
//  where
//  * good_affine_K_pts_xy is a sequence of points (x,y) in Y(K) such that (x,y) is 
//    good in Y(K_i) for i=1,2, where K_i is the completion at the prime p_i above p.
//    under both embeddings of K into Qp (terminology as in Balakrishnan-Tuitman).
//  * The boolean complete is true iff the computation proves that 
//    good_affine_rat_pts_xy contains all points of this form.
//  * sol_list contains a list of pairs ((x1,y1), (x2,y2)) such that (x_i,
//    y_i) in Y(K_i) is a good point and the pair correponds to a common root of the
//    quadratic Chabauty functions constructed from correspondence_data.
//  * zero_list contains the common roots of the quadratic Chabauty
//    functions across all good residue polydisks (=pairs of good residue
//    disks, parametrized using a local coordinate x-x(P) for some P in
//    the disk).
//  * double_zero_list is a list of pairs [l,m] such that the residue
//    polydisk corresponding to the lth disk in X(K1) and the mth disk in
//    X(K2) contains a multiple root for all quadratic Chabauty functions.
//    These disks can be accessed using Qp_points(data_i).
//  * bad_affine_K_pts_xy is a list of points (x,y) in Y(K) such that (x,y)
//    is bad in Y(K_1) or Y(K_2).
//  * data1 is the Coleman data at p1 used in the algorithm.
//  * data2 is the Coleman data at p2 used in the algorithm.
//  
//


  // ==========================================================
  // ===                   CHECK INPUT                      ===
  // ==========================================================
  
  require IsPrime(p): "p must be prime number";


  // ==========================================================
  // ===                  INITIALIZATION                    ===
  // ==========================================================

  // Increase precision if it's too small compared to p-adic precision
  while  prec - 2*Log(p, prec) le N-5 do  // 5 comes from the constant c in the lower bound TODO: add ref.
    prec +:= 1; 
  end while;
    
  K := BaseRing(BaseRing(Q));
  d := Degree(K);
  require d eq 2: "K must be quadratic"; // TODO: Generalize
  v1 := data1`v;
  require Norm(v1) eq p: "p should split in K.";
  v2 := data2`v;
  K1,loc1 := Completion(K,v1);
  K2,loc2 := Completion(K,v2);
  OK := MaximalOrder(K);
  Qp := pAdicField(p,N);
  r,Delta,s := auxpolys(Q);


  // ==========================================================
  // ===               SYMPLECTIC BASIS                  ===
  // ==========================================================
  vprint QCMod, 2: " Computing a symplectic basis of H^1";
  
  h1basis, g, r, W0 := H1Basis(Q, v1); 
  assert W0 eq IdentityMatrix(Rationals(), Degree(Q)); // TODO: Generalize
  _,_,_,_ := H1Basis(Q,v2); 
  //The second iteration is just a check for Tuitman's conditions being satisfied at both places.

  if #basis0*#basis1 gt 0 then // Use the given basis
    h1basis := basis0 cat basis1;
  end if;
  vprintf QCMod, 3: " genus = %o.\n", g;
  if IsZero(rho) then 
    rho := g;       //If not given, we assume that the Picard number is equal to the genus
  end if;
  
  // h1basis is a basis of H^1 such that the first g elements span the regular
  // differentials. Construct a symplectic basis by changing the last g elements of 
  // h1basis.

  standard_sympl_mat := ZeroMatrix(K,2*g,2*g);
  for i in [1..g] do
    standard_sympl_mat[i,g+i] := 1; standard_sympl_mat[g+i,i] := -1;
  end for;

  vprint QCMod, 3: " Computing the cup product matrix";
  cpm_prec := 2*g;
  if assigned cpm then delete cpm; end if;
  repeat 
    try 
      // This takes too long, because it works over the splitting field at
      // infinity
      //cpm := CupProductMatrix(h1basis, Q, g, r, W0 : prec := cpm_prec);
      cpm := CupProductMatrix(h1basis, Q, g, r, W0 : prec := cpm_prec, split := false);
    catch e;
      cpm_prec +:= g;
      vprint QCMod, 4: "Precision %o too low for cup product computation. Increasing to %o", cpm_prec-g, cpm_prec;
    end try;
  until assigned cpm;
  vprint QCMod, 3: " Cup product matrix", cpm;
  if cpm ne standard_sympl_mat then 
    coefficients := SymplecticBasisH1(cpm); // Create coefficients of a symplectic basis in terms of h1basis
    new_complementary_basis := [&+[coefficients[i,j]*h1basis[j] : j in [1..2*g]] : i in [1..g]];
    sympl_basis := [h1basis[i] : i in [1..g]] cat new_complementary_basis;
    if not &and[&and[Valuation(c, v1) ge 0 : c in Coefficients(w[1])] : w in sympl_basis] then
      error "The computed symplectic basis is not integral. Please try a different prime or a different basis.";
    end if; 
    if not &and[&and[Valuation(c, v2) ge 0 : c in Coefficients(w[1])] : w in sympl_basis] then
      error "The computed symplectic basis is not integral. Please try a different prime or a different basis.";
    end if;
    vprintf QCMod, 3: " Symplectic basis of H^1:\n%o\n", sympl_basis;
    basis0 := [[sympl_basis[i,j] : j in [1..Degree(Q)]] : i in [1..g]]; // basis of regular differentials
    basis1 := [[sympl_basis[i,j] : j in [1..Degree(Q)]] : i in [g+1..2*g]];  // basis of complementary subspace
  end if;
  if Type(data1) eq RngIntElt then 
    // useY means we compute a basis of H1(Y), rather than H1(X) or H1(U)
    data1 := ColemanData(Q, v1, N : useY:=true,  basis0:=basis0, basis1:=basis1, basis2:=basis2);
  end  if;
  vprintf QCMod, 2: " Computed Coleman data at p=%o wrt symplectic basis to precision %o.\n", v1, N;

  if Type(data2) eq RngIntElt then 
    data2:= ColemanData(Q, v2, N : useY:=true,  basis0:=basis0, basis1:=basis1, basis2:=basis2);
  end if;  
  vprintf QCMod, 2: " Computed Coleman data at p=%o wrt symplectic basis to precision %o.\n", v2, N;

  prec := Max([prec, 40, tadicprec(data1, 1)]);
  prec := Max(tadicprec(data2, 1), tadicprec(data1, 1));
  S<t>    := LaurentSeriesRing(Qp,prec);
  S1<z1>  := LaurentSeriesRing(Qp,prec);
  S12<z2> := LaurentSeriesRing(S1,prec);

  function to_S12(f, i)
    // f is an element of S = Qp[t]
    // coerce f into S12 by replacing t by z_i
    if i eq 1 then return S12!f; end if;
    return z2^Valuation(f)*S12!Coefficients(f);
  end function;

  // ==========================================================
  // ===                    POINTS                       ===
  // ==========================================================
  search_bound := 1000;

  // small K-rational points as Tuitman points under the 2 embeddings
  // PointSearch in Magma is defined over the rationals, so enter points manually.
  Kpoints_1 := Q_points(data1,search_bound : known_points := known_points);
  Kpoints_2 := Q_points(data2,search_bound : known_points := known_points); 
  // TODO: Maybe get rid of this and compute the K-points using our code?

  Nfactor := 1.5; // Additional precision for root finding in Qp_points
  computed_Qp_points := false;
  repeat 
    try 
      Qppoints_1 := Qp_points(data1 : Nfactor := Nfactor); // One Q_p-point for every residue disk.
      computed_Qp_points := true;
    catch e; 
      Nfactor +:= 0.5;
    end try;
  until computed_Qp_points;

  computed_Qp_points := false;
  repeat 
    try 
      Qppoints_2 := Qp_points(data2 : Nfactor := Nfactor); // One Q_p-point for every residue disk.
      computed_Qp_points := true;
    catch e; 
      Nfactor +:= 0.5;
    end try;
  until computed_Qp_points;

  // Affine points where Frobenius lift isn't defined:
  bad_affine_Qppoints_1 := [P : P in Qppoints_1 | is_bad(P, data1) and not P`inf];
  bad_Kpoints_1 := [P : P in Kpoints_1 | is_bad(P, data1) and not P`inf];
  bad_affine_Qpindices_1 := [i : i in [1..#Qppoints_1] | is_bad(Qppoints_1[i], data1) and not Qppoints_1[i]`inf];
  
  bad_affine_Qppoints_2 := [P : P in Qppoints_2 | is_bad(P, data2) and not P`inf];
  bad_Kpoints_2 := [P : P in Kpoints_2 | is_bad(P, data2) and not P`inf];
  bad_affine_Qpindices_2 := [i : i in [1..#Qppoints_2] | is_bad(Qppoints_2[i], data2) and not Qppoints_2[i]`inf];
  // Affine points where Frobenius lift is defined:
  good_Kpoints_1 := [P : P in Kpoints_1 | not is_bad(P, data1) and not P`inf];
  good_K_Qp_indices_1 := [FindQpointQp(P,Qppoints_1) : P in good_Kpoints_1];
  numberofpoints_1 := #Qppoints_1;

  good_Kpoints_2 := [P : P in Kpoints_2 | not is_bad(P, data2) and not P`inf];
  good_K_Qp_indices_2 := [FindQpointQp(P,Qppoints_2) : P in good_Kpoints_2];
  numberofpoints_2 := #Qppoints_2;
  assert numberofpoints_1 eq numberofpoints_2;

  // Find xy-coordinates of the small affine K-points.
  // Use LLL for this.
  good_coordinates_1 := [xy_coordinates(P,data1) : P in good_Kpoints_1];
  good_affine_K_pts_xy := [[alg_approx_Qp(P[1], v1), alg_approx_Qp(P[2], v1)] : P in good_coordinates_1]; 
  bad_coordinates_1 := [xy_coordinates(P,data1) : P in bad_Kpoints_1];
  // TODO: This might not always work for very bad points. Not a problem
  // in our example.
  bad_affine_K_pts_xy := [[alg_approx_Qp(P[1], v1), alg_approx_Qp(P[2], v1)] : P in bad_coordinates_1]; 

  vprintf QCMod, 2: "\n Good affine K-points:\n%o\n", good_affine_K_pts_xy;
  vprintf QCMod, 2: "\n Bad affine K-points:\n%o\n", bad_affine_K_pts_xy;

  //if ISA(Type(base_point), RngIntElt) and IsZero(base_point) then  // No base point given, take the first possible one.
  assert ISA(Type(base_point), RngIntElt) and IsZero(base_point);  // No base point given, take the first possible one.
  global_base_point_index := 1;
  bK_1 := good_Kpoints_1[global_base_point_index];
  bK_2 := good_Kpoints_2[global_base_point_index]; // base point as Kpoint
  bK_xy := good_affine_K_pts_xy[global_base_point_index];  // xy-coordinates of base point
  //else 
  //  bQ := set_point(base_point[1], base_point[2], data1); // base point given
  //  bK_xy := base_point;
  //  global_base_point_index := Index(good_affine_K_pts_xy, base_point);
  //end if;
  local_base_point_index_1 := FindQpointQp(bK_1,Qppoints_1);
  local_base_point_index_2 := FindQpointQp(bK_2,Qppoints_2);       // Index of global base point in list of local points.

  FF<y>   := function_field(Q);
  x := BaseRing(FF).1;
  bpt   := CommonZeros([x-bK_xy[1], y-bK_xy[2]])[1]; // Base point as place on the function field
  vprintf QCMod, 2: "\n Using the base point %o.\n", bK_xy;
  good_affine_K_pts_xy_no_bpt := Remove(good_affine_K_pts_xy, global_base_point_index); 

  ks_1 := Exclude(good_K_Qp_indices_1, local_base_point_index_1);  // indices in Qppoints of good affine 
  ks_2 := Exclude(good_K_Qp_indices_2, local_base_point_index_2);  // K-points with base point removed
                                                             
  // compute Teichmueller representatives of good points
  teichpoints_1 := [**]; teichpoints_2 := [**];
  for i in [1..numberofpoints_1] do
    teichpoints_1[i] := is_bad(Qppoints_1[i],data1) select 0  else teichmueller_pt(Qppoints_1[i],data1); // No precision loss
  end for;
  for i in [1..numberofpoints_2] do
    teichpoints_2[i] := is_bad(Qppoints_2[i],data2) select 0  else teichmueller_pt(Qppoints_2[i],data2); // No precision loss
  end for;

  // ==========================================================
  // ===                  CORRESPONDENCES                 ===
  // ==========================================================

  vprint QCMod, 2: "\n Computing correspondences";

  // Want rho-1 independent `nice` correspondences.
  // Construct them using powers of Hecke operator
  //q := IsZero(hecke_prime) select p else hecke_prime;
  // Commented out, since we always take q=p.

  //if Type(correspondence_data) eq RngIntElt  then 
  // TODO: Make this work.
  //  correspondences, Tq, corr_loss := HeckeCorrespondenceQC(data1,q,N : basis0:=basis0,basis1:=basis1,use_polys:=use_polys);
  //else
  correspondences := correspondence_data[2]; 
  Tp:=correspondence_data[1];  
  corr_loss:=correspondence_data[3];
  //end if;

  Ncorr := N + Min(corr_loss, 0);
  // correspondences and Tq are provably correct to O(p^Ncorr), at least if q = p. We
  // represent them via rational approximations.
  Qpcorr := pAdicField(p, Ncorr);
  mat_space := KMatrixSpace(Qpcorr, 2*g, 2*g);
  vprintf QCMod, 3: "\nHecke operator at %o acting on H^1:\n%o\n", p, Tp;
  //if IsDiagonal(Tq) or Degree(CharacteristicPolynomial(Tq)) lt 2*g then
    //error "p-Adic approximation of Hecke operator does not generate the endomorphism algebra. Please pick a different prime. ";
  //end if;
  //if q ne p then
    //printf "\n WARNING: Using Hecke operator T_%o, but %o isn't our working prime %o. The result will not be provably correct.\n", q, q, p; 
  //end if;  

  correspondences_Qp1:=[QpMatrix(M,N,v1): M in correspondences];  

  correspondences_Qp2:=[QpMatrix(M,N,v2): M in correspondences];
  //if #use_polys eq 0 then
  // Check if Hecke operator generates. Need to do this using p-adic arithmetic.
  if Dimension(sub<mat_space | ChangeUniverse(correspondences_Qp1, mat_space)>) lt rho-1 then
    error "Powers of Hecke operator don't suffice to generate the space of nice correspondences";
  end if;
  if Dimension(sub<mat_space | ChangeUniverse(correspondences_Qp2, mat_space)>) lt rho-1 then
    error "Powers of Hecke operator don't suffice to generate the space of nice correspondences";
  end if;
  //end if;
  //end if;
    
  vprintf QCMod, 3: "\n Nice correspondences:\n%o\n\n", correspondences;
  number_of_correspondences := #correspondences;
  vprintf QCMod, 2: "\n number_of_correspondences:\n%o\n\n", number_of_correspondences;

  Tp_small := ExtractBlock(Tp,1,1,g,g);     // Hecke operator at p on H^0(X,Omega^1)
  char_poly_Tp := CharacteristicPolynomial(Tp_small);  
  Qp_ext := quo<PolynomialRing(Qp) | PolynomialRing(Rationals())!char_poly_Tp>;
  //The characteristic polynomial is defined over Z, so don't actually need any embedding.
  Salpha := quo<PolynomialRing(S) | PolynomialRing(Rationals())!char_poly_Tp>;
  S12alpha := quo<PolynomialRing(S12) | PolynomialRing(Rationals())!char_poly_Tp>;

  function to_S12alpha(f, i)
    // f is an element of Salpha = Qp[t]/(char_poly_Tp)
    // coerce it into S12alpha by replacing t by z
    eltseq_f12 := [to_S12(e, i) : e in Eltseq(f)];
    return &+[eltseq_f12[i]*S12alpha.1^(i-1) : i in [1..#eltseq_f12]];
  end function;
  
  // Compute an End0(J)-equivariant splitting of the Hodge filtration.
  
  //if IsZero(eqsplit) then
  /* TODO: Fix this
    if unit_root_splitting then 
      // Compute the unit root splitting 
      FQp := ChangeRing(ChangeRing(Submatrix(data1`F,1,1,2*g,2*g), Rationals()),Qp); // Frobenius over Qp
      char_poly_frob := CharacteristicPolynomial(FQp);
      fact := Factorisation(char_poly_frob);
      assert #fact ge 2;
      non_unit_root_char_poly := &*[t[1]^t[2] : t in fact | &and[Valuation(Coefficient(t[1],i)) gt 0 : i in [0..Degree(t[1])-1]]];
      assert Degree(non_unit_root_char_poly) eq g;
      Mp := EchelonForm(ChangeRing(Evaluate(non_unit_root_char_poly, FQp), pAdicField(p, N-2))); 
      assert Rank(Mp) eq g;
      // basis of the unit root subspace wrt symplectic basis
      W_wrt_simpl := Transpose(Submatrix(ChangeRing(Mp, Rationals()), 1,1,g,2*g));
      // the splitting matrix is the unique matrix leaving the holomorphic
      // differentials invariant and vanishing along the unit root subspace.
      W_lower := ExtractBlock(W_wrt_simpl, g+1, 1, g, g);
      W_upper_minus := [-Vector(RowSequence(W_wrt_simpl)[i]) : i in [1..g]];
      split := Transpose(Matrix(Solution(W_lower, W_upper_minus)));
      eqsplit := BlockMatrix(2, 1, [IdentityMatrix(Rationals(),g), split]);
    else 
    */
      eqsplit := equivariant_splitting(Tp);
    //end if; // unit_root_splitting
  //end if; // IsZero(eqsplit)

  //vprintf QCMod, 3: "eqsplit is %o", eqsplit;

  eqsplit1 := QpMatrix(eqsplit,Ncorr,v1);
  eqsplit2 := QpMatrix(eqsplit,Ncorr,v2);
  minvaleqsplit1 := minvalp(eqsplit, v1);
  minvaleqsplit2 := minvalp(eqsplit, v2);

  // Test equivariance of splitting 
  big_split := BlockMatrix(1,2,[eqsplit,ZeroMatrix(Rationals(),2*g,g)]);
  assert IsZero(big_split*Transpose(Tp) - Transpose(Tp)*big_split);     // Test equivariance
  vprintf QCMod, 3: "\n equivariant splitting:\n%o\n", eqsplit;
  
  

  //Sum of these quantiites below will need to account for both primes, we will fix them 
  //when they actually show up in Hodge/Frobenius/power series. 

  F1_lists := [* *]; // functions vanishing in rational points, one for each corresp
  F2_lists := [* *]; // functions vanishing in rational points, one for each corresp
  local_height_lists_1 := [* *]; // local height as power series 
  global_height_lists_1 := [* *]; // global height as power series 
  E1_E2_lists_1 := [* *]; // E1 tensor E2 as power series
  E1_lists_1 := [* *]; 
  E2_lists_1 := [* *]; 
  local_height_lists_2 := [* *]; // local height as power series 
  global_height_lists_2 := [* *]; // global height as power series 
  E1_E2_lists_2 := [* *]; // E1 tensor E2 as power series
  E1_lists_2 := [* *]; 
  E2_lists_2 := [* *]; 
  NE1E2Ps := Ncorr;  // Precision of E1 tensor E2 of auxiliary points
  Nhts := Ncorr; // Precision of local heights of auxiliary points
  Nexpansions := []; // Precision of power series expansion of local heights 
  c1s1 := []; c1s2 := [];
  valetas1 := []; valbetafils1 := [];
  maxdeggammafils1 := []; minvalgammafils1 := []; 
  valetas2 := []; valbetafils2 := [];
  maxdeggammafils2 := []; minvalgammafils2 := []; 
  dim := d^2*g;
  //if #height_coeffs eq 0 then  or not use_log_basis then 
    heights1 := [* *];    // local heights of auxiliary points. Different correspondences allowed (might cut down the # of necessary rational pts).
    heights2 := [* *];    // local heights of auxiliary points. Different correspondences allowed (might cut down the # of necessary rational pts).
    basis_found := false;
    super_space := VectorSpace(Qp, dim);
    //super_space := VectorSpace(Qp, g);
    E1_E2_subspace := sub<super_space | [Zero(super_space)]>;
    E1_E2_Ps := [ ]; // E1 tensor E2 of auxiliary points
  //end if;
  

  for l := 1 to number_of_correspondences do
  //for l := 1 to 1 do
    Z := correspondences[l];

    // ==========================================================
    // ===                     HODGE                       ===
    // ==========================================================
    
    vprintf QCMod: " Computing Hodge filtration for correspondence %o.\n", l;
    FF := function_field(Q); // function field of curve over K
    infplaces:=InfinitePlaces(FF);
    generic := false;
    if #infplaces eq 1 then
      finf := CharacteristicPolynomial(ResidueClassField(infplaces[1]).1);
      split := SplittingField(finf);
      finfsplit := ChangeRing(finf, split); 
      if SymmetricGroup(Degree(finf)) eq GaloisGroup(finf) then
        generic := true;
      end if;
    end if;



    if assigned betafil1 then delete betafil1; end if;
    hodge_prec := 5; 
    repeat
      try
        if generic then
          eta1,betafil1,gammafil1,hodge_loss1 := HodgeDataGeneric(data1, Z, bpt, hodge_prec);
        else
          eta1,betafil1,gammafil1,hodge_loss1 := HodgeDataSplittingField(Q,g,W0,data1`basis,Z,bpt : r:=r, prec:=hodge_prec);
        end if;
      catch e;
        e;
        hodge_prec +:= 5;
      end try;
    until assigned betafil1;

    repeat
      try
        if generic then
          eta2,betafil2,gammafil2,hodge_loss2 := HodgeDataGeneric(data2, Z, bpt, hodge_prec);
        else
          eta2,betafil2,gammafil2,hodge_loss2 := HodgeDataSplittingField(Q,g,W0,data2`basis,Z,bpt : r:=r, prec:=hodge_prec);
        end if;
      catch e;
        hodge_prec +:= 5;
      end try;
    until assigned betafil2;
    Nhodge := Ncorr + Min(Min(0, hodge_loss1),hodge_loss2); // Correct to Nhodge

    vprintf QCMod, 2: " eta =  %o,%o.\n", eta1,eta2; 
    vprintf QCMod, 2: " beta_fil  =  %o,%o.\n", betafil1,betafil2; 
    vprintf QCMod, 2: " gamma_fil =  %o,%o.\n\n", gammafil1,gammafil2; 

    Append(~valetas1, minvalp(eta1, v1));
    Append(~valbetafils1, minvalp(betafil1, v1));
    Append(~maxdeggammafils1, Max([Degree(a) : a in Eltseq(gammafil1)]));
    Append(~minvalgammafils1, 
        Min([Min([0] cat [Valuation(c, v1) : c in Coefficients(a)]) : a in Eltseq(gammafil1)]));

    Append(~valetas2, minvalp(eta2, v2));
    Append(~valbetafils2, minvalp(betafil2, v2));
    Append(~maxdeggammafils2, Max([Degree(a) : a in Eltseq(gammafil2)]));
    Append(~minvalgammafils2, 
        Min([Min([0] cat [Valuation(c, v2) : c in Coefficients(a)]) : a in Eltseq(gammafil2)]));

    // ==========================================================
    // ===                  FROBENIUS                      ===
    // ==========================================================

    b01 := teichmueller_pt(bK_1,data1);
    b02 := teichmueller_pt(bK_2,data2);
    vprintf QCMod: " Computing Frobenius structure for correspondence %o.\n", l;
    // xy-coordinates of Teichmueller points in discs of base point under
    // the 2 embeddings into Qp. Here we approximate elements of Qp using
    // integers. This is required for FrobeniusStructure, which
    // approximates p-adics using IntegerRing(p^n).
    b0pt1 := [Rationals()!c : c in xy_coordinates(b01, data1)];
    b0pt2 := [Rationals()!c : c in xy_coordinates(b02, data2)]; 
    Z1 := QpMatrix(Z, Nhodge, v1);
    Z1 := ChangeRing(Z1, Rationals());
    Z2 := QpMatrix(Z, Nhodge, v2);
    Z2 := ChangeRing(Z2, Rationals());
    G1, NG1 := FrobeniusStructure(data1,Z1,eta1,b0pt1 : N:=Nhodge); 
    G2, NG2 := FrobeniusStructure(data2,Z2,eta2,b0pt2 : N:=Nhodge); 
    G_list1 := [**]; G_list2 := [**]; 
    // evaluations of G at Teichmuellers of all good points (0 if bad)

    for i := 1 to numberofpoints_1 do
      if is_bad(Qppoints_1[i],data1) then
        G_list1[i]:=0;
      else
        P  := teichpoints_1[i]; // P is the Teichmueller point in this disk
        pt := [IntegerRing()!c : c in xy_coordinates(P, data1)]; 
        G_list1[i] := eval_mat_R(G1, pt, r, v1); // P is good, so no precision loss. 
      end if;
    end for;
    G_list2 := [**]; 
    // evaluations of G at Teichmuellers of all good points (0 if bad)
    for i := 1 to numberofpoints_2 do
      if is_bad(Qppoints_2[i],data2) then
        G_list2[i]:=0;
      else
        P  := teichpoints_2[i]; // P is the Teichmueller point in this disk
        pt := [IntegerRing()!c : c in xy_coordinates(P, data2)]; 
        G_list2[i] := eval_mat_R(G2, pt, r, v2); // P is good, so no precision loss. 
      end if;
    end for;
    Ncurrent := Min(Min(Nhodge, NG1),NG2);

    PhiAZb_to_b01, Nptb01 := ParallelTransport(bK_1,b01,Z1,eta1,data1 : 
                                                              prec:=prec,N:=Nhodge);
    for i := 1 to 2*g do
      PhiAZb_to_b01[2*g+2,i+1] := -PhiAZb_to_b01[2*g+2,i+1];
    end for;

    PhiAZb_to_b02, Nptb02 := ParallelTransport(bK_2,b02,Z2,eta2,data2 : 
                                                          prec:=prec,N:=Nhodge);
    for i := 1 to 2*g do
      PhiAZb_to_b02[2*g+2,i+1] := -PhiAZb_to_b02[2*g+2,i+1];
    end for;

    PhiAZb1 := [**]; // Frobenius on the phi-modules A_Z(b,P) (0 if P bad)
    PhiAZb2 := [**]; // Frobenius on the phi-modules A_Z(b,P) (0 if P bad)

    Ncurrent := Min(Min(Ncurrent, Nptb01),Nptb02);
    Nfrob_equiv_iso := Ncurrent;

    minvalPhiAZbs1 := [0 : i in [1..numberofpoints_1]];
    minvalPhiAZbs2 := [0 : i in [1..numberofpoints_2]];

    for i := 1 to numberofpoints_1 do

      if G_list1[i] eq 0 then 
        PhiAZb1[i] := 0;
      else 
        pti1, Npti1 := ParallelTransport(teichpoints_1[i],Qppoints_1[i],
                                                Z1,eta1,data1:prec:=prec,N:=Nhodge);
        isoi1, Nisoi1 := frob_equiv_iso(G_list1[i],data1,Ncurrent); 
        MNi1 := Npti1 lt Nisoi1 select Parent(pti1) else Parent(isoi1);
        PhiAZb1[i] := MNi1!(pti1*PhiAZb_to_b01*isoi1);
        Nfrob_equiv_iso1 := Min(Nfrob_equiv_iso, minprec(PhiAZb1[i]));
        minvalPhiAZbs1[i] := minval(PhiAZb1[i]);
      end if;
    end for;

    for i := 1 to numberofpoints_2 do

      if G_list2[i] eq 0 then 
        PhiAZb2[i] := 0;
      else 
        pti2, Npti2 := ParallelTransport(teichpoints_2[i],Qppoints_2[i],Z2,eta2,data2:prec:=prec,N:=Nhodge);
        isoi2, Nisoi2 := frob_equiv_iso(G_list2[i],data2,Ncurrent); 
        MNi2 := Npti2 lt Nisoi2 select Parent(pti2) else Parent(isoi2);
        PhiAZb2[i] := MNi2!(pti2*PhiAZb_to_b02*isoi2);
        Nfrob_equiv_iso2 := Min(Nfrob_equiv_iso, minprec(PhiAZb2[i]));
        minvalPhiAZbs2[i] := minval(PhiAZb2[i]);
      end if;
    end for;

    Ncurrent := Min(Nfrob_equiv_iso1, Nfrob_equiv_iso2);

    Append(~c1s1, Min(minvalPhiAZbs1));
    Append(~c1s2, Min(minvalPhiAZbs2));


    PhiAZb_to_z1 := [**]; 
    // Frobenius on the phi-modules A_Z(b,z) for z in residue disk of P (0 if P bad)
    for i := 1 to numberofpoints_1 do
      PhiAZb_to_z1[i] := G_list1[i] eq 0 select 0 else 
        ParallelTransportToZ(Qppoints_1[i],Z1,eta1,data1:prec:=prec,N:=Nhodge)*PhiAZb1[i]; 
    end for;


    PhiAZb_to_z2 := [**]; 
    // Frobenius on the phi-modules A_Z(b,z) for z in residue disk of P (0 if P bad)
    for i := 1 to numberofpoints_2 do
      PhiAZb_to_z2[i] := G_list2[i] eq 0 select 0 else 
        ParallelTransportToZ(Qppoints_2[i],Z2,eta2,data2 : 
                                                prec:=prec,N:=Nhodge)*PhiAZb2[i]; 
    end for;


    gammafil_listb_to_z1 := [* 0 : k in [1..numberofpoints_1] *]; 
    // evaluations of gammafil at local coordinates for all points 
    vprintf QCMod, 3: "Computing expansions of gamma_fil for embedding 1.\n";
    for i := 1 to numberofpoints_1 do
      if G_list1[i] ne 0 then
        gammafil_listb_to_z1[i] := expand_algebraic_function(Qppoints_1[i], 
                                            gammafil1, data1, Nhodge, prec);
      end if;
    end for;


    gammafil_listb_to_z2 := [* 0 : k in [1..numberofpoints_2] *]; 
    // evaluations of gammafil at local coordinates for all points 
    vprintf QCMod, 3: "Computing expansions of gamma_fil for embedding 2.\n";
    for i := 1 to numberofpoints_2 do
      if G_list2[i] ne 0 then
        gammafil_listb_to_z2[i] := expand_algebraic_function(Qppoints_2[i], 
                                            gammafil2, data2, Nhodge, prec);
      end if;
    end for;

    // ==========================================================
    // ===                     HEIGHTS                        ===
    // ==========================================================
    minvalchangebasis1 := 0; minvalchangebasis2 := 0;
    //if #height_coeffs eq 0 or not use_log_basis then // Compute heights of auxiliary points.

      if not basis_found then  // Find a point with non-zero E1 to write down a 
                               // basis of the Lie algebra. 
                               // To minimize precision loss, want small valuation of
                               // determinant of change of basis matrix.
        min_val_det_i := Ncurrent;
        for i := 1 to #good_affine_K_pts_xy_no_bpt do
          // E1(sigma1(Pi))
          Qpti1 := i lt global_base_point_index select good_Kpoints_1[i]
                              else good_Kpoints_1[i+1];

          pti1, Npti1 := ParallelTransport(Qppoints_1[ks_1[i]], Qpti1, 
                                                Z1,eta1,data1:prec:=prec,N:=Nhodge);

          MNi1 := Npti1 lt Precision(BaseRing(PhiAZb1[ks_1[i]])) select Parent(pti1) 
                                                      else Parent(PhiAZb1[ks_1[i]]);
          PhiP1 := MNi1!(pti1*PhiAZb1[ks_1[i]]);
          // TODO: Use E1_NF (not essential)
          E1Pi1 := Vector(BaseRing(PhiP1),g,[PhiP1[j+1,1] : j in [1..g]]);
          NE1Pi1 := Min([Ncurrent, minprec(E1Pi1)]);

          // E1(sigma2(Pi))
          Qpti2 := i lt global_base_point_index select good_Kpoints_2[i]
                              else good_Kpoints_2[i+1];
          pti2, Npti2 := ParallelTransport(Qppoints_2[ks_2[i]], Qpti2, 
                                            Z2,eta2,data2:prec:=prec,N:=Nhodge);

          MNi2 := Npti2 lt Precision(BaseRing(PhiAZb2[ks_2[i]])) select Parent(pti2) 
                                                      else Parent(PhiAZb2[ks_2[i]]);
          PhiP2 := MNi2!(pti2*PhiAZb2[ks_2[i]]);
          E1Pi2 := Vector(BaseRing(PhiP2),g,[PhiP2[j+1,1] : j in [1..g]]);
          NE1Pi2 := Min([Ncurrent, minprec(E1Pi2)]);
          NE1Pi := Min(NE1Pi1, NE1Pi2);

          //vector with E1(sigma1(Pi)),E1(sigma2(Pi))
          //E1Pi := Vector(BaseRing(PhiP2),d*g,Eltseq(E1Pi1) cat Eltseq(E1Pi2));
          
          basisH0star_i1 := [];
          basisH0star_i2 := [];
          Tp_small1 := ChangeRing(QpMatrix(Tp_small,Precision(BaseRing(E1Pi1)),v1),
                                    BaseRing(E1Pi1));
          Tp_small2 := ChangeRing(QpMatrix(Tp_small,Precision(BaseRing(E1Pi2)),v2),
                                    BaseRing(E1Pi2));
          for i := 0 to g-1 do
            // basis for H^0(Omega^1)^* generated by powers of iota(Tp) acting on E1(P)
            Append(~basisH0star_i1, Eltseq(E1Pi1*Tp_small1^i)); 
            Append(~basisH0star_i2, Eltseq(E1Pi2*Tp_small2^i)); 
          end for; 
          // If this breaks, it might be due to different base rings.
          val_det_i := Min([Valuation(Determinant(Matrix(M))) : 
                                    M in [basisH0star_i1, basisH0star_i2]]);
          if val_det_i lt min_val_det_i then
            // Better point found
            basis_found := true;
            min_val_det_i := val_det_i; min_i := i; 
            NH0star := NE1Pi;
            E1P1 := E1Pi1; E1P2 := E1Pi2; 
            basisH0star1 := basisH0star_i1;
            basisH0star2 := basisH0star_i2;
          end if;
          if IsZero(val_det_i) then break; end if;
        end for;
        if min_val_det_i ge Ncurrent then  // precision loss too high to obtain meaningful result.
          error "No good basis for H^0(Omega^1)^* generated by powers of iota(Tp) acting on E1(P) found";
        end if;
      end if; // basis_found


      changebasis1 := Matrix(basisH0star1)^(-1);
      changebasis2 := Matrix(basisH0star2)^(-1);
      minvalchangebasis1 := minval(changebasis2);
      minvalchangebasis2 := minval(changebasis2);
      changebases := [changebasis1, changebasis2];

      vprintf QCMod, 4: " Using point %o to generate.\n", 
                                        good_affine_K_pts_xy_no_bpt[min_i];

   // end if; 
  //end for;  // k := 1 to numberofpoints 

    // Compute heights of auxiliary points.
    // if #height_coeffs eq 0 or not use_log_basis then 
    // heights contains the list of heights of auxiliary points 
    // if #height_coeffs eq 0 then 
//      if Dimension(E1_E2_subspace) lt dim then  // add E1_E2(P) to known subspace until dimension is dim.
//            to fit the height pairing
      i := 1;
      repeat 
        // E1_tensor_E2(P1)
        Qpti1 := i lt global_base_point_index select good_Kpoints_1[i]
                    else good_Kpoints_1[i+1];
        pti1, Npti1 := ParallelTransport(Qppoints_1[ks_1[i]], Qpti1, Z1,eta1,data1 :                                          prec:=prec,N:=Nhodge);
        MNi1 := Npti1 lt Precision(BaseRing(PhiAZb1[ks_1[i]])) select Parent(pti1) 
                    else Parent(PhiAZb1[ks_1[i]]);
        Phii1 := MNi1!(pti1*PhiAZb1[ks_1[i]]);
        Ni1 := Min([Ncurrent, Precision(BaseRing(Phii1)), minprec(Phii1)]);
        Qpti2 := i lt global_base_point_index select good_Kpoints_2[i]
                    else good_Kpoints_2[i+1];

        pti2, Npti2 := ParallelTransport(Qppoints_2[ks_2[i]], Qpti2, Z2,eta2,data2 :                                          prec:=prec,N:=Nhodge);
        MNi2 := Npti2 lt Precision(BaseRing(PhiAZb2[ks_2[i]])) select Parent(pti2) 
                    else Parent(PhiAZb2[ks_2[i]]);
        Phii2 := MNi2!(pti2*PhiAZb2[ks_2[i]]);
        Ni2 := Min([Ncurrent, Precision(BaseRing(Phii2)), minprec(Phii2)]);
        Ni := Min(Ni1, Ni2);
        Qpi := pAdicField(p, Ni);
        Qpix := PolynomialRing(Qpi);
        Qp_ext := quo< Qpix | Qpix!PolynomialRing(Rationals())!char_poly_Tp>;
        Phiis := [Phii1, Phii2]; 
        betafils := [QpSequence(Eltseq(betafil1),Ni,v1), 
                      QpSequence(Eltseq(betafil2),Ni,v2)]; 
        
        E1_P := E1_NF(Phiis, changebases, Qp_ext); 
        E2_P := E2_NF(Phiis, betafils, changebases, Qp_ext); 
        E1_E2_P_Qpext := E1_tensor_E2_NF(E1_P, E2_P);
        //E1_E2_P_Qpext contains d^2 elements of Qp_ext
        E1_E2_P_Qp := &cat[Eltseq(E) : E in E1_E2_P_Qpext]; 
        //E1_E2_P_Qp contains d^2*g elements of Qp
        NE1E2P := Min(Ni,minprec(E1_E2_P_Qp));
        NLA := Integers()!Min([Precision(BaseRing(E1_E2_subspace)), NE1E2P]);
        // p^NLA is the precision for the linear algebra computation.
        new_super_space := VectorSpace(pAdicField(p, NLA), dim);
        old_basis := ChangeUniverse(Basis(E1_E2_subspace), new_super_space); 
        new_E1_E2_subspace := sub<new_super_space | old_basis cat 
                                              [new_super_space!E1_E2_P_Qp]>;
        //if Dimension(new_E1_E2_subspace) gt Dimension(E1_E2_subspace) then
        if Dimension(new_E1_E2_subspace) gt Dimension(E1_E2_subspace) or 
          Dimension(E1_E2_subspace) eq dim then  
          // We really only use first condition. The second one is there so that we 
          // can test whether the pairing we solve for is actually the height 
          // pairing. This is done by computing E1_E2 and the heights for all 
          // available points.
          if Dimension(new_E1_E2_subspace) gt Dimension(E1_E2_subspace) then //and i le 9 then
            vprintf QCMod, 3: " Using point %o at correspondence %o to fit the height pairing.\n", good_affine_K_pts_xy_no_bpt[i], l;
          else 
            vprintf QCMod, 3: " Not using point %o at correspondence %o to fit the height pairing; already have full space.\n", good_affine_K_pts_xy_no_bpt[i], l;
          end if;
          E1_E2_subspace := new_E1_E2_subspace; 
          
          vprintf QCMod, 4: " New dimension = %o.\n", Dimension(E1_E2_subspace);

          x1, y1 := Explode(xy_coordinates(Qpti1, data1));
          gammafilP_1 := eval_list(Eltseq(gammafil1), x1, y1, v1, Ni1);
          //vprintf QCMod, 4: " gammafil_P1=%o,\n", gammafilP_1;
          height_P_1 := height(Phii1,QpSequence(Eltseq(betafil1),Ni1,v1),gammafilP_1,eqsplit1,data1);
          NhtP1 := AbsolutePrecision(height_P_1); 
          
          Append(~heights1, height_P_1); // height of A_Z(b, P)
          vprintf QCMod, 4: " Added height for point %o and correspondence %o to heights1; new size of heights1 is %o\n",
                                good_affine_K_pts_xy_no_bpt[i], l, #heights1;
          x2, y2 := Explode(xy_coordinates(Qpti2, data2));
          gammafilP_2 := eval_list(Eltseq(gammafil2), x2, y2, v2, Ni2);
          //vprintf QCMod, 4: " gammafil_P2=%o,\n", gammafilP_2;
          height_P_2 := height(Phii2,QpSequence(Eltseq(betafil2),Ni2,v2),gammafilP_2,eqsplit2,data2);
          NhtP2 := AbsolutePrecision(height_P_2); 
          Append(~heights2, height_P_2); // height of A_Z(b, P)

          Append(~E1_E2_Ps, E1_E2_P_Qp);
          Nhts := Min([Nhts, NhtP1, NhtP2]);
          NE1E2Ps := Min(NE1E2Ps, NE1E2P);
        else
          vprintf QCMod, 3: " Not using point %o at correspondence %o to fit the height pairing because of dependence.\n", 
                              good_affine_K_pts_xy_no_bpt[i], l;

        end if;
        i +:= 1;
      //until Dimension(E1_E2_subspace) eq d*g or i gt #ks_1; 
      until i gt #ks_1; 
    //end if; // #height_coeffs eq 0


    vprintf QCMod, 3: "Computing expansions of local heights and of E1 and E2.\n";
    local_height_list_1 := [*0 : k in [1..numberofpoints_1]*];
    E1_list_1 := [*0 : k in [1..numberofpoints_1]*];
    E2_list_1 := [*0 : k in [1..numberofpoints_1]*];
    local_height_list_2 := [*0 : k in [1..numberofpoints_2]*];
    E1_list_2 := [*0 : k in [1..numberofpoints_2]*];
    E2_list_2 := [*0 : k in [1..numberofpoints_2]*];


    for k := 1 to numberofpoints_1 do
      if G_list1[k] ne 0 then

        local_height_list_1[k] := height(PhiAZb_to_z1[k],betafils[1],gammafil_listb_to_z1[k],eqsplit1,data1);
        local_height_list_2[k] := height(PhiAZb_to_z2[k],betafils[2],gammafil_listb_to_z2[k],eqsplit2,data2);
        Phiks := [PhiAZb_to_z1[k], PhiAZb_to_z2[k]];
        E1_k := E1_NF(Phiks, changebases, Salpha); 
        E2_k := E2_NF(Phiks, betafils, changebases, Salpha); 
        E1_list_1[k] := E1_k[1];
        E2_list_1[k] := E2_k[1];
        E1_list_2[k] := E1_k[2];
        E2_list_2[k] := E2_k[2];

//        if use_log_basis then 
//          E1_list_1[k] := [PhiAZb_to_z1[k,j,1] : j in [2..g+1]];
//          E2_list_1[k] := [PhiAZb_to_z1[k,2*g+2,g+1+j] - loc1(betafil1[j]) : j in [1..g]]; 
//        else 
          //E1_E2_list_1[k] := E1_tensor_E2(PhiAZb_to_z1[k],QpSequence(Eltseq(betafil1),N,v1),changebasis1,data1,Salpha);
//       end if;
//        if use_log_basis then 
//          E1_list_1[k] := [PhiAZb_to_z1[k,j,1] : j in [2..g+1]];
//          E2_list_1[k] := [PhiAZb_to_z1[k,2*g+2,g+1+j] - loc1(betafil1[j]) : j in [1..g]]; 
//        else 
          //E1_E2_list_2[k] := E1_tensor_E2(PhiAZb_to_z2[k],QpSequence(Eltseq(betafil2),N,v2),changebasis2,data2,Salpha);
//       end if;

      end if;
    end for;  // k := 1 to numberofpoints 
    
   
    Append(~local_height_lists_1, local_height_list_1);
    Append(~local_height_lists_2, local_height_list_2);
    //Append(~E1_E2_lists_1, E1_E2_list_1);
    Append(~E1_lists_1, E1_list_1); 
    // actually, all E1_list_1's should be the same, since E1 doesn't depend on Z
    // TODO: Add this as sanity check
    Append(~E2_lists_1, E2_list_1);
    // Append(~E1_E2_lists_2, E1_E2_list_2);
    Append(~E1_lists_2, E1_list_2);
    Append(~E2_lists_2, E2_list_2);
    Append(~Nexpansions, Ncurrent);

  end for; //for l to number_of_correspondences
           //

  //vprintf QCMod, 4: " E1_E2_Ps1=%o,\n", E1_E2_Ps1;
  //vprintf QCMod, 4: " E1_E2_Ps2=%o,\n", E1_E2_Ps2;

  //if #height_coeffs eq 0 and Dimension(E1_E2_subspace) lt dim then
  if Dimension(E1_E2_subspace) lt dim then
    error "Not enough K-points on the curve!"; // to span the symmetric square of the Mordell-Weil group";
  end if;

  //if #height_coeffs eq 0 then 
  // Write the height pairing as a linear combination of the basis of symmetric bilinear
  // pairings dual to the E1_E2-basis of the auxiliary points. 
  E1_E2_Ps_matrix := Matrix(pAdicField(p, NE1E2Ps), dim, dim, [E1_E2_Ps[i] : i in [1..dim]]); 
  mat := E1_E2_Ps_matrix^(-1) ;
  matprec := minprec(mat);
  Qpht := pAdicField(p, Min([matprec, NE1E2Ps, Nhts]));
  heights_vector1 := Matrix(Qpht, dim,1, [heights1[i] : i in [1..dim]]);
  heights_vector2 := Matrix(Qpht, dim,1, [heights2[i] : i in [1..dim]]);
  height_coeffs1 := ChangeRing(mat, Qpht)*heights_vector1;
  height_coeffs2 := ChangeRing(mat, Qpht)*heights_vector2;
  /*
  heights_cyc := [heights1[i]+heights2[i] : i in [1..#heights1]];
  heights_anti := [heights1[i]-heights2[i] : i in [1..#heights1]];
  heights_vector_cyc := Matrix(Qpht, dim,1, [heights1[i]+heights2[i] : i in [1..dim]]);
  heights_vector_anti := Matrix(Qpht, dim,1, [heights1[i]-heights2[i] : i in [1..dim]]);
  height_coeffs_cyc := ChangeRing(mat, Qpht)*heights_vector_cyc;
  height_coeffs_anti := ChangeRing(mat, Qpht)*heights_vector_anti;
  */
  //vprintf QCMod, 4: " height_coeffs1=\n%o,\n", height_coeffs1;
  //vprintf QCMod, 4: " height_coeffs2=\n%o,\n", height_coeffs2;
  //vprintf QCMod, 4: " height_coeffs_cyc=\n%o,\n", height_coeffs_cyc;
  //vprintf QCMod, 4: " height_coeffs_anti=\n%o,\n", height_coeffs_anti;
  // Precision of height_coeffs
  Nhtcoeffs := minprec(Eltseq(height_coeffs1) cat Eltseq(height_coeffs2)); 
  vprint QCMod, 2: "\n checking height_coeffs\n";
  for j := 1 to #heights1 do
    // Check that the global height pairing is computed correctly as a
    // bilinear pairing in terms of the E1-E2-basis. 
    //diffj1 := &+[Eltseq(height_coeffs_cyc)[i]*Eltseq(E1_E2_Ps[j])[i] : i in [1..dim]] - heights_cyc[j];
    //diffj2 := &+[Eltseq(height_coeffs_anti)[i]*Eltseq(E1_E2_Ps[j])[i] : i in [1..dim]] - heights_anti[j];
    diffj1 := &+[Eltseq(height_coeffs1)[i]*Eltseq(E1_E2_Ps[j])[i] : i in [1..dim]] - heights1[j];
    diffj2 := &+[Eltseq(height_coeffs2)[i]*Eltseq(E1_E2_Ps[j])[i] : i in [1..dim]] - heights2[j];
    assert Valuation(diffj1) ge Nhtcoeffs; // 
    assert Valuation(diffj2) ge Nhtcoeffs; 
    //vprintf QCMod, 4: " difference for j= %o and the first height is %o\n", j, diffj1;
    //vprintf QCMod, 4: " difference for j= %o and the second height is %o\n", j, diffj2;
  end for;
  if #heights1 gt dim then // Otherwise nothing was checked
    vprint QCMod, 2: "\n Height coefficients are correct!";
  end if;

//  end if;
                                                                             
  c3_1 := minval(Eltseq(height_coeffs1));
  c3_2 := minval(Eltseq(height_coeffs2));
  min_root_prec := N;  // smallest precision of roots of QC function


  // Find expansions of the quadratic Chabauty functions

  zero_list := [[] : i in [1..numberofpoints_1]];
  double_zero_list := [ ];
  sol_list  := [ ];
  vprintf QCMod, 3: "Computing expansions of quadratic Chabauty functions..\n";
  for k := 1 to number_of_correspondences do
  
    // Define some constants and functions for sanity checks on valuations
    // of the QC functions. These are used in the precision analysis.

    c2_1 := Min([0, valbetafils1[k], minvaleqsplit1, valbetafils1[k]+ minvaleqsplit1]); 
    i0_1:= 0;
    i0_threshold_1 := Min([valetas1[k], valbetafils1[k]/2, (minvalgammafils1[k]-c2_1)/2]);
    repeat 
      i0_1 +:= 1;
    until -Floor(log(p,i0_1)) le i0_threshold_1;

    function valF1(i) 
      // lower bound on valuations of coefficients in entries of F1_list
      assert i ge i0_1;
      valgammafili_1 := i le maxdeggammafils1[k] select minvalgammafils1[k] else 0;
      return -2*Floor(log(p,i)) + c1s1[k] + Min(c2_1, c1s1[k]+c3_1+2*minvalchangebasis1);
    end function;

    c2_2 := Min([0, valbetafils2[k], minvaleqsplit2, valbetafils2[k]+ minvaleqsplit2]); 
    i0_2:= 0;
    i0_threshold_2 := Min([valetas2[k], valbetafils2[k]/2, (minvalgammafils2[k]-c2_2)/2]);
    repeat 
      i0_2 +:= 1;
    until -Floor(log(p,i0_2)) le i0_threshold_2;
    "i0_1", i0_1; 
    "i0_2", i0_2; 

    function valF2(i) 
      // lower bound on valuations of coefficients in entries of F2_list
      assert i ge i0_2;
      valgammafili_2 := i le maxdeggammafils2[k] select minvalgammafils2[k] else 0;
      return -2*Floor(log(p,i)) + c1s2[k] + Min(c2_2,c1s2[k]+c3_2+2*minvalchangebasis2);
    end function;


    F1_list := [**]; F2_list := [**];
    // Go through residue polydisks D(l) x D(m), with respective
    // parameters z1 and z2
    for l := 1 to numberofpoints_1 do
      if G_list1[l] eq 0 then 
        F1_list[l] := 0; F2_list[l] := 0;
      else
        F1_list[l] := [**]; F2_list[l] := [**];
        for m := 1 to numberofpoints_2 do
          if G_list1[m] eq 0 then
            F1_list[l][m] := 0; F2_list[l][m] := 0;
          else
            // Find expansion of E1 tensor E2 in D(l)xD(m)
            E1s := [to_S12alpha(E1_lists_1[k][l],1), to_S12alpha(E1_lists_2[k][m],2)];
            E2s := [to_S12alpha(E2_lists_1[k][l],1), to_S12alpha(E2_lists_2[k][m],2)];
            E1_E2_S12alpha := E1_tensor_E2_NF(E1s, E2s);
            E1_E2_S12 := &cat[Eltseq(E) : E in E1_E2_S12alpha]; 
            global_height_1 := &+[height_coeffs1[j,1]*E1_E2_S12[j]:j in [1..dim]];
            global_height_2 := &+[height_coeffs2[j,1]*E1_E2_S12[j]:j in [1..dim]];
            hv1 := to_S12(local_height_lists_1[k,l], 1);
            hv2 := to_S12(local_height_lists_2[k,m], 2);
            f1 := global_height_1 - hv1; // 1st QC function in D(l)xD(m)
            f2 := global_height_2 - hv2; // 2nd QC function in D(l)xD(m)
                                         //
            // Check that lower bounds for valuations of QC functions are
            // not violated.
            for s := 0 to Degree(f1) do
              for t := 0 to Degree(Coefficient(f1,s)) do
                u := s+t;
                if u ge i0_1 then
                  assert Valuation(Coefficient(Coefficient(f1, s), t)) ge valF1(u);
                end if;
              end for;
            end for;
            for s := 0 to Degree(f2) do
              for t := 0 to Degree(Coefficient(f2,s)) do
                u := s+t;
                if u ge i0_2 then
                  assert Valuation(Coefficient(Coefficient(f2, s), t)) ge valF2(u);
                end if;
              end for;
            end for;

            F1_list[l][m] := f1;
            F2_list[l][m] := f2;
          end if;
        end for; // m := 1 to numberofpoints 
      end if;
    end for; // l := 1 to numberofpoints 
    Append(~F1_lists, F1_list);
    Append(~F2_lists, F2_list);
    // So Fi_lists[k][l][m] contains the quadratic Chabauty function for
    // - idele class character i 
    // - correspondence k
    // - in the residue polydisk D(l) x D(m) 

    Nend := Integers()!Min(Nexpansions[k], Nhtcoeffs); // Precision used for root finding 
    vprintf QCMod: " The quadratic Chabauty function for correspondence %o is correct to precision %o^%o.\n",  k, p, Nend;
    Qp_small   := pAdicField(p,Nend); 
    Qpt1<t1> := PowerSeriesRing(Qp_small,prec);
    Qpt12<t2> := PowerSeriesRing(Qpt1,prec);
    Qps12<s1,s2> := PolynomialRing(Qp_small, 2);

    function make_power_series(f)
      // f = f(z1,z2) is an element of S12.  
      // Coerce into Qpt12
      series := Qpt12!&+[(p*t2)^j*Evaluate(Coefficient(Qpt12!f,j), p*t1) 
                                          : j in [Valuation(f)..Degree(f)]];
      prect2 := AbsolutePrecision(series);
      prect1 := Min([AbsolutePrecision(c) : c in Coefficients(series)]);
      return series, Min(prect2, prect1);
    end function;

    function make_poly(f)
      // f = f(t1,t2) is an element of Qpt12.  
      // Truncate into Qps12
      coefs := Coefficients(f);
      poly := Qps12!0;
      for j in [1..#coefs] do
        coefsj := Coefficients(coefs[j]);
        if #coefsj gt 0 then
          poly +:= s2^(j-1) * s1^(Valuation(coefs[j])) * 
                            &+[s1^(i-1)*coefsj[i]:i in [1..#coefsj]];
        end if;
      end for;

      min_val := Min([minval(Coefficients(c)) : c in coefs | c ne 0]);
      return p^(-min_val)*s2^Valuation(f)*poly, min_val;
    end function;
    //
    // ==========================================================
    // ===                 FIND ZEROES                     ===
    // ==========================================================

    vprint QCMod, 2: " Find common zeroes of the quadratic Chabauty functions";
    for i := 1 to numberofpoints_1 do
      if G_list1[i] ne 0 then
        for m := 1 to numberofpoints_2 do
          if G_list1[m] ne 0 then

            g1, precg1 := make_power_series(F1_list[i,m]);
            g2, precg2 := make_power_series(F2_list[i,m]);
              
            bound_val_coeffs1 := valF1(precg1) + 2*precg1;
            bound_val_coeffs2 := valF2(precg2) + 2*precg2;
            if Min(bound_val_coeffs1, bound_val_coeffs2) lt N then  
              error "Root finding won't yield meaningful results. Lower p-adic precision or increase t-adic precision";
            end if;

            g1_poly, min_val1 := make_poly(g1);
            g2_poly, min_val2 := make_poly(g2);
            // Commented this out, since it takes longe and is less stable
            // than two_variable_padic_system_solver
            // roots, droots := hensel_lift_n([g1_poly,g2_poly], p, Nend-1);
            //if droots gt 0 then
            if k eq 1 then  // first correspondence: compute roots
              g1_poly_1 := g1_poly;
              g1_poly_2 := g2_poly;
              vprintf QCMod, 3: " Find zeroes of the first quadratic Chabauty function in the polydisk %o,%o\n", i,m;
              roots, droots := two_variable_padic_system_solver(g1_poly, g2_poly, p, 
                                                        Nend-1, Nend-1 :safety :=1);
              zero_list[i,m] := roots;
              if droots gt 0 then  // there are multiple roots in this polydisk
                Append(~double_zero_list, [i,m]);
              end if;
            else // second correspondence
              // Check if functions vanish at zeroes of first function.
              vprintf QCMod, 3: " Check which zeroes of the first quadratic Chabauty function in the polydisk %o,%o are also zeroes of the second quadratic Chabauty function \n", i,m;
              if #zero_list[i,m] gt 0 then

                if [k,i,m] notin double_zero_list  then  
                  // Hensel condition satisfied for first correspondence
                  // for all roots in this polydisk

                  for r in zero_list[i,m] do

                    val1 := Valuation(Evaluate(g1_poly, r));
                    val2 := Valuation(Evaluate(g2_poly, r));
                    if val1 ge Nend-3 and val2 ge Nend-3 then
                      // common root to prec Nend-3, as safety measure
                      // find affine local coordinates 
                      Pp1 := Qppoints_1[i];
                      xt1, bt1 := local_coord(Pp1,prec,data1);
                      Pp2 := Qppoints_2[m];
                      xt2, bt2 := local_coord(Pp2,prec,data2);
                      // Commented out the following lines, since W0=I_4 in our
                      // example. 
                      // TODO: Need to uncomment for other examples, and need to
                      // change the base ring of W0 from F(x) to Qp(x).
                      //W0invxt1 := Evaluate(W0^(-1), xt1);
                      //b_vector1 := Matrix(Parent(xt1), Degree(Q), 1, bt1);
                      //yt1 := &+[W0invxt1[2,j]*b_vector1[j,1] : j in [1..Degree(Q)]];
                      //W0invxt2 := Evaluate(W0^(-1), xt2);
                      //b_vector2 := Matrix(Parent(xt2), Degree(Q), 1, bt2);
                      //yt2 := &+[W0invxt2[2,j]*b_vector2[j,1] : j in [1..Degree(Q)]];
                      yt1 := bt1[2];
                      yt2 := bt2[2];
                      pt1 := [Qp_small!Evaluate(c, p*r[1]) : c in [xt1, yt1]];
                      pt2 := [Qp_small!Evaluate(c, p*r[2]) : c in [xt2, yt2]];
                      double_root := false;

                      Append(~sol_list, <[pt1,pt2], double_root>); 
                    end if; // val1 ge Nend-3 and val2 ge Nend-3 then
                  end for; // r in zero_list[1][i,m] 
                else
                  // multiple roots in this disc for correspondence 1. 
                  // find roots for correspondence 2. 
                  vprintf QCMod, 3: " The QC functions wrt the first correspondence have a root in polydisk %o,%o not satisfying the Hensel condition. Find roots of QC functions wrt the second correspondence. \n", i,m;
                  roots2, droots2 := two_variable_padic_system_solver(
                                g1_poly, g2_poly, p, Nend-1, Nend-1 :safety :=1);
                  assert droots2 eq 0; // Hensel condition satisfied
                  for r in roots2 do
                    // Now check if QC fns for first correspondence vanish
                    // at roots of QC fns for second correspondence
                    val1 := Valuation(Evaluate(g1_poly_1, r));
                    val2 := Valuation(Evaluate(g1_poly_2, r));
                    if val1 ge Nend-3 and val2 ge Nend-3 then
                      // common root to prec Nend-3, as safety measure
                      // find affine local coordinates 
                      Pp1 := Qppoints_1[i];
                      xt1, bt1 := local_coord(Pp1,prec,data1);
                      Pp2 := Qppoints_2[m];
                      xt2, bt2 := local_coord(Pp2,prec,data2);
                      // TODO: For other examples fix W0 issue, same as above.
                      yt1 := bt1[2];
                      yt2 := bt2[2];
                      pt1 := [Qp_small!Evaluate(c, p*r[1]) : c in [xt1, yt1]];
                      pt2 := [Qp_small!Evaluate(c, p*r[2]) : c in [xt2, yt2]];
                      double_root := false;
                      Append(~sol_list, <[pt1,pt2], double_root>); 
                    end if; // val1 ge Nend-3 and val2 ge Nend-3 then
                  end for; // r in roots2
                  
                end if; // [k,i,m] notin double_zero_list  then  
              end if; // #zero_list[i,m] gt 0 then
            end if; // k eq 1 
          end if; // G_list1[m] ne 0
        end for;  // m:=1 to numberofpoints 
      end if; // G_list1[i] ne 0
    end for;  // i:=1 to numberofpoints 
  end for;  // k := 1 to number_of_correspondences do
  // TODO: Fix this. What's minrootprec here?
  //vprintf QCMod, 2: " All roots of the quadratic Chabauty function(s) are correct to precision at least %o^%o.\n", p, min_root_prec;
  //end for;
  // ==========================================================
  // ===                 SANITY CHECK                       ===
  // ==========================================================

  /*
     Check that both QC functions vanish at the images of the known K-points
   */
  global_pts_local := [* *];
  vprintf QCMod, 2: "\n Check that QC-functionss vanish at known K-points.\n";
  for i := 1 to number_of_correspondences do
    F1_list := F1_lists[i];
    F2_list := F2_lists[i];
    for j in [1..#good_Kpoints_1] do
      point := good_affine_K_pts_xy[j];
      vprintf QCMod, 3: "\n Check that QC-functions vanish at known K-point %o for correspondence %o.\n  ", point, i; 
      P1 := good_Kpoints_1[j]; 
      ind1 := FindQpointQp(P1, Qppoints_1); 
      P2 := good_Kpoints_2[j]; 
      ind2 := FindQpointQp(P2, Qppoints_2); 
      Append(~global_pts_local, <P1,P2>); 
      if ind1 gt 0 and G_list1[ind1] ne 0 and ind2 gt 0 and G_list2[ind2] ne 0 then
        P1p := Qppoints_1[ind1];
        P2p := Qppoints_2[ind2];
        if (not is_bad(P1p,data1)) and (not is_bad(P2p,data2)) and (not P1`inf) then		
          coord1 := (P1`x - P1p`x);
          coord2 := (P2`x - P2p`x);
          f1 := F1_list[ind1,ind2];
          f2 := F2_list[ind1,ind2];
          eval1 := Qp_small!&+[(coord2)^l*Evaluate(Coefficient(f1,l), coord1) : l in [Valuation(f1)..Degree(f1)]];
          eval2 := Qp_small!&+[(coord2)^l*Evaluate(Coefficient(f2,l), coord1) : l in [Valuation(f2)..Degree(f2)]];
          assert Valuation(eval1) ge Nend; // First QC function vanishes
          assert Valuation(eval2) ge Nend; // Second QC function vanishes
          vprintf QCMod, 3: " The valuation of QCfun 1 for correspondence %o at point %o is %o.\n", i, j,  Valuation(eval1);
          vprintf QCMod, 3: " The valuation of QCfun 2 for correspondence %o at point %o is %o.\n", i, j,  Valuation(eval2);
        end if;
      end if;
    end for;
  end for; //  i := 1 to number_of_correspondences 
  vprintf QCMod, 2: "\n QC-functions vanish at all known K-points!\n";
  complete := #sol_list eq #good_affine_K_pts_xy and #double_zero_list eq 0;

  return good_affine_K_pts_xy, complete, sol_list, zero_list, double_zero_list, global_pts_local, bad_affine_K_pts_xy, data1, data2; //, F1_lists, F2_lists, Qppoints_1, Qppoints_2;

end intrinsic;

intrinsic HeckeOperatorGenerates(S::ModSym, p::RngIntElt)
  -> BoolElt
  {Check that the Hecke operator Tp generates the Hecke algebra}
  // S is a space of cusp forms
  Tp := HeckeOperator(S, p);
  return not IsDiagonal(Tp) and Degree(MinimalPolynomial(Tp)) eq Dimension(S) div 2;
end intrinsic;



