This is Magma code to carry out quadratic Chabauty for nice (modular) curves X/K, where K 
is an imaginary quadratic field  and X is a curve of of genus > 1 over K satisfying
the following conditions :

* rank = 2*genus 
* enough rational points on the curve to solve for the height pairing. 
* p is a rational prime, which splits in K, and 
  * X has good reduction at all primes of K dividing p, and Assumption 1 in `Counting points 
    on curves using a map to P1`, Math. Comp. 2016 Tuitman is satisfied for all primes above p.
  * the closure of Jac(X)(K) in Jac(X)(K \otimes Qp) has finite index
  * the Hecke operator at p generates the Hecke algebra.

The theory is described in `RATIONAL POINTS ON THE NON-SPLIT CARTAN MODULAR CURVE OF
LEVEL 27 AND QUADRATIC CHABAUTY OVER NUMBER FIELDS' by  Jennifer S. Balakrishnan,
L. Alexander Betts, Daniel Hast, Aashraya Jha, J. Steffen Müller. This follows the papers
`Quadratic Chabauty for modular curves: Algorithms and Examples` and `Explicit Chabauty-Kim 
for the Split Cartan Modular Curve of Level 13` by J.S. Balakrishnan, N. Dogra, J.S. Müller,
J. Tuitman and J. Vonk.

The main example is a quotient X_H' of the curve X^+_{ns}(27), defined over
F = Q(zeta_3). To compute the F-points on X_H', run the code in the file
Xns27/Fpoints_XH.m from the top level folder.

Most of the code consists of wrappers built on top of https://github.com/steffenmueller/QCMod/. Much of the code is written so that previous computations of p-adic heights
and power series can now be done over completions of number fields. 

Contributions were made by Jennifer Balakrishnan, Daniel Hast, Aashraya Jha, and Steffen Müller.

List of files:
-- QCMod.spec, coleman.spec: spec files for attaching relevant intrinsics for running qc_modular and aspects
   of coleman integration required in quadratic Chabauty respectively. 
-- qc_modular.m: Contains
   - QCModAffine: Main function, takes a plane affine curve (not necessarily 
    smooth) with integer coefficients, monic in y, a prime p, data about the action of nice correspondences
    on H^1_dR(X) and outputs the rational points in those disks where Tuitman's Frobenius lift is defined. 
    Also outputs additional information, such as additional p-adic solutions which don't look rational.
    Includes numerous optional arguments.    
-- hodge.m: Computes the Hodge filtration using the algorithm described in section 4 of 
    `Explicit Chabauty-Kim for the Split Cartan Modular Curve of Level 13`
-- hodge_generic.m: Wrapper of hodge.m, which enables us to do a lot of the linear algebra required to compute
    the Hodge filtration in a smaller ground field.
-- frobenius.m: Computes the Frobenius structure using the algorithm described in section 4 of 
    `Explicit Chabauty-Kim for the Split Cartan Modular Curve of Level 13` extended to work over number fields
-- hensel.m: Computes the common zeroes of a system of bivariate power series. 
-- heights.m: Computes Nekovar heights as described in 
    `Explicit Chabauty-Kim for the Split Cartan Modular Curve of Level 13` and various
    related functions.
-- hecke_correspondence.m: Computes a Hecke operator using Eichler-Shimura and nice
    correspondences. 
-- symplectic_basis.m: Given a basis of H^1_dR of a smooth projective curve such that the
    first g elements generate regular differentials, computes the cup product and a
    symplectic basis with respect to the cup product.
-- misc.m: various functions, such as an implementation of reconstruction of global numbers from 
    p-adic numbers using LLL, rank computations using Kolyvagin-Logachev, equivariant splittings of
    the Hodge filtration of H^1 and coefficients mod p^N of p-adic points under Abel-Jacobi in
    terms of generators. 
-- hecke_correspondence.m: Given the Frobenius data calculated via ColemanData, it will calculate the 
    Hecke_correspondence at a place dividing p. Wrapper to do this directly with inputs the curve, the place
    of the number field and precision also provided.    
-- applications.m, auxpolys.m, coho.m,, froblift.m, reductions.m, singleintegrals.m: 
    Due to Jan Tuitman, computes Frobenius lifts and  Coleman integrals, based on 
      - Tuitman, `Counting points on curves using a map to P1`, Math. Comp. 2016
      - Tuitman, `Counting points on curves using a map to P1, II`, Finite Fields Appl, 2017
      - Balakrishnan-Tuitman, `Jennifer S. Balakrishnan and Jan Tuitman. Explicit Coleman
        integration for curves`, Math. Comp. 
    Modifications made to make functions work over number fields, primarily due to Daniel Hast and Steffen Müller.
-- string-replace.m: Intrinsics to replace strings, also in CHIMP. Useful for printing magma data. 

-- Xns27: Contains code to compute Q(zeta_3) points of X^+_{ns}(27) and a smooth plane quartic X_H' which appears as a 
    quotient of X^+_{ns}(27) in Theorem 1.1 of
    `RATIONAL POINTS ON THE NON-SPLIT CARTAN MODULAR CURVE OF LEVEL 27 AND QUADRATIC CHABAUTY OVER NUMBER FIELDS' by  
    Jennifer S. Balakrishnan, L. Alexander Betts, Daniel Hast, Aashraya Jha, J. Steffen M üller. 
    
    We also compute the j-invariants of all the points computed using code provided to us by Jeremy Rouse in test_j_invariants.m . 

    The PrecomputedData folder contains data required to run QCModAffine for the curve XH. It includes the action of the Hecke 
    operator on the de-Rham cohomology of a prime above 13, and also some Coleman data required for computations in QCModAffine.

    The misc folder contains some facts about the curve which were used along the way in QC computations, but are not needed to run QCModAffine. 
    Xns27/misc/Curve_endos.m contains code to compute the geometric endomorphisms of the curve, and hence the Neron--Severi group as well. The file 
    Xns27/misc/Discoid.m contains a calculation which helps us determine the local heights at the prime above 3. The file Xns27/Misc/Group_H.sage contains
    a calculation comparing two different definitions of modular curve for the subgroup H, and show that in this case they are actually isomorphic.
    The folder XHmodel contains code provided by Jeremy Rouse which was used to compute the model of XH first presented in RSZB.

-- data: Most of this is data for the smooth plane quartic X_H'. Contains  coleman_good_patch.m and 
   Hecke_good_patch_400 which are used in running QCModAffine for X_H'. Contains some files which we used to collect data,
   and the misc folder contains some facts about the curve which were used along the way in QC computations, but are not needed to run QCModAffine. 
   Curve_endos contains code to compute the geometric endomorphisms of the curve, and hence the Neron--Severi group as well. 
   The folder XH model contains code provided by Jeremy Rouse which was used to compute the model of XH first presented in RSZB.

--tests: Some tests while extending code to run for number fields. Used examples in https://github.com/steffenmueller/QCMod/. 

If you have questions or suggestions or if you find bugs, let me know.

Aashraya Jha, Boston University
aashjha@bu.edu

