AJ, 6/11: This is a work in progress.

This is Magma code to carry out quadratic Chabauty for nice (modular) curves X/K, where K 
is an imaginary quadratic field  and X is a curve of of genus > 1 over K satisfying
the following conditions :

* rank = 2*genus 
* enough rational points on the curve to solve for the height pairing. 
* p is a rational prime, which splits in K, and 
  * X has good reduction at all primes of K dividing p, and Assumption 1 in `Counting points 
    on curves using a map to P1`, Math. Comp. 2016 Tuitman is satisfied for all primes above p.
  * the closure of Jac(X)(K) in Jac(X)(K\otimes Qp) has finite index
  * the Hecke operator at p generates the Hecke algebra.

The theory is described in `RATIONAL POINTS ON THE NON-SPLIT CARTAN MODULAR CURVE OF
LEVEL 27 AND QUADRATIC CHABAUTY OVER NUMBER FIELDS' by  Jennifer S. Balakrishnan,
L. Alexander Betts, Daniel Hast, Aashraya Jha, J. Steffen M üller. This follows the papers
`Quadratic Chabauty for modular curves: Algorithms and Examples` and `Explicit Chabauty-Kim 
for the Split Cartan Modular Curve of Level 13` by J.S. Balakrishnan, N. Dogra, J.S. Müller,
J. Tuitman and J. Vonk.

Most of the code consists of wrappers around earlier version QCMod on Steffen 
Müller's repository. A lot of the code is written so that previous computations of p-adic heights
and power series can now be done over completions of number fields. 

Contributions were made by Jennifer Balakrishnan, Daniel Hast, Aashraya Jha, and Steffen Müller.

List of files:
-- qc_modular.m: Contains
   - QCModAffine: Main function, takes a plane affine curve (not necessarily 
      smooth) with integer coefficients, monic in y, and a prime p and outputs the rational points 
      in those disks where Tuitman's Frobenius lift is defined. We also require the action of a nice
      correspondence on H^1_dR(X). Also outputs additional information, such as additional p-adic solutions which don't look rational.
      Includes numerous optional arguments.
-- hodge.m: Computes the Hodge filtration using the algorithm described in section 4 of 
    `Explicit Chabauty-Kim for the Split Cartan Modular Curve of Level 13`
--hodge_generic.m: Wrapper of hodhe.m, which enables us to do a lot fo the linear algebra of 
    required to compute the Hodge filtration in a smaller ground field.
-- frobenius.m: Computes the Frobenius structure using the algorithm described in section 4 of 
    `Explicit Chabauty-Kim for the Split Cartan Modular Curve of Level 13`
-- hensel.m: Computes     
-- heights.m: Computes Nekovar heights as described in 
    `Explicit Chabauty-Kim for the Split Cartan Modular Curve of Level 13` and various
    related functions.
-- hecke_correspondence.m: Computes a Hecke operator using Eichler-Shimura and nice
    correspondences. 
-- symplectic_basis.m: Given a basis of H^1_dR of a smooth projective curve such that the
    first g elements generate regular differentials, computes the cup product and a
    symplectic basis with respect to the cup product.
-- misc.m: various functions, such as an implementation of rational reconstruction of p-adic
    numbers using LLL, rank computations using Kolyvagin-Logachev, equivariant splittings of
    the Hodge filtration of H^1 and coefficients mod p^N of p-adic points under Abel-Jacobi in
    terms of generators. Contains a function for 
-- applications.m, auxpolys.m, coho.m, coleman.m, froblift.m, reductions.m, singleintegrals.m: 
    Due to Jan Tuitman, computes Frobenius lifts and  Coleman integrals, based on 
      - Tuitman, `Counting points on curves using a map to P1`, Math. Comp. 2016
      - Tuitman, `Counting points on curves using a map to P1, II`, Finite Fields Appl, 2017
      - Balakrishnan-Tuitman, `Jennifer S. Balakrishnan and Jan Tuitman. Explicit Coleman
        integration for curves`, Math. Comp.
    with some minor modifications. 
-- Examples: Contains code to find Q(zeta_3) points of Xns+(27) and a smooth plane quartic X which appears as
   a quotient of Xns+(27). 

If you have questions or suggestions or if you find bugs, let me know.

Aashraya Jha, Boston University
aashjha@bu.edu

