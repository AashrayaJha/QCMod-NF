This is [Magma](https://magma.maths.usyd.edu.au/magma/) code to carry out quadratic Chabauty for nice (modular) curves X/K, where K 
is an imaginary quadratic field  and X is a curve of of genus > 1 over K satisfying
the following conditions:

* rank = 2*genus;
* enough rational points on the curve to solve for the height pairing;
* p is a rational prime, which splits in K, and 
  * X has good reduction at all primes of K dividing p, and Assumption 1 in Tuitman 2016 [\[6\]](#References) is satisfied for all primes above p;
  * the closure of $\mathrm{Jac}(X)(K)$ in $\mathrm{Jac}(X)(K \otimes \mathbf{Q}_p)$ has finite index;
  * the Hecke operator at p generates the Hecke algebra.

The theory is described in the paper "Rational points on the non-split Cartan modular curve of level 27 and quadratic Chabauty over number fields" by [Jennifer S. Balakrishnan](https://math.bu.edu/people/jbala/), [L. Alexander Betts](https://lalexanderbetts.net/), [Daniel Rayor Hast](https://github.com/HastD), [Aashraya Jha](https://sites.google.com/view/aashrayajha), and [J. Steffen Müller](https://www.rug.nl/staff/steffen.muller/) [\[1\]](#References).

This follows the papers
"Explicit Chabauty-Kim for the Split Cartan Modular Curve of Level 13"
and
"Quadratic Chabauty for modular curves: algorithms and examples"
by Balakrishnan, Dogra, Müller, Tuitman and Vonk [\[2, 3\]](#References).

The main example is a quotient $X_H'$ of the curve $X^+_{ns}(27)$, defined over
$F = \mathbf{Q}(\zeta_3)$. To compute the $F$-points on $X_H'$, run the code in the file
`Xns27/Fpoints_XH.m` from the top level folder.

Most of the code consists of wrappers built on top of https://github.com/steffenmueller/QCMod/. Much of the code is written so that previous computations of p-adic heights
and power series can now be done over completions of number fields. 

Contributions were made by Jennifer Balakrishnan, Daniel Hast, Aashraya Jha, and Steffen Müller.

# List of files
- `QCMod.spec`, `coleman.spec`: spec files for attaching relevant intrinsics for running `qc_modular` and aspects of Coleman integration required in quadratic Chabauty respectively. 
- `qc_modular.m`: Contains
   - `QCModAffine`: Main function, takes a plane affine curve (not necessarily 
    smooth) with integer coefficients, monic in y, a prime p, data about the action of nice correspondences
    on $H^1_{dR}(X)$ and outputs the rational points in those disks where Tuitman's Frobenius lift is defined. 
    Also outputs additional information, such as additional p-adic solutions which don't look rational.
    Includes numerous optional arguments.    
- `hodge.m`: Computes the Hodge filtration using the algorithm described in section 4 of [\[2\]](#References).
- `hodge_generic.m`: Wrapper of `hodge.m`, which enables us to do a lot of the linear algebra required to compute
    the Hodge filtration in a smaller ground field.
- `frobenius.m`: Computes the Frobenius structure using the algorithm described in section 4 of [\[2\]](#References), extended to work over number fields.
- `hensel.m`: Computes the common zeroes of a system of bivariate power series. 
- `heights.m`: Computes Nekovar heights as described in [\[2\]](#References) and various
    related functions.
- `hecke_correspondence.m`: Computes a Hecke operator using Eichler–Shimura and nice
    correspondences. 
- `symplectic_basis.m`: Given a basis of $H^1_{dR}$ of a smooth projective curve such that the
    first g elements generate regular differentials, computes the cup product and a
    symplectic basis with respect to the cup product.
- `misc.m`: various functions, such as an implementation of reconstruction of global numbers from 
    p-adic numbers using LLL, rank computations using Kolyvagin–Logachev, equivariant splittings of
    the Hodge filtration of $H^1$ and coefficients mod $p^N$ of p-adic points under Abel–Jacobi in
    terms of generators. 
- `hecke_correspondence.m`: Given the Frobenius data calculated via `ColemanData`, it will calculate the 
    Hecke correspondence at a place dividing p. Wrapper to do this directly with inputs the curve, the place
    of the number field and precision also provided.    
- `applications.m`, `auxpolys.m`, `coho.m`, `froblift.m`, `reductions.m`, `singleintegrals.m`: 
    Due to Jan Tuitman, computes Frobenius lifts and Coleman integrals, based on Tuitman 2016, Tuitman 2017, and Balakrishan–Tuitman 2020 [\[3, 6, 7\]](#References).
Modifications made to make functions work over number fields, primarily due to Daniel Hast and Steffen Müller.
- `string-replace.m`: Intrinsics to replace strings, also in [CHIMP](https://github.com/edgarcosta/CHIMP). Useful for printing Magma data. 

- `Xns27`: Contains code to compute $\mathbf{Q}(\zeta_3)$ points of $X^{+}_ {ns}(27)$ and a smooth plane quartic $X_H'$ which appears as a quotient of $X^{+}_{ns}(27)$ in Theorem 1.1 of our paper [\[1\]](#References).
    
    - We also compute the j-invariants of all the points computed using code provided to us by Jeremy Rouse in `test_j_invariants.m`. 

    - The `PrecomputedData` folder contains data required to run `QCModAffine` for the curve $X_H'$. It includes the action of the Hecke 
    operator on the de Rham cohomology of a prime above 13, and also some Coleman data required for computations in `QCModAffine`.

    - The `misc` folder contains some facts about the curve which were used along the way in quadratic Chabauty computations, but are not needed to run `QCModAffine`.
    - `Xns27/misc/Curve_endos.m` contains code to compute the geometric endomorphisms of the curve, and hence the Neron–Severi group as well.
    - The file `Xns27/misc/Discoid.m` contains a calculation which helps us determine the local heights at the prime above 3.
    -  The file Xns27/Misc/Group_H.sage contains a calculation comparing two different definitions of modular curve for the subgroup H, and show that in this case they are actually isomorphic.
    - The folder `XHmodel` contains code provided by Jeremy Rouse which was used to compute the model of $X_H'$ first presented by Rouse, Sutherland, and Zureick-Brown [\[5\]](#References).
   - `Curve_endos` contains code to compute the geometric endomorphisms of the curve, and hence the Neron–Severi group as well. 

# References

1. J. S. Balakrishnan, L. A. Betts, D. R. Hast, A. Jha, and J. S. Müller. "Rational points on the non-split Cartan modular curve of level 27 and quadratic Chabauty over number fields". In preparation.
2. J. S. Balakrishnan, N. Dogra, J. S. Müller, J. Tuitman, and J. Vonk. "Explicit Chabauty–Kim for the split Cartan modular curve of level 13". *Ann. of Math.* 189.3 (2019), pp. 885–944. DOI: [10.4007/annals.2019.189.3.6](https://doi.org/10.4007/annals.2019.189.3.6)
3. J. S. Balakrishnan, N. Dogra, J. S. Müller, J. Tuitman, and J. Vonk. "Quadratic Chabauty for modular curves: algorithms and examples”. *Compositio Math.* 159.6 (2023), pp. 1111–1152. DOI: [10.1112/S0010437X23007170](https://doi.org/10.1112/S0010437X23007170)
4. J. S. Balakrishnan and J. Tuitman. "Explicit Coleman integration for curves". *Math. Comp.* 89.326 (2020), pp. 2965–2984. DOI: [10.1090/mcom/3542](https://doi.org/10.1090/mcom/3542)
5. J. Rouse, A. V. Sutherland, and D. Zureick-Brown. "$\ell$-adic images of Galois for elliptic curves over $\mathbf{Q}$ (and an appendix with John Voight)". *Forum of Mathematics, Sigma*. Vol. 10. Cambridge University Press, 2022, pp. 1–63. DOI: [10.1017/fms.2022.38](https://doi.org/10.1017/fms.2022.38)
6. J. Tuitman. "Counting points on curves using a map to $\mathbf{P}^1$, I". *Math. Comp.* 85.298 (2016), pp. 961–981. DOI: [10.1090/mcom/2996](https://doi.org/10.1090/mcom/2996)
7. J. Tuitman. "Counting points on curves using a map to $\mathbf{P}^1$, II". *Finite Fields Appl.* 45 (2017), pp. 301–322. DOI: [10.1016/j.ffa.2016.12.008](https://doi.org/10.1016/j.ffa.2016.12.008)

If you have questions or suggestions or if you find bugs, let me know.

Aashraya Jha, Boston University
aashjha@bu.edu

