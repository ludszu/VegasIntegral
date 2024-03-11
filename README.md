# Coulomb matrix elements (CME) for gated MoS2 quantum dots (QD) with Vegas algorithm

## Coulomb integrals

To solve an interacting problem of multiple electrons in a QD, scattering CME are needed. First, using a tight-binding (TB) model, single particle (SP) problem is solved, 
which yields SP solutions. TB problem is formulated in the basis of atomic orbitals $\phi(\boldsymbol{r})$, approximated by Slater-type orbitals (STO):

$\phi_{STO}^{lmn}(\boldsymbol{r},\zeta)=Y_l^m(\varphi,\theta)R_n(r,\zeta)$,

where $Y_l^m(\varphi,\theta)$ are spherical harmonics and $R_n(r,\zeta)$ are riadial functions, and $\zeta$ depends on the chemical element.

CME needed in this problem are 6-dimensional integrals of the form

$\langle ij|V_ee|kl\rangle=\int \phi_i^* (\boldsymbol{r}) \phi_j^* (\boldsymbol{r'}) \frac{e^2}{\kappa |\boldsymbol{r}-\boldsymbol{r'}|} \phi_k(\boldsymbol{r'})\phi_l(\boldsymbol{r}) d\boldsymbol{r} d\boldsymbol{r'}$

and $i,j,k,l$ label atomic orbitals.

## Vegas algorithm

This type of 6-dimensional integtal is difficult to calculate with a straight-forward implementation. For corect results, an adapted mesh is needed. Here,
a Gnu Scientific Library (GNU) function \textit{fgsl_monte_vegas} is used (with the fortran interface FGSL). This routine uses Monte Carlo method to compute the integral.
The integrand function is probled many times for random points from a defined mesh, based on importance sampling. This means the points are selected 
based on the probability distribution of the integrand - where the contributions to the integral are the largest. Vegas algorithm starts with a preliminary coarse mesh
and redefines it more accurately in final steps.

## References

[1] GNU Vegas Monte Carlo integration https://www.gnu.org/software/gsl/doc/html/montecarlo.html#vegas

[2] L. Szulakowska, M. Cygorek, M. Bieniek, and P. Hawrylak, “Valley- and spin-polarized broken-symmetry states of interacting electrons in gated Mo S 2 quantum dots,” 
Phys. Rev. B, vol. 102, no. 24, p. 245410, Dec. 2020, doi: 10.1103/PhysRevB.102.245410.

