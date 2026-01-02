# FEM_Gravity_Potential_Solver

The code solves the Poisson equation for gravitational potential:
$$\frac{d^{2}\Phi}{dx^{2}} = 4\pi G\rho(x)$$
where the gravitational constant is assumed to be $G=1$.

The problem is defined on the domain $x \in [0, 3]$ with the following Dirichlet boundary conditions:
$\Phi(0) = 5$ \
$\Phi(3) = 2$

The density function $\rho(x)$ is defined piecewise: \
$\rho(x) = -10$ for $x \in [0, 1]$ \
$\rho(x) = 1$ for $x \in (1, 2]$  \
$\rho(x) = -10$ for $x \in (2, 3]$


### Numerical Method

**Weak Formulation:** The problem is transformed into a weak form $B(w, v) = L(v)$ using a linear shift function $\tilde{\Phi}(x) = -x + 5$ to handle non-zero boundary conditions.

**Basis Functions:** Linear "hat" functions are used for the finite element space.

**Numerical Integration:** Integrals are calculated using the Gauss-Legendre quadrature (2-point rule) to ensure high precision.

### FeaturesDynamic Mesh: 
- The number of elements $n$ is passed as a command-line parameter. 
- Visualization: Automated plot generation of the potential distribution $\Phi(x)$ using Matplotlib.
- Detailed Output: Results for each node are displayed in the terminal.