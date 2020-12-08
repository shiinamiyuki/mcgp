# Monte Carlo Geometry Processing
An implementation of "Monte Carlo Geometry Processing: A Grid-Free Approach to PDE-Based Methods on Volumetric Domains"

## Features
- Randomized Poisson Solver based on the *walk on spheres* algorithm:
  - Supports 3-D interior/exterior Dirichlet problems and 2-D interior Dirichlet problems
  - Boundary can be represented by either meshes or signed distance fields (SDFs)
  - Free gradient approximations
  - Allows for local evaluations of a solution
- Variance Reduction:
  - Control variates
  - Multiple importance sampling
- Visualization:
  - Solving over arbitrary region of interest
  - Adaptive sampling & interpolation
- Applications:
  - Helmholtz decomposition


## Idea of the Paper
This paper introduces the _walk on spheres_ (WoS) algorithm for solving linear elliptic partial differential equations (PDEs) to the geometry processing (GP) community, and illustrates its potential in GP applications by presenting lots of examples. The two main analytical apparatus, the Kakutani's principle and the mean value property, allow one to write the solution $u(x)$ to the PDE in the form of a recursive integral formula
\[u(x) = \frac{1}{|\partial B(x)|}\int_{\partial B(x)} u(y)\,dy + \int_{B(x)} f(y) G(x,y)\, dy,\] which can be easily approximated by a Monte Carlo estimate. Similarly, the gradient and higher-order derivatives can be approximated by the same manner with negligible overheads. 

The authors also introduce two variance reduction techniques for the WoS algorithm. In the first technique, control variates are used to reduce the variance of both the solution and the gradient. The most interesting part is that, the variance of the solution is controlled by the running estimate of the gradient, and the variance of the gradient is controlled by the running estimate of the solution. In other words, they reinforce each other. __TODO__ importance sampling. 

Although one of the advantage of the WoS algorithm is that one can evaluate the solution locally, it's often desired to solve a PDE on the whole domain and the WoS algorithm will behave poorly in this setting if without extra handling. Therefore, the author also introduced __TODO__

## Implementation
To our best knowledge, as of today, this is the only open source MCGP with decent amount of features implemented. Although our library accepts triangle mesh inputs, the core algorithm is based on SDFs. Moreover, all it takes are a SDF that defines the domain, and the corresponding boundary condition function.


## Potential Problems in the Original Paper
The following discussions may contain mistakes, so please read it at your own discretion.
- In section 2.3, the author defines Poisson's equation to be $\nabla u = f$. However, all derivations of related formulas are based on the definition $-\nabla u = f$.
- In section 2.6, the authors claim that the exterior problem can be solved by applying Russian roulette to conditonally terminate path far from the domain. However, even for a simple poisson equation this would not work as the Russian roulette itself generates and unbiased estimate and potentially an infinitely long walk. Thus it requires infinitely many walks to converge.
- In section 3.1, it says that "The WoS estimator just adds a single sample of the latter integral
for each step of the walk". But we think we should add a single sample for *each walk* instead of *each step of the walk*.
- In Appendix B, the author provides a way to importance sample 3D harmonic Green's function by first uniformly sampling $\hat{y}$ on a unit sphere, and then sample $r$ from a distribution proportional to $r^2sin(\theta)$. However, this turns out to be exactly the same as uniformly sampling the unit ball since importance sampling $r^2sin(\theta)$ yields in $r=u^{1/3}$, where $u \sim Uniform[0,1)$. Also the reference [Devroye 1986, Section 9.4] for this technique cannot be found (there are only 8 sections).

## Appendix: What is Special about PDEs in GP
I was curious about what's special about PDE problems in the GP setting. Now I have some answers:
- Cauchy problems (e.g., diffusion curves)
- Treatments of different boundary representations (e.g., polygon soup, implicit surfaces, NURBS, meshes)
- People are not too picky about accuracy of the solution (in the mathematical sense)
- Size of the mesh (i.e., the problem) is fixed (unless you do remeshing) and usually fairly large
