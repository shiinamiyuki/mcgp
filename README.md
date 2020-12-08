# Monte Carlo Geometry Processing
An implementation of "Monte Carlo Geometry Processing: A Grid-Free Approach to PDE-Based Methods on Volumetric Domains"

## Features
- Randomized Poisson Solver based on the *walk on spheres* algorithm:
  - Supports 3-D interior/exterior Dirichlet problems and 2-D interior Dirichlet problems
  - Boundary can be represented by either meshes or signed distance fields (SDFs)
  - Cost-free gradient approximations
  - Allows for local evaluations of a solution
- Variance Reduction:
  - Control variates
  - Multiple importance sampling
- Visualization:
  - Solving over arbitrary region of interest
  - Adaptive Sampling & Interpolation
- Applications:
  - Helmholtz decomposition

## Potential Problems in the Original Paper
The following discussions may contain mistakes, so please read it at your own discretion.
- In section 2.3, the author defines Poisson's equation to be $\nabla u = f$. However, all derivations of related formulas are based on the definition $-\nabla u = f$.
- In section 3.1, it says that "The WoS estimator just adds a single sample of the latter integral
for each step of the walk". But we think we should add a single sample for *each walk* instead of *each step of the walk*.
