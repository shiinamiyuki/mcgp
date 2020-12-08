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
  - TODO
- Applications:
  - Helmholtz decomposition
