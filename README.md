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
- Visualization:
  - Solving over an arbitrary region of interest
  - Adaptive sampling & interpolation
- Applications:
  - Helmholtz decomposition


## Idea of the Paper
This paper introduces the _walk on spheres_ (WoS) algorithm for solving linear elliptic partial differential equations (PDEs) to the geometry processing (GP) community, and illustrates its potential in GP applications by presenting lots of examples. The two main analytical apparatus, the Kakutani's principle and the mean value property, allow one to write the solution $u(x)$ to the PDE in the form of a recursive integral formula
$$u(x) = \frac{1}{|\partial B(x)|}\int_{\partial B(x)} u(y)\,dy + \int_{B(x)} f(y) G(x,y)\, dy,$$ which can be easily approximated by a Monte Carlo estimate. Similarly, the gradient and higher-order derivatives can be approximated in the same manner with negligible overheads. 

The authors also introduce two variance reduction techniques for the WoS algorithm. In the first technique, control variates are used to reduce the variance of both the solution and the gradient. The most interesting part is that, the variance of the solution is controlled by the running estimate of the gradient, and the variance of the gradient is controlled by the running estimate of the solution. In other words, they reinforce each other. 

Although one of the advantages of the WoS algorithm is that one can evaluate the solution locally, it's often desired to solve a PDE on the whole domain and the WoS algorithm will behave poorly in this setting if without extra handling. Therefore, the author also introduced adaptive sampling and interpolation using moving least squares. The author proposed a simple adaptive sampling scheme that adds samples if the 1st order Taylor expansion at that point cannot approximate the solution well. 



## Build 
```
git clone --recursive http://github.com/shiinamiyuki/mcgp.git
cd mcgp && mkdir build
cd build && cmake -DCMAKE_BUILD_TYPE=Release ..
make -j `nproc`
```
We use a custom fork of libigl which uses [embree](https://www.embree.org/api.html) for high-performance BVH building and closest point query. Passing -DLIBIGL_WITH_EMBREE=ON to turn on this feature (OFF by default). It will take around 15 minutes for the first build.

## Sample Code
The `tutorials` folder contains many runnable examples. Make sure the working directory is `./build`. Note that some of the tutorials will write images to your working directory.


## Implementation
Our core algorithm is based on SDFs, since they allow one to query the distance between a point and the boundary in $\mathcal{O}(1)$ time, which leads to the most elegant and efficient way to do WoS. Moreover, all our library takes are an SDF that defines the domain, the corresponding boundary condition function, plus a source function defined in the PDE. And for triangle meshes inputs, we constructed a bounding volume hierarchy to speed up the signed distance query.

The effect of control variates is not substantial -- it reduces the variance by around 10% based on numerical experiments. Its implementation is fairly straightforward, except for a possible pitfall that is listed at the bottom of our report. 

When interpolating the solution using moving least squares, instead of fitting a polynomial for all sample points, we only fit the polynomial for k=16 nearest points. This is because MLS depends on sample distances which is unpredictable in adaptive sampling. In practice, this knn interpolation scheme works well.

For most of the demos, we visualize the solution by solving it on a cross-section of a 3d domain. The solution (which is a scalar field) is represented by a greyscale image. For the gradient and Helmholtz decomposition (see below), we directly plot the vector field.

We also implemented the 3-D Helmholtz decomposition as an application of MCGP in graphics. Note that they didn't use the most common definition of the Helmholtz decomposition
$$X = \nabla u + \nabla\times A,$$ instead, they assume that $u$ and $A$ satisfy homogeneous Dirichlet boundary condition, which results in
$$X = \nabla u + \nabla\times A + Y.$$ 
Note that $u$ and $A$ satisfies the scalar Poisson's equation $\nabla^2 u = \nabla\dot X$ and the vector Poisson's equation $\nabla^2 A = \nabla\times X$.

## Demos

#### 2-D Poisson's Equation
In this example, we solve a Poisson's equation defined inside a unit disk $D$:
$$-\nabla^2 u = 8\pi^2\cos(2\pi x)\sin(2\pi x) \quad \text{for } x\in D, $$
$$u(x,y) = cos(2\pi x)\sin(2\pi x)\quad \text{for } x\in \partial D.$$
The numerical results are as follows. The left and middle one are the WoS approximations when the number of steps is 1024 and 4096, respectively. The rightmost one is the ground truth.

<img src="https://github.com/shiinamiyuki/mcgp/blob/main/images/poi2d_1024.png" alt="1024"/><img src="https://github.com/shiinamiyuki/mcgp/blob/main/images/poi2d_4096.png" alt="4096"/><img src="https://github.com/shiinamiyuki/mcgp/blob/main/images/poi2d_groundtruth.png" alt="truth"/>

#### Helmholtz Decomposition
Note that the decomposition relies on the correct approximation of the gradient of the solutions, so it will demonstrate both the correctness of our gradient approximation, and one of the potential applications of the WoS algorithm. We also implemented a vector field visualizer to illustrate our result. In the following example, we set the input vector field to be $(x^2 \cos(y), xyz, e^{xy})$.

<img src="https://github.com/shiinamiyuki/mcgp/blob/main/images/curlfree1.jpg" alt="Curl free" width="300"/><img src="https://github.com/shiinamiyuki/mcgp/blob/main/images/curlfree2.jpg" alt="Curl free" width="304"/>

<img src="https://github.com/shiinamiyuki/mcgp/blob/main/images/divfree1.jpg" alt="Div free" width="320"/><img src="https://github.com/shiinamiyuki/mcgp/blob/main/images/divfree2.jpg" alt="Div free" width="299"/>

#### Adaptive Sampling
We compare adaptive sampling and uniform sampling in solving this Laplace equation with high-frequency boundary conditions. Both solutions are interpolated using MLS and are run for the same amount of walks (adaptive sampling actually uses fewer walks!). Note uniform sampling performs badly on the boundary, where the solution varies most.

<img src="https://github.com/shiinamiyuki/mcgp/blob/main/images/uniform.png" alt="Uniform Sampling" width="200"/><img src="https://github.com/shiinamiyuki/mcgp/blob/main/images/adaptive.png" alt="Adaptive Sampling" width="200"/><img src="https://github.com/shiinamiyuki/mcgp/blob/main/images/boundary condition.png" alt="Boundary Condition" width="200"/>

## Possible Improvements
Here is a list of features that we wish we could implement:
- Implement MCGP on CUDA
- Use a low-discrepancy sequence for random number generation (Quasi-Monte Carlo). However, QMC enforces some constraints on usage to prevent correlation between samples and such integration will make our library not as intuitive as it is now.
- By incorporating *A Simple and Robust Mutation Strategy for the Metropolis Light Transport Algorithm, Kelemen et al. 2002*, one can implement MCGP using Markov Chain Monte Carlo (MCMC). Though the effectiveness of this approach is unknown.

## Potential Problems in the Original Paper
The following discussions may contain mistakes, so please read it at your discretion.
- In section 2.3, the author defines Poisson's equation to be $\nabla^2 u = f$. However, all derivations of related formulas are based on the definition $-\nabla^2 u = f$.
- In section 2.6, the authors claim that the exterior problem can be solved by applying Russian roulette to conditionally terminate paths far from the domain. However, even for a simple Poisson's equation, this would not work as the Russian roulette itself generates an unbiased estimate and potentially an infinitely long walk. Thus it requires infinitely many walks to converge.
- In section 3.1, it says that "The WoS estimator just adds a single sample of the latter integral
for each step of the walk". But we think we should add a single sample for *each walk* instead of *each step of the walk*.
- In section 6.6, the author claimed that $A$ is the solution to the vector Poisson equation $\nabla^2 A =\nabla\times X$, but based on its reference, the formula should be $\nabla^2 A = -\nabla\times X$.
- In Appendix B, the author provides a way to importance sample 3D harmonic Green's function by first uniformly sampling $\hat{y}$ on a unit sphere, and then sample $r$ from a distribution proportional to $r^2sin(\theta)$. However, this turns out to be exactly the same as uniformly sampling the unit ball since importance sampling $r^2sin(\theta)$ yields in $r=u^{1/3}$, where $u \sim Uniform[0,1)$. Also, the reference [Devroye 1986, Section 9.4] for this technique cannot be found (there are only 8 sections).

## Appendix: What is Special about PDEs in GP
I was curious about what's special about PDE problems in the GP setting. Now I have some answers:
- Cauchy problems (e.g., diffusion curves)
- Treatments of different boundary representations (e.g., polygon soup, implicit surfaces, NURBS, meshes)
- People are not too picky about the accuracy of the solution (in the mathematical sense)
- Size of the mesh (i.e., the problem) is fixed (unless you do remeshing) and usually fairly large
