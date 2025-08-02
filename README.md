# Prism

Prism is a C++20 library for solving partial differential equations (PDEs) using finite volume methods. It is designed to be simple and easy to use, with a focus on computational fluid dynamics (CFD). Prism has many capabilities such as:

- Handling unstructured polyhedral meshes (currently the only mesh reader implemented is for Ideas-UNV meshes).
- Mesh traversal with ease.
- Native support for scalar, vector and tensor fields.
- Support for non-orthogonal correction for diffusion scheme.
- Central difference, Upwind, Second-order upwind, QUICK and high-resolution (TVD) schemes.
- Support for explicit and implicit source terms.
- Support for implicit transiet schemes such as Backward Euler (first order) and Adam-Moulton (second order).
- Support for user defined boundary conditions (with many default boundary conditions available such as Fixed, FixedGradient, Symmetric, Outlet, ...).
- Exporting results to VTU format (right now supporting meshes with hexahedral, tetrahedral and pyramidal cells only).
- and much more...

Prism's main goal is to be simple, modular and easy to use. The following example shows how to solve steady state advection equation:

```cpp
    // solve for temperature advection: ∇.(UT) - ∇.(κ ∇T) = 0
    // where ρ is the density and U is the velocity vector.

    // first template parameter specifies the type for advective field and
    // the second for the transport field (temperature scalar field)
    // ∇.(UT)
    using div = scheme::convection::LinearUpwind<field::Velocity, field::Scalar>;

    // first template parameter specifies the type for diffusion coefficient and
    // the second for the transport field (temperature scalar field)
    // - ∇.(κ ∇T)
    using laplacian = scheme::diffusion::NonCorrected<field::UniformScalar, field::Scalar>;

    auto eqn = eqn::Transport(div(U, T),            //   ∇.(UT)
                              laplacian(kappa, T)   // - ∇.(κ ∇T) = 0
    );

    // solve
    auto solver = solver::BiCGSTAB<field::Scalar>();
    solver.solve(eqn);
```

The following example shows how pressure correction equation is implemented in prism, the complete code can be found in `examples\SIMPLESolver\main.cpp`

```cpp
    // here we create the type for diffusion (laplacian operator) with over-relaxed non-orthogonal corrector
    using laplacian_p = diffusion::Corrected<
                                        // diffusion coefficient type
                                        field::Tensor,
                                        // non-orthogonal correction type
                                        diffusion::nonortho::OverRelaxedCorrector,
                                        // transport field type
                                        field::Pressure
                                        >;
    // source term for divergence of velocity field with negative sign (sources are added to rhs by default)
    using div_U = source::Divergence<Sign::Negative, field::Velocity>;

    // assemble pressure correction equation
    auto pEqn = eqn::Transport<field::Pressure>(laplacian_p(D, P_prime), // - ∇.(D ∇P_prime)
                                                div_U(mDot)              // == - (∇.U)
    );
```

Prism is in early development stages and is not yet ready for production use. However, you can check out the examples folder for some simple usage examples.

## Examples

### 2D SIMPLE Algorithm - backward facing step

check `examples/SIMPLESolver/main.cpp` for complete implementation:

Momentum equation:

$$\nabla . (\rho UU) - \nabla. (\mu \nabla U) = - \nabla P$$

Pressure correction equation:

$$ -\nabla . (\mathcal{D} \nabla {P}^{\prime}) = - \nabla . U $$

<p align="center">
    <img alt="backwardFacingStep" src="https://github.com/eigemx/prism/blob/main/screenshots/backwardFacingStep.png?raw=true" width="90%">
</p>

### Heat Equation - 2D torus with fixed value boundary conditions

check `examples/laplacianSolver/main.cpp` for complete implementation

Heat equation:

$$ - \nabla . (\kappa . \nabla T) = 0 $$

<p align="center">
    <img alt="torus" src="https://github.com/eigemx/prism/blob/main/screenshots/torus.png?raw=true" width="90%">
</p>

### Poisson Equation - 2D slice with function form source

check `examples/poissonSolver/main.cpp` for complete implementation

Poisson equation:

$$ - \nabla . (\kappa . \nabla P) = \mathcal{f} $$

$$ f = 2 {\pi}^2 \sin(\pi x) \cos(\pi y) $$

<p align="center">
    <img alt="poissonSolverSlice" src="https://github.com/eigemx/prism/blob/main/screenshots/poisson.png?raw=true" width="90%">
</p>

## TODO

- Add documentation.
- Implement SIMPLEC, PRIME, PISO and PIMPLE solvers (currently only SIMPLE solver is implemented in example directory).
- Improve field::Vector and field::Tensor implementations using Eigen reference/map types.
- Implement k-epsilon and k-omega turbulence models.
- Implement vectorized operations overall in the codebase.

## Clone and build

> [!CAUTION]
>
> This is just a fun side project of mine. Use it at your own risk.

```bash
git clone --recurse-submodules -j8 https://github.com/eigemx/prism.git
cd prism
./vcpkg/bootstrap-vcpkg.sh -disableMetrics
./vcpkg/vcpkg install
cmake --preset=default
cd build
make
```

Once build is complete, check `bin` folder for the built applications, and also check `tests/cases` directory for many cases to play with and understand boundary field files in Prism.

```bash
./build/bin/SIMPLESolver tests/cases/pipeCoarse/mesh.unv
./build/bin/laplacianSolver tests/cases/torus/mesh.unv
./build/bin/advectionSolver tests/cases/duct/mesh.unv
./build/bind/poissonSolver tests/cases/poisson/mesh.unv
```

## How to contribute

If you want to contribute to Prism, please fork the repository and create a pull request. If you have any questions, feel free to open an issue.
