# Prism
> [!CAUTION]
>
> This is just a fun side project of mine. It is probably broken. Use it at your own risk.

Prism is a C++20 library for solving partial differential equations (PDEs) using finite volume methods. It is designed to be simple and easy to use, with a focus on computational fluid dynamics (CFD). Prism has many capabilities such as:
- Handling unstructured polyhedral meshes (currently only Ideas-UNV meshes are supported).
- Mesh traversal with ease.
- Native support for scalar, vector and tensor fields.
- Support for non-orthogonal correction for diffusion scheme.
- Central difference, Upwind, Second-order upwind, and QUICK convection schemes.
- Support for explicit and implicit source terms.
- Support for user defined boundary conditions (with many default boundary conditions available such as Fixed, FixedGradient, Symmetric, Outlet, ...).
- Support for implicit and explicit under relaxation.
- Exporting results to VTU format (right now supporting meshes with hexahedral, tetrahedral and pyramidal cells only).
- and much more...

Prism's main goal is to be simple, modular and easy to use. The following example shows how to solve steady state advection equation:
```cpp
    // solve for temperature advection: ∇.(ρUT) - ∇.(κ ∇T) = S
    // where ρ is the density and U is the velocity vector, and S is an arbitraty constant source
    using div = scheme::convection::SecondOrderUpwind<field::UniformScalar, field::Scalar>;
    using laplacian = scheme::diffusion::NonCorrected<field::UniformScalar, field::Scalar>;

    auto eqn = eqn::Transport(div(rho, U, T),     //   ∇.(ρUT)
                              laplacian(kappa, T) // - ∇.(κ ∇T) = 0
    );

    // solve
    auto solver = solver::BiCGSTAB<field::Scalar>();
    solver.solve(eqn, 100, 1e-5, 1);
```

Prism is in early development stages and is not yet ready for production use. However, you can check out the examples folder for some simple usage examples.

## TODO
- Add documentation.
- Implement SIMPLE & PISO solvers (basic SIMPLE solver is in example directory, but not yet functional).
- Implement backward Euler scheme for transient problems.
- Improve field::Vector and field::Tensor implementations using Eigen reference types.
- Implement k-epsilon and k-omega turbulence models.

## Clone and build
```bash
git clone --recurse-submodules -j8 https://github.com/eigemx/prism.git
cd prism
./vcpkg/bootstrap-vcpkg.sh -disableMetrics
./vcpkg/vcpkg install
cmake --preset=default
cd build
make
```

## How to contribute
If you want to contribute to Prism, please fork the repository and create a pull request. If you have any questions, feel free to open an issue.

