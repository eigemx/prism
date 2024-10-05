add_library(prism
    linear.cpp
    
    equation/noslip.cpp

    numerics/solver.cpp

    field/ifield.cpp
    field/tensor.cpp
    field/scalar.cpp
    field/velocity.cpp
    field/velocity_boundary.cpp
    field/pressure.cpp
    field/pressure_boundary.cpp
    
    export.cpp
    exceptions.cpp

    # mesh
    mesh/face.cpp
    mesh/cell.cpp
    mesh/boundary.cpp
    mesh/pmesh.cpp 
    mesh/unv.cpp
    mesh/trie.cpp
    mesh/reorder.cpp

    scheme/nonortho.cpp

    # gradient
    gradient/boundary.cpp
    gradient/gradient.cpp
    gradient/green_gauss.cpp
    gradient/least_squares.cpp

    # operations
    operations/rhie_chow.cpp
)
