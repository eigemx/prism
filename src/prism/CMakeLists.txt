add_library(prism
    linear.cpp
    
    field/ifield.cpp
    field/boundary.cpp
    field/tensor.cpp
    field/scalar.cpp
    field/velocity.cpp
    field/pressure.cpp
    
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

    nonortho/nonortho.cpp

    # gradient
    gradient/green_gauss.cpp
    gradient/least_squares.cpp
    gradient/gradient.cpp
    gradient/boundary.cpp

    schemes/convection.cpp
    schemes/convection_boundary.cpp

    # schemes
    schemes/backward_euler.cpp

    # operations
    operations/rhie_chow.cpp
)
