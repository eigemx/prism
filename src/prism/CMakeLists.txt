add_library(prism
    linear.cpp
    
    equation/noslip.cpp

    field/ifield.cpp
    field/tensor.cpp
    field/scalar.cpp
    field/velocity.cpp
    field/velocity_boundary.cpp
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

    schemes/nonortho.cpp

    # gradient
    gradient/boundary.cpp

    # operations
    operations/rhie_chow.cpp
)
