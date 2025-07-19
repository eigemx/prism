add_library(prism
    linear.cpp
    
    equation/noslip.cpp
    equation/symmetry.cpp
    equation/outlet.cpp

    numerics/solver.cpp

    field/ifield.cpp
    field/tensor.cpp
    field/scalar.cpp
    field/vector.cpp
    field/pressure.cpp
    field/pressure_boundary.cpp
    
    export.cpp
    exceptions.cpp

    mesh/face.cpp
    mesh/cell.cpp
    mesh/boundary.cpp
    mesh/pmesh.cpp 
    mesh/unv.cpp
    mesh/trie.cpp
    mesh/reorder.cpp
    mesh/utilities.cpp

    scheme/nonortho.cpp

    gradient/boundary.cpp
    gradient/gradient.cpp
    gradient/green_gauss.cpp
    gradient/least_squares.cpp

    operations/rhie_chow.cpp
)
