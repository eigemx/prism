add_library(prism-core
    # core headers
    core.h
    types.h
    print.h
    meta.h
    fvscheme.h
    equation.h
    field.h
    export.h

    # core sources
    field.cpp
    equation.cpp
    export.cpp
    
    # mesh module headers
    mesh/face.h
    mesh/cell.h
    mesh/boundary.h
    mesh/pmesh.h
    mesh/unv.h
    mesh/trie.h
    mesh/reorder.h

    mesh/face.cpp
    mesh/cell.cpp
    mesh/boundary.cpp
    mesh/pmesh.cpp 
    mesh/unv.cpp
    mesh/trie.cpp
    mesh/reorder.cpp

    # diffusion module
    diffusion/diffusion.h
    diffusion/diffusion.cpp

    # gradient module
    gradient/gradient.h
    gradient/gradient.cpp

    # solver module
    solver/solver.h
    solver/solver.cpp
)
