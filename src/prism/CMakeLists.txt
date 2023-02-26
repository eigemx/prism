add_library(prism-core
    # core
    types.h
    print.h
    meta.h
    visitor.h
    fvscheme.h
    equation.h

    visitor.cpp
    fvscheme.cpp
    equation.cpp
    
    # mesh module
    mesh/face.cpp
    mesh/cell.cpp
    mesh/boundary.cpp
    mesh/pmesh.cpp 
    mesh/unv.cpp
    
    mesh/face.h
    mesh/cell.h
    mesh/boundary.h
    mesh/pmesh.h
    mesh/unv.h

    diffusion/diffusion.h
    diffusion/diffusion.cpp

    gradient/gradient.h
    gradient/gradient.cpp
)
