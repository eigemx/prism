add_library(prism-core
    # core
    core.h
    types.h
    print.h
    meta.h
    fvscheme.h
    equation.h

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
)
