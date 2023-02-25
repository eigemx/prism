add_library(prism-core
    # core
    types.h
    print.h
    version.h
    
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
)
