add_library(prism-core
    # core
    visitor.cpp
    types.h
    print.h
    visitor.h
    
    # mesh module
    mesh/face.cpp
    mesh/cell.cpp
    mesh/boundary.cpp
    mesh/pmesh.cpp 
    mesh/inspector.cpp
    mesh/unv.cpp
    
    mesh/face.h
    mesh/cell.h
    mesh/boundary.h
    mesh/pmesh.h
    mesh/inspector.h
    mesh/unv.h
)
