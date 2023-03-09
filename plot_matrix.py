import numpy as np
import matplotlib.pyplot as plt
from math import sqrt

if __name__ == '__main__':
    mat = np.fromfile('matrix_u.csv', sep=' ')
    m = int(sqrt(mat.shape[0]))
    
    A = mat.reshape((m, m))
    plt.imshow(A)
    plt.show()
    quit()
    
    # read b vector from vector_u.csv
    b = np.fromfile('vector_u.csv', sep='\n')
    b = b.reshape((m, 1))
    
    # solve Ax = b
    sol = np.linalg.solve(A, b)
    
    # export solution to solution_u.csv
    np.savetxt('solution_u.csv', sol, delimiter=' ', fmt='%.6f')
    
    