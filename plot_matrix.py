import numpy as np
import matplotlib.pyplot as plt
from math import sqrt

if __name__ == '__main__':
    mat = np.fromfile('matrix.csv', sep=' ')
    m = int(sqrt(mat.shape[0]))
    
    mat = mat.reshape((m, m))
    
    plt.imshow(mat, interpolation='nearest', cmap='gray')
    plt.show()

    