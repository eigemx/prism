import numpy as np
import matplotlib.pyplot as plt
import sys

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python plot_matrix.py <file>")
        sys.exit(1)

    A = np.loadtxt(sys.argv[1], delimiter=",")  

    plt.imshow(A, cmap='inferno')
    plt.show()