import numpy as np
import matplotlib.pyplot as plt
from math import sqrt
from scipy.sparse.linalg import spsolve

if __name__ == "__main__":
    mat = np.fromfile("matrix.csv", sep=" ")
    m = int(sqrt(mat.shape[0]))

    A = mat.reshape((m, m))
    plt.imshow(A, cmap="inferno")
    plt.show()
    # quit()

    # check if A is symmetric
    is_symmetric = np.allclose(A, A.T)

    # check if A is positive definite
    is_pos_def = np.all(np.linalg.eigvals(A) > 0)

    # check if A is singular
    is_singular = np.linalg.det(A) == 0

    print("A is symmetric: ", is_symmetric)
    print("A is positive definite: ", is_pos_def)
    print("A is singular: ", is_singular)

    # read b vector from vector_u.csv
    b = np.fromfile("rhs.csv", sep="\n")
    b = b.reshape((m, 1))

    # solve Ax = b
    p_inv = np.linalg.pinv(A)
    sol = np.dot(p_inv, b)

    # export solution to solution_u.csv
    np.savetxt("solution.csv", sol, delimiter=" ", fmt="%.6f")
