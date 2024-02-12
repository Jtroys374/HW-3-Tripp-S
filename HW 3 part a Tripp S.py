def cholesky_decomposition(A):
    n = len(A)
    L = [[0.0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i + 1):
            if i == j:
                sum_val = sum(L[i][k] ** 2 for k in range(j))
                L[i][j] = (A[i][i] - sum_val) ** 0.5
            else:
                sum_val = sum(L[i][k] * L[j][k] for k in range(j))
                L[i][j] = (A[i][j] - sum_val) / L[j][j]

    return L


def forward_substitution(L, b):
    n = len(L)
    y = [0.0] * n

    for i in range(n):
        y[i] = (b[i] - sum(L[i][j] * y[j] for j in range(i))) / L[i][i]

    return y


def backward_substitution(L_transpose, y):
    n = len(L_transpose)
    x = [0.0] * n

    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - sum(L_transpose[i][j] * x[j] for j in range(i + 1, n))) / L_transpose[i][i]

    return x


def doolittle_LU_decomposition(A):
    n = len(A)
    L = [[0.0] * n for _ in range(n)]
    U = [[0.0] * n for _ in range(n)]

    for i in range(n):
        L[i][i] = 1.0

        for j in range(i, n):
            U[i][j] = A[i][j] - sum(L[i][k] * U[k][j] for k in range(i))

        for j in range(i + 1, n):
            L[j][i] = (A[j][i] - sum(L[j][k] * U[k][i] for k in range(i))) / U[i][i]

    return L, U


def forward_substitution_LU(L, b):
    n = len(L)
    y = [0.0] * n

    for i in range(n):
        y[i] = b[i] - sum(L[i][j] * y[j] for j in range(i))

    return y


def backward_substitution_LU(U, y):
    n = len(U)
    x = [0.0] * n

    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - sum(U[i][j] * x[j] for j in range(i + 1, n))) / U[i][i]

    return x


def is_symmetric_positive_definite(A):
    n = len(A)

    # Check for symmetry
    for i in range(n):
        for j in range(i + 1, n):
            if A[i][j] != A[j][i]:
                return False, "Matrix is not symmetric."

    # Check for positive definiteness
    try:
        cholesky_decomposition(A)
    except ValueError:
        return False, "Matrix is not positive definite."

    return True, "Matrix is symmetric and positive definite."


def solve_matrix_equation(A, b):
    symmetric_positive_definite, message = is_symmetric_positive_definite(A)
    if symmetric_positive_definite:
        print("Using Cholesky method.")
        L = cholesky_decomposition(A)
        y = forward_substitution(L, b)
        x = backward_substitution(L, y)
    else:
        print("Using Doolittle method.")
        L, U = doolittle_LU_decomposition(A)
        y = forward_substitution_LU(L, b)
        x = backward_substitution_LU(U, y)

    return x


# Problem 1
A1 = [[1, -1, 3, 2], [-1, 5, -5, -2], [3, -5, 19, 3], [2, -2, 3, 21]]
b1 = [15, -35, 94, 1]

print("Problem 1:")
solution1 = solve_matrix_equation(A1, b1)
print("Solution vector:", solution1)

# Problem 2
A2 = [[4, 2, 4, 0], [2, 2, 3, 2], [4, 3, 6, 3], [0, 2, 3, 9]]
b2 = [0, 2, 3, 9]

print("\nProblem 2:")
solution2 = solve_matrix_equation(A2, b2)
print("Solution vector:", solution2)
