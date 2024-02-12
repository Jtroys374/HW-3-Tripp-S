import math


def t_distribution_probability(df, z):
    if df <= 0:
        raise ValueError("Degrees of freedom must be greater than 0.")

    if z == 0:
        return 0.5

    def gamma_function(x):
        return math.gamma(x)

    def integral_function(t, df):
        return ((1 + t ** 2 / df) ** (- (df + 1) / 2)) / (math.sqrt(df) * gamma_function(df / 2) * 2 ** (df / 2 - 1))

    n = 1000  # Number of subdivisions for Simpson's rule
    a = 0
    b = z
    h = (b - a) / n

    integral_sum = (integral_function(a, df) + integral_function(b, df))

    for i in range(1, n, 2):
        integral_sum += 4 * integral_function(a + i * h, df)

    for i in range(2, n - 1, 2):
        integral_sum += 2 * integral_function(a + i * h, df)

    return integral_sum * h / 3


def main():
    df = int(input("Enter degrees of freedom: "))
    z_values = [float(input(f"Enter z value {i + 1}: ")) for i in range(3)]

    for z in z_values:
        probability = t_distribution_probability(df, z)
        print(f"For z = {z} and {df} degrees of freedom, the probability is: {probability:.6f}")


if __name__ == "__main__":
    main()
