"""Small exact polynomials and sufficient Bernstein positivity certificates."""

import math
from fractions import Fraction
from itertools import permutations


def polynomial(coefficients):
    result = tuple(Fraction(x) for x in coefficients)
    while len(result) > 1 and result[-1] == 0:
        result = result[:-1]
    return result or (Fraction(0),)


def add(a, b):
    return polynomial(
        (a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
        for i in range(max(len(a), len(b)))
    )


def scale(a, value):
    return polynomial(x * value for x in a)


def subtract(a, b):
    return add(a, scale(b, -1))


def multiply(a, b):
    result = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i + j] += x * y
    return polynomial(result)


def evaluate(a, x):
    value = Fraction(0)
    for coefficient in reversed(a):
        value = value * x + coefficient
    return value


def determinant(matrix):
    n = len(matrix)
    total = polynomial([0])
    for order in permutations(range(n)):
        inversions = sum(order[i] > order[j] for i in range(n) for j in range(i + 1, n))
        term = polynomial([(-1) ** inversions])
        for i, j in enumerate(order):
            term = multiply(term, matrix[i][j])
        total = add(total, term)
    return total


def adjugate(matrix):
    n = len(matrix)
    return [
        [
            scale(
                determinant(
                    [
                        [matrix[row][col] for col in range(n) if col != i]
                        for row in range(n)
                        if row != j
                    ]
                ),
                (-1) ** (i + j),
            )
            for j in range(n)
        ]
        for i in range(n)
    ]


def remove_stationary_endpoint_factors(a):
    """Divide by 1-x and 1+x while exact; both are positive on (-1,1)."""
    a = polynomial(a)
    for endpoint in (Fraction(1), Fraction(-1)):
        while len(a) > 1 and evaluate(a, endpoint) == 0:
            quotient = [Fraction(0)] * (len(a) - 1)
            quotient[-1] = a[-1]
            for i in range(len(quotient) - 2, -1, -1):
                quotient[i] = a[i + 1] + endpoint * quotient[i + 1]
            # Synthetic division gives x-endpoint; flip x-1 to 1-x.
            a = scale(quotient, -1 if endpoint == 1 else 1)
    return a


def bernstein(a, left, right):
    """Coefficients on [left,right], computed with exact rational arithmetic."""
    left, right = Fraction(left), Fraction(right)
    if left >= right:
        raise ValueError('Require left < right')
    degree = len(a) - 1
    transformed = [
        sum(a[j] * math.comb(j, i) * left ** (j - i) for j in range(i, degree + 1))
        * (right - left) ** i
        for i in range(degree + 1)
    ]
    return tuple(
        sum(transformed[i] * Fraction(math.comb(k, i), math.comb(degree, i)) for i in range(k + 1))
        for k in range(degree + 1)
    )


def positive_on(a, left, right):
    """Certify strict positivity, excluding only stationary endpoints -1,+1."""
    coefficients = bernstein(a, left, right)
    return (
        all(x >= 0 for x in coefficients)
        and any(x > 0 for x in coefficients)
        and (left == -1 or coefficients[0] > 0)
        and (right == 1 or coefficients[-1] > 0)
    )
