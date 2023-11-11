import sympy as sp

x, y, C_1, C_2, e = sp.symbols("x y C_1 C_2 e")

# Define the equation
equation = sp.Eq(y,10*C_1*e**(-2*x)+20*C_2*e**(3*x))

# Extract the coefficient of C_1
expression = equation.lhs - equation.rhs  # Move everything to the left-hand side
coefficient_C_1 = expression.coeff(C_2)

print(f"The constant factor of C_1 is: {coefficient_C_1}")

# Convert the polynomial expression to a Sympy expression
expr = sp.sympify(coefficient_C_1)

# Find the coefficient of the leading term of the polynomial expression
constant = expr.args[0]

# Print the constant
print(constant)