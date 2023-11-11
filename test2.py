import sympy as sp

x, y, z, c_1, c_2, e, y_prime = sp.symbols("x, y, z, c_1, c_2, e, y_prime")

# Define the equations correctly
eq1 = "3 * y = 3 * c_1 * e**(2*x) + 3 * c_2 * e**(-3*x)"
eq2 = "y_prime = 2 * c_1 * e**(2*x) - 3 * c_2 * e**(-3*x)"

# Create SymPy equations
equation1 = sp.Eq(eq1, 0)
equation2 = sp.Eq(eq2, 0)

# Add the two equations
combined_equation = sp.Eq(equation1.lhs + equation2.lhs, equation1.rhs + equation2.rhs)

# Print the combined equation
sp.pretty(combined_equation)
print(combined_equation)



