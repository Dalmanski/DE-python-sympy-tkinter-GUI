import sympy as sp

# Define the variables
x, C = sp.symbols('x C')

# Define the differential equation
equation = sp.Eq(sp.diff(sp.Symbol('y')(x), x) - sp.Symbol('y')(x)**2 * sp.exp(-2*x), 0)

# Solve the differential equation
solution = sp.dsolve(equation)

# Extract the solution expression
y_solution = solution.rhs

# Display the solution
print("Solution:")
print(f"y = {y_solution}")
print(f"y = {sp.pretty(y_solution)}")




