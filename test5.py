from sympy import symbols, Function, dsolve, Eq, exp, tan, sympify
import sympy as sp


# Define the variables
x, C = symbols('x C')
y = Function('y')(x)

# Define the differential equation

equation = Eq(y.diff(x) - y**2 * exp(-2*x), 0)
#equation = Eq(tan(x) * y.diff(x) - y, 0)
#equation = Eq(y.diff(x) - 3*x**2*y, 0)
print(equation)
print(sp.pretty(equation))

# Solve the differential equation
solution = dsolve(equation)

# Extract the solution expression
y_solution = solution.rhs

# Display the solution
print("Solution:")
print(f"y = {y_solution}")
print(f"y = {sp.pretty(y_solution)}")









