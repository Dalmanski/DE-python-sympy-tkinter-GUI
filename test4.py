from sympy import symbols, Function, dsolve, Eq, exp, sympify

# Define the variables
x, C = symbols('x C')
y = Function('y')(x)

# Input your differential equation
user_input = "dy/dx - y**2 * exp(-2*x) = 0"

# Replace 'dy/dx' with y'
equation_str = user_input.replace('dy/dx', "y'")

# Parse the user-provided equation
equation = Eq(sympify(equation_str), 0)

# Solve the differential equation
solution = dsolve(equation)

# Extract the solution expression
y_solution = solution.rhs

# Display the solution
print("\nOriginal Equation:")
print(equation)
print()

print("Solution:")
print(f"y = {y_solution}")

