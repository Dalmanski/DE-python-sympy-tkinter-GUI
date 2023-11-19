from sympy import symbols, Function, dsolve, Eq, exp, simplify, sympify, tan
import sympy as sp

# Define the variables
x, C = symbols('x C')
y = Function('y')(x)


problem = "y.diff(x) - y**2 * exp(-2*x) = 0"
#problem = "tan(x)*y.diff(x)-y=0"
#problem = "y.diff(x)-3*x**2*y=0"
#problem = "dy/dx-3*x**2*y=0"
#if "*dy/dx" in problem:
problem = problem.replace('dy/dx','')
problem = problem.replace('y.diff(x)','')
problem = problem.replace('y','y(x)')


# Initialize variables to store left-hand side and right-hand side of equations
eqlhs, eqrhs = problem.split('=')

# Parse the user-provided equations
eqlhs = sympify(eqlhs.strip())

# Define the right-hand side separately
eqrhs = sympify(eqrhs.strip())

#if "*dy/dx" in problem:
    #equation = Eq(eqlhs, eqrhs)
#else:
equation = Eq(y.diff(x) + eqlhs, eqrhs)

print("Original Equation:")
print(equation)
print(sp.pretty(equation))

# Solve the differential equation
solution = dsolve(equation)

# Extract the solution expression
y_solution = solution.rhs

# Display the solution
print("\nSolution:")
print(f"y = {y_solution}")
print(f"y = {sp.pretty(y_solution)}")






