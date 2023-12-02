

import sympy as sp
from sympy import symbols, Function, Eq, dsolve, sympify, simplify_logic, exp, Derivative

# Define the variables
x, y, k, t, a, c, C_1, C_2, sin, cos, tan, e, y_prime, y_prime2 = sp.symbols("x y k t a c C_1 C_2 sin cos tan e y_prime y_prime2") 

y = Function('y')(x)

mode = 'add'
#problem = "y.diff(x) - y**2 * exp(-2*x) = 0"
problem = "tan(x)*y.diff(x)-y=0"
#problem = "y.diff(x)-3*x**2*y=0"
#problem = "dy/dx-3*x**2*y=0"
#if "*dy/dx" in problem:
problem = problem.replace(' ','')
problem = problem.replace('dy/dx','y.diff(x)')
if 'y.diff(x)*' in problem or '*y.diff(x)' in problem:
    mode = 'multiply'
print(f'{mode} the dx/dy')
problem = problem.replace('dy/dx','')
problem = problem.replace('y.diff(x)','')
problem = problem.replace('y','y(x)')


# Initialize variables to store left-hand side and right-hand side of equations
eqlhs, eqrhs = problem.split('=')
print(f"before eqlhs = {eqlhs}, eqrhs = {eqrhs}")

# Parse the user-provided equations
eqlhs = sympify(eqlhs.strip())

    
# Define the right-hand side separately
eqrhs = sympify(eqrhs.strip())

#if "*dy/dx" in problem:
    #equation = Eq(eqlhs, eqrhs)
#else:
if mode == 'multiply':
    equation = Eq(Derivative(y, x) * eqlhs, eqrhs)
    str_eqlhs = str(equation.lhs)
    if 'y(x)*' in str_eqlhs:
        str_eqlhs = str_eqlhs.replace('y(x)*','y(x)+')
        eqlhs = sympify(str_eqlhs)
        equation = Eq(eqlhs, eqrhs)
else:
    equation = Eq(y.diff(x) + eqlhs, eqrhs)


print(f"eqlhs = {eqlhs}, eqrhs = {eqrhs}")
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






