from sympy import symbols, Function, dsolve, Eq, exp, tan, sympify, parsing
import sympy as sp
from sympy.abc import x



# Define the variables
x, C, sin = symbols('x C sin')
#y = Function('y')(x)
y = sp.Function('y')
#def parse_with_diff(expr_str):
#    """Custom parsing function that handles diff(x) symbol."""
#    expr_parts = expr_str.split()
#    for i, part in enumerate(expr_parts):
#        if part == "diff(x)":
#            expr_parts[i] = sp.Derivative(sp.Symbol("y"), x)
#    return sp.sympify(" ".join(expr_parts))

#equationStr = "dy_dx - y**2 * exp(-2*x) = 0"
#eqlhs, eqrhs = equationStr.split('=')

#eqlhs = sp.sympify(eqlhs)
#eqrhs = sp.sympify(eqrhs)


#equation = sp.Eq(eqlhs, eqrhs)
#dy_dx = sp.simplify(sp.diff(equation, x))
equation = Eq(y.diff(x), sin(5)*x)
#equation = Eq((1*x)*y.diff(x)-x*y,x+x**2)
#equation = Eq(tan(x) * y.diff(x) - y, 0)
#equation = Eq(y.diff(x) - 3*x**2*y, 0)
print(equation)
print(sp.pretty(equation))

# Solve the differential equation
solution = sp.dsolve(equation)

# Extract the solution expression
y_solution = solution.rhs

# Display the solution
print("Solution:")
print(f"y = {y_solution}")
print(sp.pretty(solution))









