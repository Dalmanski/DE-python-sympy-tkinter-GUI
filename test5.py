from sympy import symbols, Function, dsolve, Eq, exp, tan, sympify, parsing
import sympy as sp



# Define the variables
x, C = symbols('x C')
y = Function('y')(x)

#def parse_with_diff(expr_str):
#    """Custom parsing function that handles diff(x) symbol."""
#    expr_parts = expr_str.split()
#    for i, part in enumerate(expr_parts):
#        if part == "diff(x)":
#            expr_parts[i] = sp.Derivative(sp.Symbol("y"), x)
#    return sp.sympify(" ".join(expr_parts))

#equationStr = "y.diff(x) - y**2 * exp(-2*x) = 0"
#eqlhs, eqrhs = equationStr.split('=')

#eqlhs = parse_with_diff(eqlhs)
#eqrhs = sp.sympify(eqrhs)


#equation = Eq(eqlhs, eqrhs)

#equation = Eq(y.diff(x) - y**2 * exp(-2*x), 0)
equation = Eq(x*(x+1)*y.diff(x)-y-1,0)
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
print(sp.pretty(solution))









