import sympy as sp
from sympy import symbols, Function, Eq, dsolve, sympify, simplify_logic, exp, Derivative, diff, cos

# Define the variable
x = symbols('x')

# Define y as a function of x
y = Function('y')(x)

# Original expression
# dy/dx = 
# original_expression = Derivative(y, x) - y**2 * exp(-2*x)
original_expression = Eq(cos(x)**2 * diff(y, x) - x * y**2, 0)
#original_expression = "y.diff(x) - y**2 * exp(-2*x)"
#original_expression = original_expression.replace("y", "y(x)")
#original_expression = sympify(original_expression)
# Convert to the desired form
# original_expression = original_expression.subs(exp(2*x), exp(-2*x))

#original_expression = Eq(original_expression, 0)

seperableVariable = dsolve(original_expression)

# Display the result
print(original_expression)
print(sp.pretty(seperableVariable))





