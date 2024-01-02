import sympy as sp
from sympy import symbols, Function, Eq, dsolve, sympify, simplify_logic, exp, Derivative

# Define the variable
x = symbols('x')

# Define y as a function of x
y = Function('y')(x)

# Original expression
# dy/dx = 
# original_expression = Derivative(y, x) - y**2 * exp(-2*x)
original_expression = Eq(Derivative(y,x), x**2/y**2)
#original_expression = "y.diff(x) - y**2 * exp(-2*x)"
#original_expression = original_expression.replace("y", "y(x)")
#original_expression = sympify(original_expression)
# Convert to the desired form
# original_expression = original_expression.subs(exp(2*x), exp(-2*x))

#original_expression = Eq(original_expression, 0)

seperableVariable = dsolve(original_expression)

# Display the result
print(sp.pretty(seperableVariable))





