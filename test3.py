from sympy import symbols, Function, exp, Eq, simplify

# Define the variable
x = symbols('x')

# Define y as a function of x
y = Function('y')(x)

# Original expression
original_expression = -y**2 / exp(2*x)

# Convert to the desired form
original_expression = original_expression.subs(exp(2*x), exp(-2*x))

# Display the result
print(original_expression)





