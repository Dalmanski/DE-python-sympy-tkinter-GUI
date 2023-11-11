import sympy as sp

x, y, C_1, C_2, e = sp.symbols("x y C_1 C_2 e")

# Define the equation
equation = input("Enter your problem:")
#equation = sp.Eq(y,10*C_1*e**(-2*x)+20*C_2*e**(3*x))

 # Initialize variables to store left-hand side and right-hand side of equations
eqlhs, eqrhs = equation.split('=')

# Parse the user-provided equations
eqlhs = sp.sympify(eqlhs.strip())
eqrhs = sp.sympify(eqrhs.strip())

equation = sp.Eq( eqlhs, eqrhs)

print(f"The equation is: {equation}")


