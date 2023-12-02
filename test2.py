import sympy as sp

# Define the symbolic variables
x, C1, C2 = sp.symbols('x C1 C2')
y = sp.Function('y')(x)

# Define the equations
eq1 = sp.Eq(y, C1 * sp.exp(2 * x) + C2 * sp.exp(-3 * x))
eq2 = sp.Eq(sp.diff(y, x), 2 * C1 * sp.exp(2 * x) - 3 * C2 * sp.exp(-3 * x))

# Combine like terms
combined_eq = eq1 + eq2

# Solve for C2
C2 = sp.solve(combined_eq, C2)[0]

# Substitute C2 back into the first equation
eq1_simplified = eq1.subs(C2, C2)

# Print the simplified equation
print(sp.pretty(eq1_simplified))












