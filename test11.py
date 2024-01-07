from sympy import Function, classify_ode, Eq, Derivative, exp
from sympy.abc import x

f = Function('f')
eq = classify_ode(Eq(f(x).diff(x), 0), f(x))

print(eq)
#eq = classify_ode(f(x).diff(x, 2) + 3*f(x).diff(x) + 2*f(x) - 4)
eq = classify_ode(Derivative(f,x) - f**2 * exp(-2*x) , 0)

print(eq)









