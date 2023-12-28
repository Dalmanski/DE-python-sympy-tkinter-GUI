import sympy as sp
from sympy import sympify, Eq, Derivative

def find_iv_dv(equation_str):
    try:
        eqlhs, eqrhs = equation_str.split('=')
        # Parse the differential equation
        eqlhs = sympify(eqlhs.strip())
        eqrhs = sympify(eqrhs.strip())
        equation = Eq(eqlhs, eqrhs)

        # Extract variables from derivatives on the left-hand side
        lhs_derivative_vars = [var for var in eqlhs.free_symbols if isinstance(Derivative(eqlhs, var), Derivative)]

        # Extract variables from the right-hand side
        rhs_vars = [var for var in equation.rhs.free_symbols]

        # Combine all variables, prioritizing those from derivatives
        all_vars = sorted(set(lhs_derivative_vars + rhs_vars), key=lambda sym: sym.name)

        if len(all_vars) == 2:
            independent_variable, dependent_variable = all_vars
            return independent_variable, dependent_variable
    except sp.SympifyError:
        print("Error parsing the equation.")

    return None

# Example usage
equation = input("Enter your differential equation: ")
result = find_iv_dv(equation)

if result:
    iv, dv = result
    print(f"Independent variable (IV): {iv}")
    print(f"Dependent variable (DV): {dv}")
else:
    print("Unable to identify independent and dependent variables. Please check your equation.")

