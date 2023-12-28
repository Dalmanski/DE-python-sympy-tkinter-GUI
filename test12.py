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
        lhs_derivative_vars = []
        for term in eqlhs.args:
            if isinstance(term, Derivative):
                lhs_derivative_vars.extend(term.free_symbols)

        # Extract variables from the right-hand side
        rhs_vars = [var for var in equation.rhs.free_symbols]

        # Combine all variables, prioritizing those from derivatives
        all_vars_set = set(lhs_derivative_vars + rhs_vars)

        # Maintain the order of variables based on their first appearance
        encountered_vars = []
        unique_vars = []
        for var in all_vars_set:
            if var not in encountered_vars:
                encountered_vars.append(var)
                unique_vars.append(var)

        if len(unique_vars) == 2:
            independent_variable = unique_vars[0]
            dependent_variable = unique_vars[1]
            print(f"IV = {independent_variable}, DV = {dependent_variable}")
            return independent_variable, dependent_variable
    except sp.SympifyError:
        print("Error parsing the equation.")

    return None, None  # Ensure a tuple is returned in case of an error

# Example usage
equation = input("Enter your differential equation: ")
iv, dv = find_iv_dv(equation)

print(f"IV = {iv}, DV = {dv}")

    #print("Unable to identify independent and dependent variables. Please check your equation.")







