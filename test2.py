from sympy import sympify, SympifyError

def is_valid_sympy_expression(expression_str):
    try:
        sympify(expression_str)
        return True
    except SympifyError:
        return False

# Examples
expression1 = "C_1*exp(2*x) + C_2*exp(-3*x)"
expression2 = "dx"

print(f"Expression 1 is valid: {is_valid_sympy_expression(expression1)}")
print(f"Expression 2 is valid: {is_valid_sympy_expression(expression2)}")













