import sympy as sp
import string
import re

def derivative():
    # Define the variable
    x = sp.symbols('x')

    # Ask the user to input their own function
    user_function = input("Enter your function in terms of 'x': ")
    #user_function = user_function.lower()

    # Check if the user's input is empty
    if not user_function:
        print("Error: Empty input. Please try again.")
        exit(1)

    # Replace and format the user's input
    replacements = {
        '^': '**',
        'c_1': 'C_1',
        'c_2': 'C_2',
        'C_1e': 'C_1*e',
        'C_2e': 'C_2*e',
        'C_1sin': 'C_1*sin',
        'C_1cos': 'C_1*cos',
        'C_1tan': 'C_1*tan',
        'C_2sin': 'C_2*sin',
        'C_2cos': 'C_2*cos',
        'C_2tan': 'C_2*tan',
        'sinx': 'sin(x)',
        'cosx': 'cos(x)',
        'tanx': 'tan(x)',
        '{': '(',
        '}': ')',
        'Asin': 'A*sin',
        'Acos': 'A*cos',
        'Atan': 'A*tan'
    }
    for find, replace in replacements.items():
        user_function = user_function.replace(find, replace)   

    # Read variables such as a-z and numbers -1000000 to 1000000
    #for i in string.ascii_lowercase:
    #   if i+'x' in user_function:
    #        user_function = user_function.replace(i+'x', i+'*x')
    for i in string.ascii_lowercase:
        if i+'x' in user_function:
            user_function = user_function.replace(i+'x', i+'*x')
        if i+'t' in user_function:
            user_function = user_function.replace(i+'t', '('+i+'*t)')
            x = sp.symbols('t')
    for i in range(-1000000, 1000001):
        if str(i)+'x' in user_function:
            user_function = user_function.replace(str(i)+'x', str(i)+'*x')

    # Check if 'y=' is in the input, and extract the expression after 'y='
    if 'y=' in user_function:
        y = sp.sympify(user_function.split('y=')[1])
    elif 'x=' in user_function:
        y = sp.sympify(user_function.split('x=')[1])
    else:
        y = sp.sympify(user_function)

    # Calculate the derivative
    dy_dx = sp.diff(y, x)

    # Convert the derivative to a plain string and remove the '*'
    plain_string_derivative = str(dy_dx)
    plain_string_derivative = plain_string_derivative.replace('log(e)', '')
    plain_string_derivative = plain_string_derivative.replace('**', '^')
    plain_string_derivative = plain_string_derivative.replace('*', '')
    plain_string_derivative = plain_string_derivative.replace('sinx', 'sin(x)')
    plain_string_derivative = plain_string_derivative.replace('cosx', 'cos(x)')
    plain_string_derivative = plain_string_derivative.replace('tanx', 'tan(x)')

    if 'y=' in user_function:
        plain_string_derivative = 'y\' = ' + plain_string_derivative
    if 'x=' in user_function:
        plain_string_derivative = 'x\' = ' + plain_string_derivative

    # Print the plain string derivative
    print(f"The derivative of {y} is {plain_string_derivative}")

derivative()