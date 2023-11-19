import sympy as sp
import string
import re
import tkinter as tk
from tkinter import ttk
import math
import tkinter.font as tkFont
import numpy as np
from sympy import symbols, Function, dsolve, Eq, exp, simplify, sympify, tan

howManyDerive = 0
x, y, k, t, a, c, C_1, C_2, sin, cos, tan, e, y_prime, y_prime2 = sp.symbols("x y k t a c C_1 C_2 sin cos tan e y_prime y_prime2") 

def find_variable(expr):
    for term in expr.args:
        if isinstance(term, sp.Mul):
            variables = [symbol for symbol in term.free_symbols if symbol != C_1 and symbol != C_2 and symbol != k and symbol != e]
            if variables:
                return variables[0]
            
def replacement(problem):
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
        problem = problem.replace(find, replace)

    # Read variables such as a-z and numbers -1000000 to 1000000
    for i in string.ascii_lowercase:
        if i+'x' in problem:
            problem = problem.replace(i+'x', i+'*x')
        if i+'t' in problem:
            problem = problem.replace(i+'t', '('+i+'*t)')
            #x = sp.symbols('t')
    for i in range(-1000000, 1000001):
        if str(i)+'x' in problem:
            problem = problem.replace(str(i)+'x', str(i)+'*x')
        if str(i)+'t' in problem:
            problem = problem.replace(str(i)+'t', str(i)+'*t')
    return problem

def derive(problem):
    global dy_dx
    global eqlhs, eqrhs
    global equation1, equation2
    global howManyDerive
    # Define the variable
    global x, y, k, t, a, c, C_1, C_2, sin, cos, tan, e, y_prime, y_prime2 

    # Check if the user's input is empty
    #if not problem:
    #    output_label.config(text="Error: Empty input. Please try again.")
    #    return

    # Define replacements
    if howManyDerive == 0:
        
        problem = replacement(problem)
    
        print(f'on {howManyDerive}, you have {problem}')

        # Initialize variables to store left-hand side and right-hand side of equations
        eqlhs, eqrhs = problem.split('=')

        # Parse the user-provided equations
        eqlhs = sp.sympify(eqlhs.strip())
        eqrhs = sp.sympify(eqrhs.strip())

        variable = find_variable(eqrhs)
        print(f'variable = {variable}')
        x = sp.symbols(str(variable))

    if howManyDerive == 0:
        # Add the two equations
        equation1 = sp.Eq(eqlhs, eqrhs)
        # Create the symbolic expression
        symbol_expr = sp.sympify(eqrhs)
        print(f'left side = {eqlhs}')
        print(f'right side = {eqrhs}')
    elif howManyDerive == 1:
        # Create the symbolic expression
        symbol_expr = sp.sympify(problem.rhs)  
        print(f'right side = {symbol_expr}')

    # Calculate the derivative
    print(f'variable = {x}')
    dy_dx = sp.diff(symbol_expr, x)
    print(f'dy_dx = {dy_dx}')


def Derivative():
    # Retrieve user input from the entry widget
    user_function = entry.get()
    derive(user_function)
    output()


def Arbitrary_Constant():

    global eqlhs
    global equation1
    global howManyDerive
    howManyDerive = 0
    # Define the variable
    global x, y, k, t, a, c, C_1, C_2, sin, cos, tan, e, y_prime, y_prime2 
    
    elimination = entryEliminate.get()

    # Split by spaces
    split_elimination = elimination.split()

    # List to store student names
    eliminations = []
    # Add the split to the list
    for name in split_elimination:
        eliminations.append(name)
    
    for howManyDerive in range(2):
        if howManyDerive == 0:
            # Retrieve user input from the entry widget
            user_function = entry.get()
            derive(user_function)
            equation2 = sp.Eq(y_prime, dy_dx)
            equation2 = equation2.subs(sp.log(e), 1)
        
            # Add the two equations
            equation = sp.Eq(equation1.lhs, equation1.rhs)
            # Create the symbolic expression
            symbol_expr = sp.sympify(equation1.rhs)
            
            # Extract the coefficient of C_1
            expression = equation.lhs - equation.rhs  # Move everything to the left-hand side
            print(expression)
            if len(eliminations) == 2:
                coefficient_C_1 = expression.coeff(eliminations[1])
                print(f'DETECTED 2 elimination so use on {eliminations[1]}')
            else:
                coefficient_C_1 = expression.coeff(eliminations[0])
                print(f'DETECTED 1 elimination {eliminations[0]}')
            print(f'coefficient = {coefficient_C_1}')
            #coefficient_C_1 = expression.coeff(C_2)

            #                                       AYOHA NI NAUNSA MNI UY DUKA KAAYO KO
            if 'C_1' or 'C_2' in coefficient_C_1:
                # Convert the polynomial expression to a Sympy expression
                expr = sp.sympify(coefficient_C_1)
                print(f'expr = {expr}')

                # Find the coefficient of the leading term of the polynomial expression
                constant1 = expr.args[0]
            else:
                constant1 = coefficient_C_1
                print('You convert')
            
            print(f'You use {constant1}')
            
            # Add the two equations
            equation = sp.Eq(equation2.lhs, equation2.rhs)

            # Create the symbolic expression
            symbol_expr = sp.sympify(equation2.rhs)
  
            # Extract the coefficient of C_1
            expression = equation.lhs - equation.rhs  # Move everything to the left-hand side

            if len(eliminations) == 2:
                coefficient_C_1 = expression.coeff(eliminations[1])
                print(f'DETECTED 2 elimination so use {eliminations[1]}')
            else:
                coefficient_C_1 = expression.coeff(eliminations[0])
                print(f'DETECTED 1 elimination {eliminations[0]}')

            # Convert the polynomial expression to a Sympy expression
            expr = sp.sympify(coefficient_C_1)
            # Find the coefficient of the leading term of the polynomial expression
            constant2 = expr.args[0]
            print(f'You use {constant2}')

            constant = -1*(constant2 / constant1)
            print(f'So {constant}')


            equation4 = sp.Eq(constant*equation1.lhs + equation2.lhs, constant*equation1.rhs + equation2.rhs)
            print(f'equation 1: {equation1}')
            print(f'equation 2: {equation2}')
            print(f'equation 4: {equation4}')
            

        elif howManyDerive == 1:
            derive(equation2)
            equation3 = sp.Eq(y_prime2, dy_dx)
            equation3 = equation3.subs(sp.log(e), 1)
            
            # Add the two equations
            equation = sp.Eq(equation2.lhs, equation2.rhs)

            # Create the symbolic expression
            symbol_expr = sp.sympify(equation2.rhs)
  
            # Extract the coefficient of C_1
            expression = equation.lhs - equation.rhs  # Move everything to the left-hand side

            if len(eliminations) == 2:
                coefficient_C_1 = expression.coeff(eliminations[1])
                print(f'DETECTED 2 elimination so use {eliminations[1]}')
            else:
                coefficient_C_1 = expression.coeff(eliminations[0])
                print(f'DETECTED 1 elimination {eliminations[0]}')

            # Convert the polynomial expression to a Sympy expression
            expr = sp.sympify(coefficient_C_1)
            print(expr)
            # Find the coefficient of the leading term of the polynomial expression
            constant2 = expr.args[0]
            print(f'You use {constant2}')

            
            # Add the two equations
            equation = sp.Eq(equation3.lhs, equation3.rhs)

            # Create the symbolic expression
            symbol_expr = sp.sympify(equation3.rhs)
  
            # Extract the coefficient of C_1
            expression = equation.lhs - equation.rhs  # Move everything to the left-hand side

            if len(eliminations) == 2:
                coefficient_C_1 = expression.coeff(eliminations[1])
                print(f'DETECTED 2 elimination so use {eliminations[1]}')
            else:
                coefficient_C_1 = expression.coeff(eliminations[0])
                print(f'DETECTED 1 elimination {eliminations[0]}')

            # Convert the polynomial expression to a Sympy expression
            expr = sp.sympify(coefficient_C_1)

            # Find the coefficient of the leading term of the polynomial expression
            constant3 = expr.args[0]
            print(f'You use {constant3}')

            constant = -1*(constant3 / constant2)
            print(f'So {constant}')

            equation5 = sp.Eq(constant*equation2.lhs + equation3.lhs, constant*equation2.rhs + equation3.rhs)
            print(f'Equation 5: {equation5}')

    # Add the two equations
    equation = sp.Eq(equation5.lhs, equation5.rhs)
    print(equation)

    # Create the symbolic expression
    symbol_expr = sp.sympify(equation5.rhs)
    print(symbol_expr)

    # Extract the coefficient of C_1
    expression = equation.lhs - equation.rhs  # Move everything to the left-hand side

    coefficient_C_1 = expression.coeff(eliminations[0])

    #coefficient_C_1 = expression.coeff(C_1)

    print(expression)
    print(coefficient_C_1)

    # Convert the polynomial expression to a Sympy expression
    expr = sp.sympify(coefficient_C_1)
    print(expr)

    # Find the coefficient of the leading term of the polynomial expression
    constant5 = expr.args[0]
    print(f'You use {constant5}')

            
    # Add the two equations
    equation = sp.Eq(equation4.lhs, equation4.rhs)

    # Create the symbolic expression
    symbol_expr = sp.sympify(equation4.rhs)
  
    # Extract the coefficient of C_1
    expression = equation.lhs - equation.rhs  # Move everything to the left-hand side

    coefficient_C_1 = expression.coeff(eliminations[0])

    # Convert the polynomial expression to a Sympy expression
    expr = sp.sympify(coefficient_C_1)

    # Find the coefficient of the leading term of the polynomial expression
    constant4 = expr.args[0]
    print(f'You use {constant4}')

    constant = -1*(constant5 / constant4)
    print(f'So {constant}')

    Answer = sp.Eq( equation5.lhs + constant*equation4.lhs, equation5.rhs + constant*equation4.rhs)
    print(Answer)
    result_text = f"Solution:"
    result_text += f"\nEquation 1 :\n{sp.pretty(equation1)}"
    result_text += f"\nEquation 2 :\n{sp.pretty(equation2)}"
    result_text += f"\nEquation 3 :\n{sp.pretty(equation3)}"
    result_text += f"\nEquation 4 :\n{sp.pretty(equation4)}"
    result_text += f"\nEquation 5 :\n{sp.pretty(equation5)}"
    result_text += f"\nThe arbitrary Constant is: \n {sp.pretty(Answer)}"
    output_label.config(text=result_text)
    #output()

def Separable_Variables():

    # Define the variables
    x, C = symbols('x C')
    y = Function('y')(x)
    global k, t, a, c, C_1, C_2, sin, cos, tan, e, y_prime, y_prime2 

    problem = entry.get()
    problem = problem.replace('dy/dx','')
    problem = problem.replace('y.diff(x)','')
    problem = replacement(problem)
    problem = problem.replace('y','y(x)')

    #problem = "y.diff(x) - y**2 * exp(-2*x) = 0"
    #problem = "tan(x)*y.diff(x)-y=0"
    #problem = "y.diff(x)-3*x**2*y=0"
    #problem = "dy/dx-3*x**2*y=0"
    #if "*dy/dx" in problem:

    # Initialize variables to store left-hand side and right-hand side of equations
    eqlhs, eqrhs = problem.split('=')

    # Parse the user-provided equations
    eqlhs = sympify(eqlhs.strip())

    # Define the right-hand side separately
    eqrhs = sympify(eqrhs.strip())

    #if "*dy/dx" in problem:
        #equation = Eq(eqlhs, eqrhs)
    #else:
    equation = Eq(y.diff(x) + eqlhs, eqrhs)
    result_text = f"Solution:"
    result_text += f"\n\nEquation :\n{sp.pretty(equation)}"
    print("Original Equation:")
    print(equation)
    print(sp.pretty(equation))

    # Solve the differential equation
    solution = dsolve(equation)
    result_text += f"\n\nFinal answer:\n{sp.pretty(solution)}"
    # Extract the solution expression
    y_solution = solution.rhs

    # Display the solution
    print("\nSolution:")
    print(f"y = {y_solution}")
    print(f"y = {sp.pretty(y_solution)}")
    output_label.config(text=result_text)

def output():
     # Convert the derivative to a plain string and format it
    plain_string_derivative = str(dy_dx)
    plain_string_derivative = plain_string_derivative.replace('log(e)', '')
    plain_string_derivative = plain_string_derivative.replace('**', '^')
    plain_string_derivative = plain_string_derivative.replace('*', '')
    plain_string_derivative = plain_string_derivative.replace('sinx', 'sin(x)')
    plain_string_derivative = plain_string_derivative.replace('cosx', 'cos(x)')
    plain_string_derivative = plain_string_derivative.replace('tanx', 'tan(x')

    # Display the result
    #result_text = f"The derivative is: \n{variable_name}={plain_string_derivative}"
    result_text = f"The derivative is: \n {eqlhs} = {sp.pretty(plain_string_derivative)}"
    output_label.config(text=result_text)



#FRONT END, diri mo mag design2
def main():
    # diri mo himo na mura FRONT END or design
    root = tk.Tk()
    root.configure(bg='lightblue')
    root.title('Differential Equation')
    root.geometry("720x480")

    # Create a ComboBox
    style = ttk.Style()
    style.configure('TCombobox', relief='solid', highlightbackground='black', highlightcolor='black', highlightthickness=1, justify='center')
    combo = ttk.Combobox(root, values=('Derivative', 'Arbitrary Constant', 'Separable Variables', 'Option 4'), style='TCombobox')
     
    # Set items for the ComboBox
    #combo['values'] = ('Derivative', 'Abritary Constant', 'Option 3', 'Option 4')

    # Set the default value (optional)
    #combo.set('Derivative')
    combo.set('Derivative')
    selected_item = combo.get()

    global dy_dx, eqlhs

    combo.pack(fill='x')
    combo.place(relx=0.5, rely=0.5, anchor='center', y=-140)  # Set the position to (10, 10)
    combo.configure(width=50)  # Set the width to 100

    # Function to handle the selection

    def on_select(event):
        global topic
        topic = combo.get()
        if topic == "Arbitrary Constant":
            entryEliminate.grid(row=1, column=0, padx=300, pady=150)
        else:
            entryEliminate.grid_forget()
        #print(f"Selected option: {topic}")

    def button_clicked():
        selected_item = combo.get()
        if selected_item == "Derivative":
            Derivative()
        elif selected_item == "Arbitrary Constant":
            Arbitrary_Constant()
        elif selected_item == "Separable Variables":
            Separable_Variables()
        # Add more conditions for other options if needed

    # Bind the function to the <<ComboboxSelected>> event
    combo.bind("<<ComboboxSelected>>", on_select)

    # Create a custom font
    custom_font = ('Arial', 32, 'bold italic')
    custom_font2 = tkFont.Font(family="Computer Modern", size=10)  # Set your desired font family and size

    # Create a label with the specified font and text
    labelTitle = tk.Label(root, text="DIFFERENTIAL EQUATION", font=custom_font)
    labelTitle.pack(fill='x')
    labelTitle.place(relx=0.5, rely=0.5, anchor='center', y=-200)
    labelTitle.configure(width=50)

    # Create an entry widget for user input
    global entry
    entry = tk.Entry(root, relief='solid', highlightthickness=1, justify='center')
    entry.pack(fill='x')
    entry.place(relx=0.5, rely=0.5, anchor='center', y=-100)  # Set the position to (10, 10)
    entry.configure(width=100)  # Set the width to 100

     # Create an entry widget for user input arbitary constant eliminate
    global entryEliminate
    entryEliminate = tk.Entry(root, relief='solid', highlightthickness=1, justify='center')
    entryEliminate.pack(fill='x')
    entryEliminate.place(relx=0.5, rely=0.5, anchor='center', y=-80)  # Set the position to (10, 10)
    entryEliminate.grid(row=1, column=0, padx=300, pady=150)
    entryEliminate.grid_forget()

    # Create a button to trigger the action
    button = tk.Button(root, text="Compute", command=button_clicked, relief='raised', highlightthickness=1)
    button.pack(fill='x')
    button.place(relx=0.5, rely=0.5, anchor='center', y=-60)
    button.configure(width=50)

    # Create a label to display the output
    global output_label
    output_label = tk.Label(root, text="", relief='flat', highlightthickness=1, font=custom_font2)
    output_label.pack(side='top', fill='x')
    output_label.pack(fill='x')
    output_label.place(relx=0.5, rely=0.5, anchor='center', y=90)
    output_label.configure(width=50)

    root.mainloop()

if __name__ == '__main__':
    main()

