import sympy as sp
import string
import re
import tkinter as tk
from tkinter import ttk
import math
import tkinter.font as tkFont

howManyDerive = 0

def derive(problem):
    global dy_dx
    global eqlhs, eqrhs
    global equation1, equation2
    global howManyDerive
    # Define the variable
    x, y, k, t,  C_1, C_2, sin, cos, tan, e, y_prime, y_prime2 = sp.symbols("x y k t C_1 C_2 sin cos tan e y_prime y_prime2") 

    print(f'on {howManyDerive}, you have {problem}')

    # Check if the user's input is empty
    #if not problem:
    #    output_label.config(text="Error: Empty input. Please try again.")
    #    return

    # Define replacements
    if howManyDerive == 0:
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
                x = sp.symbols('t')
        for i in range(-1000000, 1000001):
            if str(i)+'x' in problem:
                problem = problem.replace(str(i)+'x', str(i)+'*x')
    

        # Initialize variables to store left-hand side and right-hand side of equations
        eqlhs, eqrhs = problem.split('=')

        # Parse the user-provided equations
        eqlhs = sp.sympify(eqlhs.strip())
        eqrhs = sp.sympify(eqrhs.strip())

    if howManyDerive == 0:
        # Add the two equations
        equation1 = sp.Eq(eqlhs, eqrhs)
        # Create the symbolic expression
        symbol_expr = sp.sympify(eqrhs)
    elif howManyDerive == 1:
        equation3 = sp.Eq(problem.lhs, problem.rhs)
        # Create the symbolic expression
        symbol_expr = sp.sympify(problem.rhs)  

    # Calculate the derivative
    dy_dx = sp.diff(symbol_expr, x)

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
    x, y, k, t, C_1, C_2, sin, cos, tan, e, y_prime, y_prime2 = sp.symbols("x y k t C_1 C_2 sin cos tan e y_prime y_prime2") 
    

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
            coefficient_C_1 = expression.coeff(C_2)

            # Convert the polynomial expression to a Sympy expression
            expr = sp.sympify(coefficient_C_1)

            # Find the coefficient of the leading term of the polynomial expression
            constant1 = expr.args[0]
            print(f'You use {constant1}')

            
            # Add the two equations
            equation = sp.Eq(equation2.lhs, equation2.rhs)

            # Create the symbolic expression
            symbol_expr = sp.sympify(equation2.rhs)
  
            # Extract the coefficient of C_1
            expression = equation.lhs - equation.rhs  # Move everything to the left-hand side
            coefficient_C_1 = expression.coeff(C_2)

            # Convert the polynomial expression to a Sympy expression
            expr = sp.sympify(coefficient_C_1)

            # Find the coefficient of the leading term of the polynomial expression
            constant2 = expr.args[0]
            print(f'You use {constant2}')

            constant = -1*(constant2 / constant1)
            print(f'So {constant}')


            equation4 = sp.Eq(constant*equation1.lhs + equation2.lhs, constant*equation1.rhs + equation2.rhs)
            print(equation1)
            print(equation2)
            print(equation4)
            

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
            coefficient_C_1 = expression.coeff(C_2)

            # Convert the polynomial expression to a Sympy expression
            expr = sp.sympify(coefficient_C_1)

            # Find the coefficient of the leading term of the polynomial expression
            constant2 = expr.args[0]
            print(f'You use {constant2}')

            
            # Add the two equations
            equation = sp.Eq(equation3.lhs, equation3.rhs)

            # Create the symbolic expression
            symbol_expr = sp.sympify(equation3.rhs)
  
            # Extract the coefficient of C_1
            expression = equation.lhs - equation.rhs  # Move everything to the left-hand side
            coefficient_C_1 = expression.coeff(C_2)

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
    coefficient_C_1 = expression.coeff(C_1)
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
    coefficient_C_1 = expression.coeff(C_1)

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


def main():
    root = tk.Tk()
    root.configure(bg='lightblue')
    root.title('Differential Equation')
    root.geometry("720x480")

    # Create a ComboBox
    style = ttk.Style()
    style.configure('TCombobox', relief='solid', highlightbackground='black', highlightcolor='black', highlightthickness=1, justify='center')
    combo = ttk.Combobox(root, values=('Derivative', 'Arbitrary Constant', 'Option 3', 'Option 4'), style='TCombobox')
     
    # Set items for the ComboBox
    #combo['values'] = ('Derivative', 'Abritary Constant', 'Option 3', 'Option 4')

    # Set the default value (optional)
    #combo.set('Derivative')
    combo.set('Arbitrary Constant')

    global dy_dx
    global eqlhs

    combo.pack(fill='x')
    combo.place(relx=0.5, rely=0.5, anchor='center', y=-140)  # Set the position to (10, 10)
    combo.configure(width=50)  # Set the width to 100

    # Function to handle the selection
    def on_select(event):
        global topic
        topic = combo.get()
        #print(f"Selected option: {topic}")

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

    # Create a button to trigger the action
    selected_item = combo.get()
    if selected_item == "Derivative":
        button = tk.Button(root, text="Compute", command=Derivative, relief='raised', highlightthickness=1)
    elif selected_item == "Arbitrary Constant":
        button = tk.Button(root, text="Compute", command=Arbitrary_Constant, relief='raised', highlightthickness=1)
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

