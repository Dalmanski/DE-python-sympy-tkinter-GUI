import sympy as sp
import string
import re
import tkinter as tk
import math
import tkinter.font as tkFont
import numpy as np
from tkinter import ttk
from sympy import symbols, Function, dsolve, Eq, exp, simplify, sympify, tan, sin, cos, Derivative, SympifyError
from PIL import Image, ImageTk

howManyDerive = 0
x, y, k, t, a, c, C_1, C_2, sin, cos, tan, e, y_prime, y_prime2 = sp.symbols("x y k t a c C_1 C_2 sin cos tan e y\' y\'\'") 


def find_variable(expr):
    for term in expr.args:
        if isinstance(term, sp.Mul):
            variables = [symbol for symbol in term.free_symbols if symbol != C_1 and symbol != C_2 and symbol != k and symbol != e]
            if variables:
                return variables[0]


def replacement(problem):
    replacements = {
        ' ': '',
        '^': '**',
        'c_1': 'C_1',
        'c_2': 'C_2',
        'e^': 'exp',
        '{': '(',
        '}': ')'
    }
        
    for find, replace in replacements.items():
        problem = problem.replace(find, replace)

    # Read variables such as a-z and numbers -1000000 to 1000000
    #for i in string.ascii_lowercase:
    #    if i+'x' in problem:
    #        problem = problem.replace(i+'x', i+'*x')
    #    if i+'t' in problem:
    #        problem = problem.replace(i+'t', '('+i+'*t)')
    #        #x = sp.symbols('t')
    #for i in range(-1000000, 1000001):
    #    if str(i)+'x' in problem:
    #        problem = problem.replace(str(i)+'x', str(i)+'*x')
    #    if str(i)+'t' in problem:
    #        problem = problem.replace(str(i)+'t', str(i)+'*t')
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


def Derivatives():
    # Retrieve user input from the entry widget
    user_function = entry.get()
    derive(user_function)
    plain_string_derivative = str(dy_dx)
    # Display the result
    result_text = f"The derivative is: \n {sp.pretty(eqlhs)} = {sp.pretty(plain_string_derivative)}"
    output_label.config(text=result_text)
    final_label.config(text=result_text)


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

    # List to store
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
    x, y, k, t, a, c, C_1, C_2, sin, cos, tan, e, y_prime, y_prime2 = sp.symbols("x y k t a c C_1 C_2 sin cos tan e y_prime y_prime2") 
    y = Function('y')(x)   
    
    mode = 'add'
    #problem = "y.diff(x) - y**2 * exp(-2*x) = 0"
    #problem = "tan(x)*y.diff(x)-y=0"
    #problem = "y.diff(x)-3*x**2*y=0"
    #problem = "dy/dx-3*x**2*y=0"
    #if "*dy/dx" in problem:
    problem = entry.get()
    print(f"Your problem = {problem}")
    problem = problem.replace(' ','')
    problem = problem.replace('y\'','y.diff(x)')
    problem = problem.replace('dy/dx','y.diff(x)')
    if 'y.diff(x)*' in problem or '*y.diff(x)' in problem:
        mode = 'multiply'
    print(f'{mode} the dx/dy')
    problem = problem.replace('y.diff(x)','dy/dx')
    problem = problem.replace('dy/dx','0')
    problem = replacement(problem)
    problem = problem.replace('y','y(x)')
    print(f"Your rephrase problem = {problem}")

    # Initialize variables to store left-hand side and right-hand side of equations
    eqlhs, eqrhs = problem.split('=')
    print(f"before eqlhs = {eqlhs}, eqrhs = {eqrhs}")

    # Parse the user-provided equations
    eqlhs = sympify(eqlhs.strip())
        
    # Define the right-hand side separately
    eqrhs = sympify(eqrhs.strip())

    #if "*dy/dx" in problem:
        #equation = Eq(eqlhs, eqrhs)
    #else:
    if mode == 'multiply':
        equation = Eq(Derivative(y, x) * eqlhs, eqrhs)
        str_eqlhs = str(equation.lhs)
        if 'y(x)*' in str_eqlhs:
            str_eqlhs = str_eqlhs.replace('y(x)*','y(x)+')
            eqlhs = sympify(str_eqlhs)
            equation = Eq(eqlhs, eqrhs)
    else:
        equation = Eq(y.diff(x) + eqlhs, eqrhs)

    result_text = f"Solution:"
    result_text += f"\n\nEquation :\n{sp.pretty(equation)}"
    print(f"eqlhs = {eqlhs}, eqrhs = {eqrhs}")
    print("Original Equation:")
    print(equation)
    print(sp.pretty(equation))

    # Solve the differential equation
    solution = dsolve(equation)
    result_text += f"\n\nFinal answer:\n{sp.pretty(solution)}"

    y_solution = solution.rhs

    # Display the solution
    print("\nSolution:")
    print(f"y = {y_solution}")
    print(f"y = {sp.pretty(y_solution)}")
    output_label.config(text=result_text)
    final_label.config(text=sp.pretty(solution))


def Growth_And_Decay(user_input):
    # Remove spaces from user input
    user_input = user_input.replace(' ', '')

    # Split user input into a list of entries
    user_data = user_input.split(',')

    # Initialize variables
    p = None
    rp = None
    rt = None
    t = None

    # Extract values from user input
    for entry in user_data:
        variable_name, value = entry.split('=')

        try:
            if variable_name == 'p':
                p = float(value)
            elif variable_name == 'rp':
                rp = float(value)
            elif variable_name == 'rt':
                rt = float(value)
            elif variable_name == 't':
                t = float(value)
            else:
                output_label.config(text=f"\nInvalid variable name: {variable_name}\n")
        except ValueError:
            output_label.config(text=f"\nInvalid value for variable {variable_name}: {value}\n")

    # Check if all required variables are provided
    if p is None or rp is None or rt is None or t is None:
        result_text = "\nThere are some of them are None...\n"
        result_text += f"\nPopulation(p) = {p}"
        result_text += f"\nRate percent(rp) = {rp}"
        result_text += f"\nRate time(rt) = {rt}"
        result_text += f"\nTime(t) = {t}\n"
        result_text += f"\nPlease input something like this p=5000, rp=15, rt=10, t=30\n"
        output_label.config(text=result_text)
        return

    # Prepare result text
    result_text = "\nGiven:\n"
    result_text += f"\nPopulation = {int(p)}"
    result_text += f"\nRate percent = {int(rp)}"
    result_text += f"\nRate time = {int(rt)}"
    result_text += f"\nTime = {int(t)}\n"

    # Define the differential equation
    def dPdt(P, t):
        return k * P

    # Solve the differential equation using separation of variables
    def solve_ode(P0, k, t):
        return P0 * np.exp(k * t)

    # Calculate population
    P0 = p
    k = np.log(1 + (rp / 100)) / rt
    result_text += f"\nk = log(1+({int(rp)}/100))/{int(rt)}"
    result_text += f"\nk = {k:.6f}\n"

    P_t = solve_ode(P0, k, t)
    result_text += f"\nP_t = {int(P0)}*e^({k:.6f}*{int(t)})"
    result_text += f"\nP_t = {P_t:.2f}\n"

    # Calculate rate of population growth
    growth_rate = dPdt(P_t, t)
    result_text += f"\nGrowth rate = {k:.6f}*{int(P_t)}"
    result_text += f"\nGrowth rate = {growth_rate:.2f}\n"

    # Update result text with calculated values
    result_text += f"\nPopulation after {int(t)} years: {int(math.ceil(P_t))}"
    result_text += f"\nRate of population growth at t = {int(t)}: {int(math.ceil(growth_rate))}"

    # Display result text
    output_label.config(text=result_text)
    final_label.config(text=f"{int(math.ceil(P_t))} population\n{int(math.ceil(growth_rate))} population/yr")


def is_valid_sympy_expression(expression_str):
    expression_str = expression_str.replace('y\'','dy/dx')
    eqlhs, eqrhs = expression_str.split('=')
    try:
        sympify(eqlhs)
        sympify(eqrhs)
        return True
    except SympifyError:
        return False
    

#FRONT END, diri mo mag design2
def main():
    # diri mo himo na mura FRONT END or design
    root = tk.Tk()
    root.title('Differential Equation')
    root.geometry("720x480")
    root.resizable(False, False)

   # Load the background image and resize it to fit the window
    background_image = Image.open("design.png")
    background_image = background_image.resize((root.winfo_reqwidth()*4, root.winfo_reqheight()*5), 3)  # 3 corresponds to ANTIALIAS

    background_photo = ImageTk.PhotoImage(background_image)

    # Create a label to display the background image
    background_label = tk.Label(root, image=background_photo)
    background_label.place(relwidth=1, relheight=1.8)

    # Create a ComboBox
    style = ttk.Style()
    style.configure('TCombobox', relief='solid', highlightbackground='black', highlightcolor='black', highlightthickness=1, justify='center')
    combo = ttk.Combobox(root, values=('Derivative', 'Arbitrary Constant', 'Separable Variables', 'Growth and Decay'), style='TCombobox')
    combo.set('Derivative')
     
    # Set items for the ComboBox
    #combo['values'] = ('Derivative', 'Abritary Constant', 'Option 3', 'Option 4')

    # Set the default value (optional)
    #combo.set('Derivative')
    combo.set('Derivative')
    selected_item = combo.get()

    global dy_dx, eqlhs

    combo.pack(fill='x')
    combo.place(relx=0.5, rely=0.5, anchor='center', y=-150)  # Set the position to (10, 10)
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
        user_function = entry.get()
        elimination = entryEliminate.get()
        if user_function:
            if selected_item != "Growth and Decay" and '=' in user_function:
                if is_valid_sympy_expression(user_function):
                    if selected_item == "Derivative":
                        Derivatives()
                    elif selected_item == "Arbitrary Constant":
                        if elimination:
                            Arbitrary_Constant()
                        else:
                            result_text = f"\nYou didn't input any elimination.\n"
                            result_text += f"\nFor example, input k or any variables\n"
                            result_text += f"\nIf there's 2 elimination, just type \"C_1 C_2\" or any 2 variables\n"
                            output_label.config(text=result_text)
                    elif selected_item == "Separable Variables":
                        Separable_Variables()
                else:
                    result_text = f"\nThis is not a correct equation.\n"
                    result_text += f"\nThe correct way to enter an equation is like this:"
                    result_text += f"\ny=C_1*exp(2*x)+C_2*exp(-3*x)"
                    result_text += f"\ntan(x) * dy/dx - y = 0\n"
                    result_text += f"\nMake sure that there's multiplication (*) between terms.\n"
                    output_label.config(text=result_text)
            elif selected_item == "Growth and Decay":
                Growth_And_Decay(user_function)
            else:
                result_text = f"\nThere's no = sign.\n"
                output_label.config(text=result_text)
        else:
            result_text = f"\nYour input is empty.\n"
            # Handle the case where the input is empty
            output_label.config(text=result_text)

        # Add more conditions for other options if needed

    # Bind the function to the <<ComboboxSelected>> event
    combo.bind("<<ComboboxSelected>>", on_select)

    # Create a custom font
    custom_font = ('Arial', 32, 'bold italic')
    custom_font2 = tkFont.Font(family="Computer Modern", size=10)  # Set your desired font family and size

    # Create a label with the specified font and text
    #labelTitle = tk.Label(root, text="DIFFERENTIAL EQUATION", font=custom_font)
    #labelTitle.pack(fill='x')
    #labelTitle.place(relx=0.5, rely=0.5, anchor='center', y=-200)
    #labelTitle.configure(width=50)

    # Create an entry widget for user input
    global entry
    entry = tk.Entry(root, relief='solid', highlightthickness=1, justify='center')
    entry.pack(fill='x')
    #entry.place(relx=0.5, rely=0.5, anchor='center', y=-100)  
    #entry.configure(width=100)  
    entry.place(relx=0.5, rely=0.5, anchor='center', x=-190, y=-100)  
    entry.configure(width=50)  

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
    button.place(relx=0.5, rely=0.5, anchor='center', y=-50)
    button.configure(width=50)

    # Create a label to display the output
    global output_label
    output_label = tk.Label(root, text="", relief='flat', highlightthickness=1, font=custom_font2)
    output_label.pack(side='top', fill='x')
    output_label.pack(fill='x')
    output_label.place(relx=0.5, rely=0.5, anchor='center', y=90)
    output_label.configure(width=50)

    # Create a label to display the final answer
    global final_label
    final_label = tk.Label(root, text="", relief='flat', highlightthickness=1, font=custom_font2)
    final_label.pack(side='top', fill='x')
    final_label.pack(fill='x')
    final_label.place(relx=0.5, rely=0.5, anchor='center',x=190, y=-100)
    final_label.configure(width=35)

    root.mainloop()

if __name__ == '__main__':
    main()

