import sympy as sp
import string
import re
import tkinter as tk
import math
import tkinter.font as tkFont
import numpy as np
from tkinter import ttk, font
from sympy import symbols, Function, dsolve, Eq, exp, simplify, sympify, tan, sin, cos, asin, acos, atan, Derivative, SympifyError
from PIL import Image, ImageTk
from scipy.optimize import fsolve


howManyDerive = 0
x, y, k, t, a, c, C_1, C_2, sin, cos, tan, asin, acos, atan, e, y_prime, y_prime2 = sp.symbols("x y k t a c C_1 C_2 sin cos tan asin acos atan e y\' y\'\'") 


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
    global x, y, k, t, a, c, C_1, C_2, sin, cos, tan, asin, acos, atan, e, y_prime, y_prime2

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
    user_function = entries[0].get()
    derive(user_function)
    # Display the result
    result_text = f"\nThe derivative is: \n{sp.pretty(eqlhs)} = {sp.pretty(dy_dx)}"
    final_label.insert(tk.END, f"{sp.pretty(eqlhs)} = {sp.pretty(dy_dx)}")
    return result_text


def Arbitrary_Constant():
    global eqlhs
    global equation1
    global howManyDerive
    howManyDerive = 0
    # Define the variable
    global x, y, k, t, a, c, C_1, C_2, sin, cos, tan, asin, acos, atan, e, y_prime, y_prime2
    elimination = entries[1].get()
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
            user_function = entries[0].get()
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
    result_text += f"\nThe Arbitrary Constant is: \n {sp.pretty(Answer)}"
    final_label.insert(tk.END, f"The Arbitrary Constant is: \n {sp.pretty(Answer)}")
    return result_text
    #output()


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


def Separable_Variables():
    # Define the variables
    y = Function('y')(x)   
    problem = entries[0].get()
    print(f"Your problem = {problem}")
    problem = problem.replace(' ','')
    problem = problem.replace('y\'','Derivative(y,x)')
    problem = problem.replace('dy/dx','Derivative(y,x)')
    eqlhs, eqrhs = problem.split('=')
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
    equation = Eq(eqlhs, eqrhs)

    result_text = f"\n<SOLUTION>\n"
    result_text += f"\nEquation :\n{sp.pretty(equation)}\n"
    print(f"eqlhs = {eqlhs}, eqrhs = {eqrhs}")
    print("Original Equation:")
    print(equation)
    print(sp.pretty(equation))
    # Solve the differential equation
    solution = dsolve(equation)
    result_text += f"\nFinal answer:\n{sp.pretty(solution)}"
    y_solution = solution.rhs
    # Display the solution
    print("\nSolution:")
    print(f"y = {y_solution}")
    print(f"y = {sp.pretty(y_solution)}")
    final_label.insert(tk.END, sp.pretty(solution))
    return result_text


def Growth_And_Decay(p1, p2, rp, rt, t1, t2, p3):
    find = None
    half_life = False

    if p1 == '?':
        result_text = "\nSorry, we didn't have a formula on finding p1"
        return result_text
    elif rp == '?':
        result_text = "\nSorry, we didn't have a formula on finding rp"
        return result_text
    elif rt == '?':
        result_text = "\nSorry, we didn't have a formula on finding rt"
        return result_text
    elif t2 == '?':
        result_text = "\nSorry, we didn't have a formula on finding t2"
        return result_text
    elif p2 == '?':
        result_text = "\nSorry, we didn't have a formula on finding t"
        return result_text
    elif p3 == '?':
        find = 'p3'
    else:
        result_text = "\nWhat should I find in? p1? rp? rt? t1? t2? p2? p3?"
        result_text += "\nIf you want to find p3 just input \"p3=?\"\n"
        return result_text
    
    if t1 == '?':
        half_life = True
 
    # Check if all required variables are provided
    if p1 is None or rp is None or rt is None or t is None or p2 is None or p3 is None:
        result_text = "\nThere are some of them are missing...\n"
        result_text += f"\nPopulation (p1) = {p1}"
        result_text += f"\nRate percent (rp) = {rp}"
        result_text += f"\nRate time (rt) = {rt}"
        result_text += f"\nTime (t2) = {t2}"
        result_text += f"\nPopulation (p2) = {p2}"
        result_text += f"\nPopulation at t={t2} (p3) = {p3}\n"
        result_text += f"\nPlease input something like this p1=5000, rp=15, rt=10, t2=30, p2=?, p3=?\n"
        return result_text

    p1 = float(p1)
    rp = float(rp)
    rt = float(rt)
    if t1 != '?':
        p3 = float(t1)
    t2 = float(t2)
    p2 = float(p2)
    if p3 != '?':
        p3 = float(p3)

    # Prepare result text
    result_text = "\n<SOLUTION>\n"
    result_text += "\nGiven:"
    result_text += f"\nPopulation (p1) = {p1}"
    result_text += f"\nRate percent (rp) = {rp}"
    result_text += f"\nRate time (rt) = {rt}"
    result_text += f"\nTime (t1) = {t1}"
    result_text += f"\nTime (t2) = {t2}"
    result_text += f"\nPopulation (p2) = {p2}"
    result_text += f"\nPopulation at t={t2} (p3) = {p3}\n"

    # Solve the differential equation using separation of variables
    def solve_ode(P0, k, t):
        return P0 * np.exp(k * t)

    # Calculate population
    if int(rp) != 0:
        P0 = p1
        k = np.log(1 + (rp / 100)) / rt
        result_text += f"\nk = log(1+({int(rp)}/100))/{int(rt)}"
        result_text += f"\nk = {k:.6f}\n"
        # If there's no percent (rp) and theres p2 on finding p3
        # Calculate the decay constant
    else:
        k = np.log(p2 / p1) / rt 
        result_text += f"\nk = log({int(p2)}/{int(p1)}))/{int(rt)}"
        result_text += f"\nk = {k:.6f}\n"

    P_t = solve_ode(p1, k, t2)
    result_text += f"\nP_t = {int(p2)}*e^({k:.6f}*{int(t2)})"
    result_text += f"\nP_t = {P_t:.2f}\n"

    # Calculate rate of population growth
    growth_rate = P_t * k
    result_text += f"\nGrowth rate = {k:.6f}*{int(P_t)}"
    result_text += f"\nGrowth rate = {growth_rate:.2f}\n"

    # Update result text with calculated values
    result_text += f"\nPopulation after {int(t2)} years: {int(math.ceil(P_t))}"
    result_text += f"\nRate of population growth at t = {int(t2)}: {int(math.ceil(growth_rate))}"

    if half_life == True:
        t1 = math.log(2) / k
        result_text += f"\nHalf-life: t1 = {t1:.2f}"
        final_label.insert(tk.END, f"{int(math.ceil(P_t))} population\nhalf-life (t1) = {abs(t1):.2f} years")
    else:
        # Display result text
        final_label.insert(tk.END, f"{int(math.ceil(P_t))} population\n{int(math.ceil(growth_rate))} population/yr")
    return result_text


def Newtons_Law_of_Cooling_Heating(T1, Tm, T2, t1, t2, T3):  
    find = None      
    if T1 == '?':
        result_text = "\nSorry, we didn't have a formula on finding T1"
        return result_text
    elif Tm == '?':
        result_text = "\nSorry, we didn't have a formula on finding Tm"
        return result_text
    elif T2 == '?':
        result_text = "\nSorry, we didn't have a formula on finding T2"
        return result_text
    elif t1 == '?':
        result_text = "\nSorry, we didn't have a formula on finding t1"
        return result_text
    elif t2 == '?':
        find = 't2'
    elif T3 == '?':
        find = 'T3'
    else:
        result_text = "\nWhat should I find in? T1? Tm? T2? t1? t2? T3?"
        result_text += "\nIf you want to find t1 just input \"t1=?\""
        return result_text

    # Check if all required variables are provided
    if T1 is None or Tm is None or T2 is None or t1 is None or t2 is None or T3 is None:
        result_text = "\nThere are some of them are missing...\n"
        result_text += f"\nInitial Temperature (T1) = {T1}"
        result_text += f"\nRoom Temperature (Tm) = {Tm}"
        result_text += f"\nAfter Temperature (T2) = {T2}"
        result_text += f"\nMinute later (t1) = {t1}"
        result_text += f"\nTime after (t2) = {t2}\n"
        result_text += f"\nPlease input like this T1=18, Tm=70, T2=31, t1=1, t2=5, t3=?"
        return result_text
    
    T1 = float(T1) 
    Tm = float(Tm)
    T2 = float(T2)
    t1 = float(t1)
    if t2 != '?':
        t2 = float(t2)
    if T3 != '?':
        T3 = float(T3)

    # Prepare result text
    result_text = "\n<SOLUTION>\n"
    result_text += "\nGiven:"
    result_text += f"\nInitial Temperature (T1) = {T1}"
    result_text += f"\nRoom Temperature (Tm) = {Tm}"
    result_text += f"\nAfter Temperature (T2) = {T2}"
    result_text += f"\nMinute later (t1) = {t1}"
    result_text += f"\nTime after (t2) = {t2}"
    result_text += f"\nTemperature after t = {t2} (T3) = {T3}\n"


    def temperature_reading(initial_temperature, room_temperature, t, rate_constant):
        # Calculate the temperature reading at time t
        temperature = room_temperature + (initial_temperature - room_temperature) * np.exp(-rate_constant * t)   
        return temperature

    # Calculate the rate constant
    # Define the constants
    equation = lambda rate_constant: Tm + (T1 - Tm) * np.exp(-rate_constant * t1)
    # Solve for the rate constant using the temperature_after_min value
    rate_constant = fsolve(lambda x: equation(x) - T2, 0.2)[0]
    result_text += f"\n{int(equation(rate_constant))} = {int(Tm)} + {int(Tm-T1)}e^k({int(t1)})"
    result_text += f"\nk = {rate_constant}"

    if find == 'T3':
        # Calculate the temperature reading after minutes
        t2 = t2  # Time in minutes
        temperature_after_t2 = temperature_reading(T1, Tm, t2, rate_constant)
        result_text += f"\n{int(Tm)} + ({int(T1)} - {int(Tm)} * e^(-{rate_constant} * {int(t1)}))\n"
        result_text += f"\nTemperature reading after t = {int(t2)}: {temperature_after_t2:.2f}°\n"
        final_label.insert(tk.END, f"Temperature(T3) = {temperature_after_t2:.2f}°")
    elif find == 't2':
        def find_time_for_temperature(target_temperature, rate_constant):
            # Define the equation for temperature at time t
            equation = lambda t: temperature_reading(T1, Tm, t, rate_constant) - target_temperature
            # Solve for time using fsolve
            time = fsolve(equation, 0.2)[0]
            return time
        # Find the time for a target temperature
        target_temperature = T3
        time_for_target_temperature = find_time_for_temperature(target_temperature, rate_constant)
        result_text += f"\n{int(Tm)} + ({int(T1)} - {int(Tm)} * e^(-{rate_constant} * {int(t1)}))\n"
        result_text += f"\nTime for target temperature ({target_temperature}°F): {math.ceil(time_for_target_temperature)} minutes\n"
        final_label.insert(tk.END, f"Time(t2) = {math.ceil(time_for_target_temperature)} minutes")
    
    #output_label.config(text=result_text)
    return result_text
  


#FRONT END, diri mo mag design2
def main():
    root = tk.Tk()
    root.title('Differential Equation')
    root.geometry("720x480")
    root.resizable(False, False)

    background_image = Image.open("design.png")
    background_image = background_image.resize((root.winfo_reqwidth() * 4, root.winfo_reqheight() * 5), 3)
    background_photo = ImageTk.PhotoImage(background_image)

    background_label = tk.Label(root, image=background_photo)
    background_label.place(relwidth=1, relheight=1.8)

    # Create a ComboBox
    style = ttk.Style()
    style.configure('TCombobox', relief='solid', highlightbackground='black', highlightcolor='black', highlightthickness=1, justify='center')
    combo = ttk.Combobox(root, style='TCombobox')
    #combo.set('Newton\'s Law of Cooling/Heating')   
    combo.set('Derivative')   
    combo['values'] = ('Derivative', 'Arbitrary Constant', 'Separable Variables', 'Growth and Decay', 'Newton\'s Law of Cooling/Heating')
    selected_item = combo.get()

    global dy_dx, eqlhs
    global entries

    combo.pack(fill='x')
    combo.place(relx=0.5, rely=0.5, anchor='center', y=-170)  # Set the position to (10, 10)
    combo.configure(width=50)  # Set the width to 100


    def on_select(event):
        global entries
        topic = combo.get()
        listbox.config(state=tk.NORMAL)
        entryText.config(state=tk.NORMAL)
        clear_listbox(listbox)
        clear_listbox(entryText)
        clear_listbox(final_label)

        if topic == "Derivative":
            result_text = "\nThis is Derivative.\n"
            result_text += "\nHere are the examples:"
            result_text += "\ny=a*x^2+b*x+c"
            result_text += "\ny=6*x^3-9*x+4"
            insert_multiline_text(listbox, result_text)
            entries = entryShow(1)
        elif topic == f"Arbitrary Constant":
            result_text = "\nThis is Arbitrary Constant. Make sure you input"
            result_text += "\non the input of elimination.\n"
            result_text += "\nExample:"
            result_text += "\ny=C_1*e^(2*x)+C_2*e^(-3*x)"
            result_text += "\nOn the 2nd input (Eliminate)"
            result_text += "\nC_1 C_2"
            insert_multiline_text(listbox, result_text)
            entries = entryShow(2)
        elif topic == "Separable Variables":
            insert_multiline_text(entryText, "INPUT:")
            result_text = "\nThis is Separable Variables.\n"
            result_text += "\nExample:"
            result_text += "\ntan(x) * dy/dx - y = 0"
            result_text += "\ndy/dx - y**2 * exp(-2*x) = 0"
            result_text += "\ndy/dx - 3*x^2*y = 0"
            result_text += "\ny'=2*x*y+3*y-4*x-6"
            result_text += "\ndy/dx = f(x)/g(y)"
            insert_multiline_text(listbox, result_text)
            entries = entryShow(1)
        elif topic == "Growth and Decay":
            insert_multiline_text(entryText, " p1:  p2:  rp:  rt:  t1:  t2:  p3:  ")
            result_text = "\nThis is Growth and Decay.\n"
            result_text += "\nInput example:"
            result_text += "\np1=5000, p2=5000, rp=15, rt=10, t1=0, t2=30, p3=?"
            result_text += "\np1=100,  p2=96,  rp=0,  rt=100, t1=0, t2=258,  p3=?"
            result_text += "\np1=100,  p2=96,  rp=0,  rt=100, t1=?, t2=258,  p3=?"
            insert_multiline_text(listbox, result_text)
            entries = entryShow(7)
        elif topic == "Newton\'s Law of Cooling/Heating":
            insert_multiline_text(entryText, " T1:  Tm:  T2:  t1:  t2:  T3:  ")
            result_text = "\nThis is Newton\'s Law of Cooling/Heating\n"
            result_text += "\nInput example:"
            result_text += "\nT1=18, Tm=70, T2=31, t1=1, t2=5, T3=?"
            result_text += "\nT1=20, T2=22, T3=90, t1=1, t2=?, Tm=100"
            result_text += "\nT1=20, T2=22, T3=98, t1=1, t2=?, Tm=100"
            insert_multiline_text(listbox, result_text)
            entries = entryShow(6)

        #print(f"Selected option: {topic}")
        listbox.config(state=tk.DISABLED)
        entryText.config(state=tk.DISABLED)
    
    # Bind the function to the <<ComboboxSelected>> event
    combo.bind("<<ComboboxSelected>>", on_select)

    global final_label
    final_label = tk.Text(root, relief='flat', highlightthickness=1, height=2)
    final_label.place(relx=0.5, rely=0.5, anchor='center',x=195, y=-97)
    final_label.config(bg="black")
    final_label.config(fg="white")
    final_label.configure(width=35)

    def insert_multiline_text(text_widget, text):
        lines = text.split('\n')
        #num_lines = int(final_label.index('end-1c').split('.')[0])
        #final_label.config(height=1, y=-100)
        #height_increment = 0
        #y_increment = 0
        #for i in range(num_lines):
        #    if i > 1:
        #        height_increment += 1
        #        y_increment += 10
        #        final_label.config(height=1 + height_increment, y=-100 + y_increment)
        for line in lines:
            text_widget.insert(tk.END, line + '\n')


    def clear_listbox(listbox):
        listbox.delete(1.0, tk.END)

    def is_valid_sympy_expression(expression_str):
        if not ',' in expression_str:
            expression_str = expression_str.replace('y\'','dy/dx')
            eqlhs, eqrhs = expression_str.split('=')
            try:
                sympify(eqlhs)
                sympify(eqrhs)
                return True
            except SympifyError:
                return False
        else:
            return False


    def button_clicked():
        global entries
        selected_item = combo.get()
        user_function = entries[0].get()
        listbox.config(state=tk.NORMAL)
        final_label.config(state=tk.NORMAL)
        clear_listbox(listbox)
        clear_listbox(final_label)
        if user_function:
            if selected_item != "Growth and Decay" and selected_item != "Newton\'s Law of Cooling/Heating" and '=' in user_function:
                if is_valid_sympy_expression(user_function):
                    if selected_item == "Derivative":
                        result_text = Derivatives()
                        insert_multiline_text(listbox, result_text)
                    elif selected_item == "Arbitrary Constant":
                        elimination = entries[1].get()
                        if elimination:
                            result_text = Arbitrary_Constant()
                            insert_multiline_text(listbox, result_text)
                        else:
                            result_text = "\nYou didn't input any elimination.\n"
                            result_text += "\nFor example, input k or any variables\n"
                            result_text += "\nIf there's 2 elimination, just type \"C_1 C_2\" or any 2 variables\n"
                            insert_multiline_text(listbox, result_text)
                    elif selected_item == "Separable Variables":
                        result_text = Separable_Variables()
                        insert_multiline_text(listbox, result_text)
                    else:
                        result_text = "\nInvalid...\n"
                        insert_multiline_text(listbox, result_text)
                elif ',' in user_function:
                    result_text = "\nThis topic doesn't require comma(,)\n"
                    result_text += "\nSwitch the topic to growth and decay or newtons law of cooling/heating."
                    insert_multiline_text(listbox, result_text)
                else:
                    result_text = "\nThis is not a correct equation.\n"
                    result_text += "\nThe correct way to enter an equation is like this:"
                    result_text += "\ny=C_1*exp(2*x)+C_2*exp(-3*x)"
                    result_text += "\ntan(x) * dy/dx - y = 0\n"
                    result_text += "\nMake sure that there's multiplication (*) between terms.\n"
                    insert_multiline_text(listbox, result_text)
            elif selected_item == "Growth and Decay":
                inp1 = entries[0].get()
                inp2 = entries[1].get()
                inp3 = entries[2].get()
                inp4 = entries[3].get()
                inp5 = entries[4].get()
                inp6 = entries[5].get()
                inp7 = entries[6].get()
                result_text = Growth_And_Decay(inp1, inp2, inp3, inp4, inp5, inp6, inp7)
                insert_multiline_text(listbox, result_text)
            elif selected_item == "Newton\'s Law of Cooling/Heating":
                inp1 = entries[0].get()
                inp2 = entries[1].get()
                inp3 = entries[2].get()
                inp4 = entries[3].get()
                inp5 = entries[4].get()
                inp6 = entries[5].get()
                result_text = Newtons_Law_of_Cooling_Heating(inp1, inp2, inp3, inp4, inp5, inp6)
                insert_multiline_text(listbox, result_text)
            else:
                result_text = f"\nThere's no = sign.\n"
                insert_multiline_text(listbox, result_text)
        else:
            result_text = f"\nYour input is empty.\n"
            # Handle the case where the input is empty
            insert_multiline_text(listbox, result_text)

        listbox.config(state=tk.DISABLED)
        final_label.config(state=tk.DISABLED)


    def create_entry(root, x_offset):
        entry = tk.Entry(root, relief='flat', highlightthickness=1, justify='center', bg="black", fg="white", width=5)
        # entry.insert(tk.END, "hi")  # Insert "hi" into the entry
        entry.pack(fill='x')
        entry.place(relx=0.5, rely=0.5, anchor='center', x=x_offset-310, y=-100)
        return entry

    def entryShow(numEntry):
    # Destroy existing entries, excluding the ComboBox
        for widget in root.winfo_children():
            if isinstance(widget, tk.Entry) and widget.winfo_name() != combo.winfo_name():
                widget.destroy()
        entries = []
        if numEntry == 1:
            entry = tk.Entry(root, relief='flat', highlightthickness=1, justify='center')
            entry.pack(fill='x')
            entry.place(relx=0.5, rely=0.5, anchor='center', x=-190, y=-100)
            entry.config(bg="black", fg="white")
            entry.configure(width=48)
            entries.append(entry)
        else:
            # Create and configure entries dynamically based on numEntry
            entries = [create_entry(root, 40 * i) for i in range(numEntry)]
        return entries

    # Call entryShow with the desired number of entries to show 
    entries = entryShow(1)

    entryText = tk.Text(root, relief='flat', highlightthickness=0, wrap='none', height=1)
    entryText.pack(fill='x')
    entryText.place(relx=0.5, rely=0.5, anchor='center', x=-190, y=-135)
    entryText.config(bg="white", fg="black")
    # Create a Font object and configure it to be bold
    bold_font = font.Font(entryText, entryText.cget("font"))
    bold_font.configure(weight="bold")
    entryText.configure(font=bold_font, width=36)
    result_text = "INPUT:"
    insert_multiline_text(entryText, result_text)

    # Create a button to trigger the action
    button = tk.Button(root, text="Compute", command=button_clicked, relief='solid', highlightthickness=1)
    button.place(relx=0.5, rely=0.5, anchor='center', y=-45)
    button.configure(width=50)

    # Define custom font
    listbox = tk.Text(root, bg="white", wrap=tk.WORD, width=55, height=11)
    listbox.place(relx=0.5, rely=0.5, anchor='center', y=105)

    result_text = "\nWelcome to Differential Equation.\n"
    result_text += "\nPick a topic."
    insert_multiline_text(listbox, result_text)

    # Disable text editing
    listbox.config(state=tk.DISABLED)
    entryText.config(state=tk.DISABLED)
    
    # Create a label to display the final answer
    

    root.mainloop()

if __name__ == '__main__':
    main()