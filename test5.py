import sympy as sp

def findDVIV(problemText):
    IV = None
    DV = None
    numerator = True
    i = 0
    while i < len(problemText):
        if problemText[i] == 'd':  # If it detects 'd'
            i += 1
            if numerator:  # If this is NUMERATOR
                DV = ""
                while i < len(problemText) and problemText[i] not in ['+', '-', '*', '/', '=']:  # Repeat until it finds a symbol
                    DV += problemText[i]
                    i += 1
                if i < len(problemText) and problemText[i] == '/':
                    numerator = False  # After numerator, it will change to the denominator target
                else:
                    DV = None
            elif not numerator:  # If this is DENOMINATOR
                IV = ""
                while i < len(problemText) and problemText[i] not in ['+', '-', '*', '/', '=']:  # Repeat until it finds a symbol
                    IV += problemText[i]
                    i += 1
                numerator = True  # After denominator, it will change to the numerator target
        i += 1
    # Create strings for displaying the variables
    return DV, IV


while True:
  equation = input("Enter your equation: ")
  dv, iv = findDVIV(equation)
  print(f"DV = {dv}, IV = {iv}")
  try_again = input("Do you want to enter another equation? (y/n): ")
  if try_again.lower() != "y":
    break














