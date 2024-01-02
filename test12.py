def solving():
    problemText = ""
    IV = []  # List for Independent Variables
    DV = []  # List for Dependent Variables
    partial = '\u2202'  # For ∂
    numerator = True
    theresPartial = False
    parenthesis = False
    reciprocalDVIV = False
    order = '1'
    degree = '0'
    preExponential = '0'  # it is used to check if the denominator exponent is the same on the numerator exponent
    partExponential = '0'  # For solving each part of their exponent
    linearity = ""

    problemText = input("Enter the problem equation: ")  # Get the problem equation
    numerator = True

    print("Solution:")

    DV.clear()
    IV.clear()

    selectComboBox = "Classification by Type"

    if selectComboBox == "Classification by Type":
        i = 0
        while i < len(problemText):
            if problemText[i] == '(':  # Check if there's a parenthesis
                parenthesis = True
            elif problemText[i] == ')':
                parenthesis = False
                if i + 1 < len(problemText) and problemText[i + 1] == '^':
                    i += 2  # Skip '^' and move to the exponent character
                    partExponential = problemText[i]
                else:
                    partExponential = '1'  # If no exponent, assume it's '1'

            if problemText[i] == 'd' or problemText[i] == partial:  # If it detects 'd' or '∂'
                if problemText[i] == partial:  # If this is partial
                    theresPartial = True

                i += 1
                if problemText[i] == '^':  # If it is exponential, get the Order if it is the highest
                    i += 1
                    order = problemText[i]
                    if numerator:  # Get the exponential from the numerator
                        preExponential = problemText[i]
                elif numerator:  # If it's not exponent
                    preExponential = '1'
                    linearity = "Linear"

                while i < len(problemText) and problemText[i] not in ['x', 'y', 'z', 'v']:
                    i += 1  # Find until 'x', 'y', 'z', or 'v'

                if numerator:  # If this is NUMERATOR
                    if problemText[i] not in DV:  # To prevent repeat answer
                        DV.append(problemText[i])

                    while i < len(problemText) and problemText[i] != '/':  # Repeat until it finds a fraction or not
                        if problemText[i] in ['+', '-']:  # If it's found '+' or '-' instead of '/'
                            numerator = False
                            break
                        i += 1

                    if i < len(problemText) and problemText[i] == '/':
                        numerator = False  # After numerator, it will change to the denominator target
                    else:  # If didn't detect fraction, it will reverse DV and IV
                        reciprocalDVIV = True
                elif not numerator:  # If this is DENOMINATOR
                    if problemText[i] not in IV:  # To prevent repeat answer
                        IV.append(problemText[i])

                    if degree == '1':  # Check if it is linear or not linear (Work on the 2nd)
                        linearity = "Linear"

                    if order == preExponential:  # Check if the exponent of the numerator is the same as the exponent of the denominator
                        degree = '1'

                    if parenthesis:  # Get the exponential on that part
                        partExponential = preExponential

                    numerator = True  # After denominator, it will change to the numerator target
            i += 1

        # Create strings for displaying the variables
        dependentVariable = ",".join(DV)
        independentVariable = ",".join(IV)

        if not reciprocalDVIV:
            print("\nDependent Variable:", dependentVariable)
            print("Independent Variable:", independentVariable)
        else:
            print("\nDependent Variable:", independentVariable)
            print("Independent Variable:", dependentVariable)

        print("Type:", "Partial" if len(DV) > 1 or len(IV) > 1 or theresPartial else "Ordinary")
        print("Order:", order)
        print("Degree:", degree)
        print("Linearity:", linearity)
    else:
        print("Please Select the Formula.")


# Call the function
solving()








