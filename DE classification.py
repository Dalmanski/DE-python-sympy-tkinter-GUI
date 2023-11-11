import sympy as sp

def classify_pde(pde):
    try:
        pde_eq = sp.Eq(sp.sympify(pde), 0)

        # Identify the dependent variable.
        dv = sp.Dummy()
        ivs = list(pde_eq.free_symbols - {dv})

        if len(ivs) >= 1:
            # The PDE is a partial differential equation (PDE).
            eqs = [eq.as_ordered_terms() for eq in sp.Eq(pde_eq.lhs, pde_eq.rhs).as_ordered_terms()]
            order = max([eq[0].as_poly(iv).degree() for eq in eqs])
            degree = sp.degree(pde_eq.rhs, dv)
            linear = all(len(term.as_ordered_terms()) <= 1 for term in pde_eq.rhs.as_ordered_terms())
            homogeneous = all(term == 0 for term in pde_eq.rhs.as_ordered_terms())

            classification = {
                "dv": dv,
                "ivs": ivs,
                "order": order,
                "degree": degree,
                "linear": linear,
                "homogeneous": homogeneous,
                "type": "PDE",
            }
            return classification

        # Return an error message if no independent variables found.
        return "Error: Invalid PDE (No independent variables found)"
    except Exception as e:
        # Return an error message for an invalid PDE.
        return "Error: Invalid PDE"

# Input the PDE interactively
pde_input = input("Enter a PDE in the format LHS = RHS: ")
classification_result = classify_pde(pde_input)

# Display the classification result
print("Classification:")
if isinstance(classification_result, dict):
    for key, value in classification_result.items():
        print(f"{key}: {value}")
else:
    print(classification_result)
