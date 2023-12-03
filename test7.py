import numpy as np
import math

def parse_user_input(user_input):
    user_data = user_input.split(',')

    p = None
    rp = None
    rt = None
    t = None

    for entry in user_data:
        variable_name, value = entry.split('=')

        if variable_name == 'p':
            p = value
        elif variable_name == 'rp':
            rp = value
        elif variable_name == 'rt':
            rt = value
        elif variable_name == 't':
            t = value
        else:
            print(f"variable {variable_name} is not define.")

    return p, rp, rt, t

user_input = "p1=5000, rp=0, rt=10, t=30, p2=96"
user_input = user_input.replace(' ','')
p1, rp, rt, t = parse_user_input(user_input)

print(f"p: {p1}")
print(f"rp: {rp}")
print(f"rt: {rt}")
print(f"t: {t}")

p1 = float(p1)
rp = float(rp)
rt = float(rt)
t = float(t)

# Define the differential equation
def dPdt(P1, t):
    return k * P1

# Solve the differential equation using separation of variables
def solve_ode(P0, k, t):
    return P0 * np.exp(k * t)

# If there's a percent (rp) on finding p3
# Calculate population
P0 = p1
k = np.log(1 + (rp / 100)) / rt
result_text = f"|k = log(1+({int(rp)}/100))/{int(rt)}"
result_text += f"|k = {k:.6f}|"
# If there's no percent (rp) and theres p2 on finding p3

P_t = solve_ode(P0, k, t)
result_text += f"|P_t = {int(P0)}*e^({k:.6f}*{int(t)})"
result_text += f"|P_t = {P_t:.2f}|"

# Calculate rate of population growth
growth_rate = dPdt(P_t, t)
result_text += f"|Growth rate = {k:.6f}*{int(P_t)}"
result_text += f"|Growth rate = {growth_rate:.2f}|"

# Update result text with calculated values
result_text += f"|Population after {int(t)} years: {int(math.ceil(P_t))}"
result_text += f"|Rate of population growth at t = {int(t)}: {int(math.ceil(growth_rate))}"

# Display result text
print(f"{int(math.ceil(P_t))} population\n{int(math.ceil(growth_rate))} population/yr")