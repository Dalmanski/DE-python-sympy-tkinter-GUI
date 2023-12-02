import numpy as np

# Define the initial population and growth rate
P0 = 5000
k = np.log(1.15) / 10  # Growth rate per year
print(k)

# Define the differential equation
def dPdt(P, t):
    return k * P

# Solve the differential equation using separation of variables
def solve_ode(P0, k, t):
    return P0 * np.exp(k * t)

# Calculate the population after 30 years
t = 30
P_t = solve_ode(P0, k, t)

# Calculate the rate of population growth at t = 30
growth_rate = dPdt(P_t, t)

# Print the results
print(f"Population after 30 years: {P_t:.2f}")
print(f"Rate of population growth at t = 30: {growth_rate:.2f}")


