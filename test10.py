import numpy as np
import math

# Function to solve the differential equation
def solve_decay(P0, k, t):
    P_t = P0 * np.exp(k * t)  # Corrected the sign for decay
    growth_rate = k * P_t  # Corrected the sign for decay
    return P_t, growth_rate

# Function to calculate the half-life
def calculate_half_life(k):
    print(math.log(2))
    print(k)
    return math.log(2) / k

# Given information
p1 = 100  # Initial amount of radium in mg
p2 = 96  # Second population
t = 100  # Time in years

# Calculate the decay constant
k = np.log(p2 / p1) / t  # Negative sign for decay

# Time for 2.58 centuries
t_2_centuries = 2.58 * t

# Solve the differential equation
population, growth_rate = solve_decay(p1, k, t_2_centuries)

# Display the results
print(f"Population after {t_2_centuries:.2f} years: {population:.6f} mg")
print(f"Rate of population decay at t = {t_2_centuries:.2f} years: {growth_rate:.6f} mg/yr")
print(f"Decay constant k: {k:.6f} per year")
print(f"Half-life: {calculate_half_life(k):.2f} years")





