import numpy as np
from scipy.optimize import fsolve

def calculate_rate_constant(initial_temperature, room_temperature, temperature_after_1min, time):
    # Define the constants
    equation = lambda rate_constant: room_temperature + (initial_temperature - room_temperature) * np.exp(-rate_constant * time)

    # Solve for the rate constant using the temperature_after_1min value
    rate_constant = fsolve(lambda x: equation(x) - temperature_after_1min, 0.2)[0]

    return rate_constant

def temperature_reading(t, rate_constant):
    # Define the constants
    initial_temperature = 18  # Initial temperature in °F
    room_temperature = 70  # Room temperature in °F

    # Calculate the temperature reading at time t
    temperature = room_temperature + (initial_temperature - room_temperature) * np.exp(-rate_constant * t)

    return temperature


# Calculate the rate constant
rate_constant = calculate_rate_constant(18, 70, 31, 1)
print("Rate constant:", rate_constant)

def find_time_for_temperature(target_temperature, rate_constant):
    # Define the equation for temperature at time t
    equation = lambda t: temperature_reading(t, rate_constant) - target_temperature

    # Solve for time using fsolve
    time = fsolve(equation, 0.2)[0]

    return time

# Find the time for a target temperature
target_temperature = 57.66
time_for_target_temperature = find_time_for_temperature(target_temperature, rate_constant)
print(f"Time for target temperature ({target_temperature}°F): {time_for_target_temperature} minutes")
