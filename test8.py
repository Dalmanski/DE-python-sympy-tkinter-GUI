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

# Calculate the temperature reading after 5 minutes
t = 5  # Time in minutes
temperature_after_5min = temperature_reading(t, rate_constant)
print("Temperature reading after 5 minutes:", temperature_after_5min, "°F")








