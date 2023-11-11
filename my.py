import tkinter as tk
from tkinter import Label, Entry, Button
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import sympy as sp

root = tk.Tk()
root.title ("Growth/Decay Calculator")

t, r, P0 = sp.symbols('t r P0')
P = sp.Function('P')(t)

def population_formula(P0, r, t):
    return P0 * np.exp(r * t)

def calculate_population():
    try:
        P0 = float(initial_entry.get())
        r = float(rate_entry.get())
        t = float(time_entry.get())

        P = population_formula(P0, r, t)

        result_label.config(text=f"Population at time {t}: {P:.2f}")
    except ValueError:
        result_label.config(text="Please enter valid numbers for the parameters.")


rate_label = Label(root, text="Rate (r):")
rate_label.pack()
rate_entry = Entry(root)
rate_entry.pack()

initial_label = Label(root, text="Initial Population (P0):")
initial_label.pack()
initial_entry = Entry(root)
initial_entry.pack()

time_label = Label(root, text="Time (t):")
time_label.pack()
time_entry = Entry(root)
time_entry.pack()

result_label = Label(root, text="")
result_label.pack()


calculate_button = Button(root, text="Calculate", command=calculate_population)
calculate_button.pack()

fig, ax = plt.subplots()
canvas = FigureCanvasTkAgg(fig, master=root)
canvas_widget = canvas.get_tk_widget()
canvas_widget.pack()



root.mainloop()