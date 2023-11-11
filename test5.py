import tkinter as tk
from tkinter import ttk
import numpy as np

def add_student():
    # Get the student name from the entry widget
    student_name = entry.get()

    # Add the student name to the list
    student_names.append(student_name)

    # Update the result label
    update_result_label()

def find_student():
    try:
        # Get the position from the entry widget
        position = int(entry_position.get())

        # Retrieve the student name based on the position
        student_name = student_names[position]

        # Display the result in the label
        result_label.config(text=f'Student at position {position}: {student_name}')
    except (ValueError, IndexError):
        result_label.config(text=f'Error: Invalid position')

def update_result_label():
    # Display the list of student names in the label
    result_label.config(text=f'Student Names: {student_names}')

# Create the main application window
root = tk.Tk()
root.title("Student Array Viewer")

# List to store student names
student_names = []

# Create an entry widget for adding students
entry = tk.Entry(root, width=30)
entry.pack(pady=5)

# Create a button to add students
add_button = tk.Button(root, text="Add Student", command=add_student)
add_button.pack(pady=5)

# Create a label to display the list of students
result_label = tk.Label(root, text="")
result_label.pack(pady=10)

# Entry widget for inputting the position
entry_position = tk.Entry(root, width=10)
entry_position.pack(pady=5)

# Button to find a student based on position
find_button = tk.Button(root, text="Find Student", command=find_student)
find_button.pack(pady=5)

# Start the Tkinter main event loop
root.mainloop()










