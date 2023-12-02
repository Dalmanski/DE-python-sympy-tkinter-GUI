import tkinter as tk

root = tk.Tk()

# Create a Frame to hold the label and scrollbar
frame = tk.Frame(root)
frame.pack()

# Create the Label
label = tk.Label(frame, text="This is some text that will be scrolled down", wraplength=300)
label.pack()

# Create the Scrollbar
scrollbar = tk.Scrollbar(frame, orient=tk.VERTICAL, command=label.yview)
scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

# Configure the Label to expand and fill the available space
label.config(height=10)

# Configure the Scrollbar to match the height of the Label
label['yscrollcommand'] = scrollbar.set

root.mainloop()

