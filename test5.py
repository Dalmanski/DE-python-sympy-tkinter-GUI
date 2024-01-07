import tkinter as tk

def on_select(event):
    print("Item selected!")

def main():
    root = tk.Tk()
    root.title("Event Binding Example")

    label = tk.Label(root, text="Click the label!")
    label.pack(padx=20, pady=20)

    # Bind the left mouse button click event to on_select function
    label.bind("<Button-1>", on_select)

    root.mainloop()

if __name__ == '__main__':
    main()






















