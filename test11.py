from customtkinter import CTk, CTkTextbox
import tkinter as tk

def create_rounded_rectangle(x1, y1, x2, y2, radius):
    points = [
        (x1 + radius, y1),
        (x2 - radius, y1),
        (x2, y1 + radius),
        (x2, y2 - radius),
        (x2 - radius, y2),
        (x1 + radius, y2),
        (x1, y2 - radius),
        (x1, y1 + radius),
    ]
    return tk.create_polygon(points, smooth=True)

root = CTk()
textbox = CTkTextbox(root)
textbox.pack()

canvas = textbox.get_canvas()
rounded_rect = create_rounded_rectangle(5, 5, 200, 30, 10)  # Adjust dimensions as needed
canvas.itemconfig(textbox.text_id, outline="", width=0)
canvas.tag_lower(rounded_rect, textbox.text_id)

root.mainloop()

