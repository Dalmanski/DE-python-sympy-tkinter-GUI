import tkinter as tk
from PIL import Image, ImageTk

label = None  # Declare label as a global variable

def load_image(root, current_index, image_paths):
    global label  # Access the global label variable
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()
    img_path = image_paths[current_index]
    img = Image.open(img_path).resize((screen_width + 385, screen_height + 200))
    photo_image = ImageTk.PhotoImage(img)
    label = tk.Label(root, image=photo_image)
    label.place(x=-2, y=-2)
    label.lower()
    label.image = photo_image  # Keep a reference to the image to prevent it from being garbage collected
    return label

def show_image(root, current_index, image_paths):
    global label  # Access the global label variable
    label.destroy()  # Destroy the existing label widget
    return load_image(root, current_index, image_paths)

def main():
    global label  # Access the global label variable
    screen_width = 800
    screen_height = 600
    image_paths = ["23232.jpg", "design2.jpg", "Untitled design1.jpg"]

    root = tk.Tk()
    root.title("Image Changer")

    button_1 = tk.Button(root, text="Image 1", command=lambda: show_image(root, 0, image_paths))
    button_2 = tk.Button(root, text="Image 2", command=lambda: show_image(root, 1, image_paths))
    button_3 = tk.Button(root, text="Image 3", command=lambda: show_image(root, 2, image_paths))

    button_1.pack(side=tk.LEFT)
    button_2.pack(side=tk.LEFT)
    button_3.pack(side=tk.LEFT)

    # Load the initial image
    label = load_image(root, 2, image_paths)
    show_image(root, 0, image_paths)

    # Run the main loop
    root.mainloop()

if __name__ == "__main__":
    main()
