from PIL import Image
import os

folder = "frames"

# Load all PNG/JPG files sorted alphabetically
files = sorted([f for f in os.listdir(folder) if f.lower().endswith((".png", ".jpg"))])

images = [Image.open(os.path.join(folder, f)) for f in files]

# Save as GIF (10 fps → 100 ms per frame)
images[0].save("output.gif",
               save_all=True,
               append_images=images[1:],
               duration=100,
               loop=0)