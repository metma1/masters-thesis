import glob
from PIL import Image

def make_gif(path, frame_folder):
    frame_path = path + "\\" + frame_folder
    frames = [Image.open(image) for image in glob.glob(f"{frame_path}\\*.png")]
    frame_one = frames[0]
    frame_one.save(frame_path + "\\" + frame_folder + ".gif", format="GIF", append_images=frames,
               save_all=True, duration=200, loop=0)