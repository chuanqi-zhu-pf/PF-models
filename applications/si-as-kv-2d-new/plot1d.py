import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from moviepy.editor import ImageSequenceClip
# from moviepy.editor import VideoFileClip
# clip = VideoFileClip("my_video.mp4", ffmpeg_exec="path/to/ffmpeg")


ns = 1000
step_arr = np.arange(0, ns*61, ns)
# frames = []
for step in step_arr:
    dfc = pd.read_csv(f"data/con/1d{step}.csv", header=None)
    arrc = np.array(dfc[0].values)
    plt.plot(arrc)
    # plt.ylim([0.06, 0.22])
    # plt.text(0.5, 0.5, "text")
    plt.savefig(f"fig/con/1d{step}.png")
    # Clear the figure for the next frame
    # plt.clf()

    # # Append the image to the list of frames
    # frames.append(f'fig/con/1d{step}.png')
    plt.close()

# clip = ImageSequenceClip(frames, fps=10)
# clip.write_videofile('output.mp4')

# dfc = pd.read_csv(f"data/phi/1d{step}.csv", header=None)
# arrc = dfc[0].values
# plt.plot(arrc)
# plt.savefig(f"fig/phi/1d{step}")
# plt.close()

# dfc = pd.read_csv(f"data/conl/1d{step}.csv", header=None)
# arrc = dfc[0].values
# plt.plot(arrc)
# plt.savefig(f"fig/conl/1d{step}")
# plt.close()

# dfc = pd.read_csv(f"data/temp/1d{step}.csv", header=None)
# arrc = dfc[0].values
# plt.plot(arrc)
# plt.savefig(f"fig/temp/1d{step}")
# plt.close()
