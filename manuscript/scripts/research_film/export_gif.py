"""Create a compact looping GIF from the research film, with a global palette."""
import argparse
import os
from pathlib import Path
import subprocess

ap = argparse.ArgumentParser()
ap.add_argument('directory', type=Path)
a = ap.parse_args()
encoder = os.environ.get('FFMPEG')
if not encoder:
    import imageio_ffmpeg
    encoder = imageio_ffmpeg.get_ffmpeg_exe()
movie = a.directory / 'mappa-grow-branch-decode.mp4'
palette = a.directory / 'palette.png'
subprocess.run([encoder, '-y', '-threads', '4', '-i', str(movie), '-vf',
                'fps=12,scale=960:-1:flags=lanczos,palettegen=max_colors=192',
                '-frames:v', '1', str(palette)], check=True)
subprocess.run([encoder, '-y', '-threads', '4', '-i', str(movie), '-i', str(palette),
                '-filter_complex_threads', '2', '-lavfi',
                'fps=12,scale=960:-1:flags=lanczos[x];[x][1:v]paletteuse=dither=bayer:bayer_scale=3',
                '-loop', '0', str(a.directory / 'mappa-grow-branch-decode.gif')], check=True)
palette.unlink()
