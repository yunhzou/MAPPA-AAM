"""Create the README GIF from the rendered 3D film."""
from pathlib import Path
import argparse
import os
import shutil
import subprocess

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('directory', type=Path)
args = parser.parse_args()
ffmpeg = os.environ.get('FFMPEG') or shutil.which('ffmpeg')
if not ffmpeg:
    import imageio_ffmpeg
    ffmpeg = imageio_ffmpeg.get_ffmpeg_exe()
subprocess.run([ffmpeg, '-v', 'error', '-y', '-i', str(args.directory / 'gold-oxygen-3d.mp4'),
                '-filter_complex', 'fps=12,scale=960:-1:flags=lanczos,split[s0][s1];'
                '[s0]palettegen=max_colors=128[p];[s1][p]paletteuse=dither=bayer:bayer_scale=3',
                '-loop', '0', str(args.directory / 'gold-oxygen-preview.gif')], check=True)
