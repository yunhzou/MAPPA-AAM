"""Export a GitHub-renderable preview and a still from the verification MP4."""
from pathlib import Path
import argparse,os,shutil,subprocess
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('directory',type=Path)
a=parser.parse_args()
ffmpeg=os.environ.get('FFMPEG') or shutil.which('ffmpeg')
if not ffmpeg:
    import imageio_ffmpeg
    ffmpeg=imageio_ffmpeg.get_ffmpeg_exe()
subprocess.run([ffmpeg,'-v','error','-y','-i',str(a.directory/'broken-molecule-verification.mp4'),
    '-filter_complex','fps=12,scale=960:-1:flags=lanczos,split[s0][s1];[s0]palettegen=max_colors=128[p];[s1][p]paletteuse=dither=bayer:bayer_scale=3',
    '-loop','0',str(a.directory/'broken-molecule-verification-preview.gif')],check=True)
subprocess.run([ffmpeg,'-v','error','-y','-ss','26','-i',str(a.directory/'broken-molecule-verification.mp4'),
    '-frames:v','1',str(a.directory/'poster.png')],check=True)
