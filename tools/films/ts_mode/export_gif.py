"""Encode the GitHub preview and validate both movie formats."""
import json,os,subprocess,sys
from pathlib import Path
from PIL import Image
import imageio_ffmpeg
p=Path(sys.argv[1]);ffmpeg=os.environ.get('FFMPEG') or imageio_ffmpeg.get_ffmpeg_exe()
mp4=p/'rp-to-ts-mode.mp4';gif=p/'rp-to-ts-mode-preview.gif'
subprocess.run([ffmpeg,'-y','-i',str(mp4),'-vf','fps=12,scale=960:-1:flags=lanczos,split[s0][s1];[s0]palettegen=max_colors=128[p];[s1][p]paletteuse=dither=bayer:bayer_scale=3','-loop','0',str(gif)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
subprocess.run([ffmpeg,'-v','error','-i',str(mp4),'-f','null','-'],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
reader=imageio_ffmpeg.read_frames(str(mp4));meta=next(reader);count=sum(1 for _ in reader);assert count==720 and meta['size']==(1440,900)
im=Image.open(gif);total=0
for i in range(im.n_frames):im.seek(i);total+=im.info['duration']
assert im.n_frames==360 and total==30000 and im.size==(960,600)
(p/'media-validation.json').write_text(json.dumps(dict(status='passed',mp4_frames=count,mp4_size=meta['size'],mp4_fps=meta['fps'],gif_frames=im.n_frames,gif_size=im.size,gif_duration_ms=total),indent=2)+'\n')
print('Validated 720 MP4 frames and 360 GIF frames.')
