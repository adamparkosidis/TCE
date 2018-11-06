from PIL import Image, ImageDraw, ImageFont
import os
import sys
import glob
import re
from images2gif import writeGif

font = ImageFont.truetype("arial.ttf", 16)

def animate_files(files, output_path, begin, end, duration=0.1):
    " Convert the list of files to an animated GIF "
    # Create a list of PIL images
    images = []
    for f in [files[i] for i in range(int(begin),int(end),1)]:
        print f
        img = Image.open(f)
        draw = ImageDraw.Draw(img)
        draw.text((0, 0),f,(0,0,0),font=font)
        images.append(img)
    # Write the file
    writeGif(output_path,
             images,
             duration=duration)
def main(args=None):
    print args
    if len(args) > 1:
        os.chdir(args[1])
    if len(args) > 2:
        begin = int(args[2])
    else:
        begin = 0
    if len(args) > 3:
        end = int(args[3])
    else:
        end = 0
    if len(args) > 4:
        adding = args[4] + '_'
    else:
        adding = ""
    files = glob.glob('./' + adding + 'plotting_*.jpg')
    files.sort(key=lambda f: int(re.findall(r'\d+', os.path.basename(f))[0] ))
    if end == 0:
        end = len(files)
    print "beginning animation from step = ", begin, "to step = ", end
    animate_files(files, './animated.gif', begin, end)
    print "animation have been created"

if __name__=='__main__':
    main(sys.argv)
