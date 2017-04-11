import numpy as np
import PIL
import PIL.Image

def visual_histogram(img, lowval=None, highval=None):
    #input: image (preferably greyscale)
    #output: visual histogram as image, binary image, 256x100
    if lowval==None:
        lowval=0
    if highval==None:
        highval=255
    hist = img.histogram()
    if len(hist)!=256:
        print "Convert to greyscale first"
        return
    maxval = max(hist)
    normheight = [int(100*i/maxval) for i in hist]
    histimg = PIL.Image.new(mode="RGB", size=(256,100))
    draw = PIL.ImageDraw.Draw(histimg)
    for i in range(256):
        draw.line((i,100) + (i,100-normheight[i]), 
                fill=(255,255*(i>lowval)*(i<highval)*(highval>lowval),255*(i>lowval)*(i<highval)))
    del draw
    return histimg

def unit_disk(n):
    #Make this faster later if necessary, probably good enough for now
    h = (n-1)/2.
    mat = np.ones((n,n))
    for i in range(n):
        for j in range(n):
            if (i-h)**2+(j-h)**2>h**2+.5:
                mat[i,j]=0
    return mat
