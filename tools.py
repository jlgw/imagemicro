import numpy as np
import PIL
import PIL.Image
import PIL.ImageChops
import scipy

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

def watershed_pts(img, point_list):
    #input: image, list of points
    #output: watershed image where each point is a basin marker, each basin is
    #given a unique shade
    #Note: point_list should be of length < 256
    markers = np.zeros(img.size).astype(np.int16)
    factor = 256/len(point_list) # just for looks and easier use
    for i,a in enumerate(point_list):
        markers[a[0],a[1]] = i*factor
    imin = np.array(img).transpose()
    imgdata = scipy.ndimage.measurements.watershed_ift(np.array(imin), markers)
    img2 = PIL.Image.new(mode="L", size=img.size)
    img2.putdata(imgdata.transpose().flatten())
    return img2

def grid(img):
    #input: image
    #output: list of points, making up an NxN grid uniform with respect
    #to the x and y axes
    xpts = [int((0.5+i)*self.img.size[0]/n) for i in range(n)]
    ypts = [int((0.5+i)*self.img.size[1]/n) for i in range(n)]
    for i,a in enumerate(xpts):
        for j,b in enumerate(ypts):
            points[i*n+j] = (a,b)
    return points
