import numpy as np
import PIL
import PIL.Image
import PIL.ImageChops
import scipy
import scipy.ndimage

def visual_histogram(img, lowval=None, highval=None):
    #input: image (preferably greyscale)
    #output: visual histogram as image, binary image, 256x100
    if lowval==None:
        lowval=0
    if highval==None:
        highval=255
    hist = img.histogram()
    if len(hist)!=256:
        print "#Non-greyscale image in display buffer, is histogram window already open?"
        return
    maxval = max(hist)
    normheight = [int(100*i/maxval) for i in hist]
    histimg = PIL.Image.new(mode="RGB", size=(256,100))
    draw = PIL.ImageDraw.Draw(histimg)
    for i in range(256):
        draw.line((i,100) + (i,100-normheight[i]), 
                fill=(255,255*(1-(i>lowval)*(i<highval)*(highval>lowval)),
                    255*(1-(i>lowval)*(i<highval))))
    del draw
    return histimg

def aver(hist, length):
    avg = len(hist)*[0]
    avval = 0
    for i in range(length, len(hist)):
        avval += hist[i]-hist[i-length]
        avg[i-length/2] = avval/length
    return avg

def highest_2nd_derivative(img):
    # Input: List of length 256
    # Output: Point of maximum approximated second derivative
    deg = 4
    d = 2
    padlength = 20
    hist = img.histogram()
    p = np.polyfit(range(256), hist, deg)
    pd = [p[n]*reduce(lambda x,y: x*y, range(deg-n-d+1,deg-n+1)) for n in range(deg-d+1)]

    leftpad = padlength*[0]
    rightpad = padlength*[0]
    pl = leftpad + [sum([pd[deg-d-n]*i**(n) for n in range(0, deg+1-d)]) 
            for i in range(256)] + rightpad
    pl = aver(pl,40)[padlength:256+padlength]
    return pl.index(max(pl[padlength:255+padlength]))



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

def grid(img, n):
    #input: image
    #output: list of points, making up an NxN grid uniform with respect
    #to the x and y axes
    points = n*n*[0]
    xpts = [int((0.5+i)*img.size[0]/n) for i in range(n)]
    ypts = [int((0.5+i)*img.size[1]/n) for i in range(n)]
    for i,a in enumerate(xpts):
        for j,b in enumerate(ypts):
            points[i*n+j] = (a,b)
    return points

PIL_filters = ["BLUR", "SMOOTH_MORE", "CONTOUR", "DETAIL", "EDGE_ENHANCE",
        "EDGE_ENHANCE_MORE", "EMBOSS", "FIND_EDGES", "SHARPEN", "SMOOTH"]

def filter(filtername):
    #Takes the name of a filter and returns the filter in question
    #Seems like there should be an easier way to do this safely
    if filtername=="BLUR":
        return PIL.ImageFilter.BLUR
    elif filtername=="SMOOTH_MORE":
        return PIL.ImageFilter.SMOOTH_MORE
    elif filtername=="CONTOUR":
        return PIL.ImageFilter.CONTOUR
    elif filtername=="DETAIL":
        return PIL.ImageFilter.DETAIL
    elif filtername=="EDGE_ENHANCE":
        return PIL.ImageFilter.EDGE_ENHANCE
    elif filtername=="EDGE_ENHANCE_MORE":
        return PIL.ImageFilter.EDGE_ENHANCE_MORE
    elif filtername=="EMBOSS":
        return PIL.ImageFilter.EMBOSS
    elif filtername=="FIND_EDGES":
        return PIL.ImageFilter.FIND_EDGES
    elif filtername=="SHARPEN":
        return PIL.ImageFilter.SHARPEN
    elif filtername=="SMOOTH":
        return PIL.ImageFilter.SMOOTH

def shade_by_size(img):
    #img.show()
    #extracts the 256 largest elements and shades them based on size
    ma, n = scipy.ndimage.measurements.label(np.array(img.convert("L")))
    label, counts = np.unique(ma, return_counts=True)
    objects = [x for (y,x) in sorted(zip(counts[1:],label[1:]))][-255::]

    #This seems inefficient, there should be a better way of doing this
    tot = np.zeros(ma.shape)
    for i,o in enumerate(objects):
        tmp = ma.copy()
        tmp[tmp!=o]=0
        tmp[tmp==o]=i+1
        tot += tmp
    nimg = PIL.Image.new(mode="L", size=img.size)
    nimg.putdata(tot.flatten())
    return nimg

def gs_fft(img):
    fftvals = np.fft.fftshift(np.fft.fft2(img))
    freq = np.abs(fftvals)
    phase = np.angle(fftvals)
    freq = np.log(freq+np.e)
    fftmax = np.max(freq)
    freq = freq*255./fftmax
    return freq, phase, fftmax

def gs_ifft(img, phase, fftmax):
    fftvals = np.array(img)
    fftvals = ((np.exp(fftvals*fftmax/255))-np.e)
    ifftvals = np.real(np.fft.ifft2(np.fft.ifftshift(fftvals*np.exp(1j*phase))))
    return ifftvals

def linbin(img, bins = 12):
    ma, n = scipy.ndimage.measurements.label(np.array(img.convert("L")))
    label, counts = np.unique(ma, return_counts=True)
    label, counts = label[1:], counts[1:]

    sortedobjects = [x for (y,x) in sorted(zip(counts,label))]
    sortedcount = [counts[i-1] for i in sortedobjects]
    #print sortedcount[-1]
    linbins = np.linspace(0,sortedcount[-1], bins)
    
    #maybe improve this
    bin_divide = list(np.digitize(sortedcount, linbins))
    bin_count = [bin_divide.count(i) for i in range(bins)]
    
    return (bin_count, linbins)

def logbin(img, bins = 12):
    #counts size of each binary structure, puts into log bins
    ma, n = scipy.ndimage.measurements.label(np.array(img.convert("L")))
    label, counts = np.unique(ma, return_counts=True)
    label, counts = label[1:], counts[1:]

    sortedobjects = [x for (y,x) in sorted(zip(counts,label))]
    sortedcount = [counts[i-1] for i in sortedobjects]
    #print sortedcount[-1]
    logbins = np.exp(np.linspace(0,np.log(sortedcount[-1]), bins))
    
    #maybe improve this
    bin_divide = list(np.digitize(sortedcount, logbins))
    bin_count = [bin_divide.count(i) for i in range(bins)]
    
    return (bin_count, logbins)
