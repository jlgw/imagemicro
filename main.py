import PIL
import PIL.Image
import PIL.ImageTk
import PIL.ImageDraw
import PIL.ImageFilter
import Tkinter
import tkFileDialog
import scipy
import scipy.ndimage
import numpy as np

#img is the primary display buffer
#should get rid of the global variables later and change to lambda functions
#should also split up into multiple files and separate UI from methods

gs_buffer_count = 4
bin_buffer_count = 4
greyscale_buffers = gs_buffer_count*[None]
bin_buffers = bin_buffer_count*[None]

logger=True

def visual_histogram(img):
    #input: image (preferably greyscale)
    #output: visual histogram as image, binary image, 256x100
    hist = img.histogram()
    if len(hist)!=256:
        print "Convert to greyscale first"
        return
    maxval = max(hist)
    normheight = [int(100*i/maxval) for i in hist]
    histimg = PIL.Image.new(mode="1", size=(256,100))
    draw = PIL.ImageDraw.Draw(histimg)
    for i in range(256):
        draw.line((i,100) + (i,100-normheight[i]), fill=255)
    del draw
    return histimg

root = Tkinter.Tk()
root.geometry('+%d+%d' % (100,100))

def update(instance, image1):
    #updates the current instance
    instance.geometry('%dx%d' % (image1.size[0],image1.size[1]))
    curr = PIL.ImageTk.PhotoImage(image1)
    label_image = Tkinter.Label(instance, image=curr)
    label_image.place(x=0,y=0,width=image1.size[0],height=image1.size[1])
    instance.mainloop()
    instance.quit()

def greyscale():
    #convert to greyscale - deprecated? we do this by default with openfile()
    global img
    img = img.convert('L')
    update(root,img)

def openfile():
    global img
    global orig
    filename = tkFileDialog.askopenfilename()
    img = PIL.Image.open(filename).convert("L")
    orig = img.copy()
    print 'img = PIL.Image.open("' + filename + '").convert("L")'
    print 'orig = img.copy()'
    update(root,img)

def load_original():
    global img
    global orig
    img = orig
    update(root,img)

def savefile():
    global img
    filename = tkFileDialog.asksaveasfilename()
    img.save(filename)
    print 'img.save("' + filename + '")'

def threshold():
    #Opens a visual histogra
    #add thresholding to this
    global img
    global bin
    old = img.copy()
    hist_img = visual_histogram(img)
    if hist_img==None:
        #Only accept greyscale image
        return
    def cancel():
        #reverts to the old image and closes the histogram prompt
        global img
        histogram_window.destroy()
        img=old
        update(root,img)
    def click(event=None, mv=None):
        #shows red binary image superimposed on input greyscale image
        global img
        global bin
        global mval
        if mv==None:
            mval = event.x # cut-off value
        else:
            # this gives the default mval (cut-off value) as the maximum position in the histogram
            # this isn't ideal but it's better than no default value at all
            mval = img.histogram().index(max(img.histogram()))
        colors = old.convert(mode="RGBA").split()
        colors[0].paste(old.point(lambda p: (255*(p>mval))))
        colors[1].paste(old.point(lambda p: 0))
        colors[2].paste(old.point(lambda p: 0))
        
        bin = PIL.Image.merge("RGBA", colors)
        bin.putalpha(old.point(lambda p: (255*(p>mval))))
        img = PIL.Image.blend(bin, old.convert(mode="RGBA"), 0.5)
        update(root,img)

    def apply():
        #applies changes and saves the binary image in the primary buffer
        global bin
        global img
        img = bin.convert(mode='L').point(lambda p: 255*(p>0)).convert(mode='1')
        if logger:
            print "img = img.point(lambda p: 255*(p>" + str(mval) + "))"
        histogram_window.destroy()
        update(root,img)

    histogram_window = Tkinter.Toplevel()
    histogram_window.bind("<ButtonPress-1>", click)
    histogram_window.geometry('%dx%d' % (hist_img.size[0],hist_img.size[1]))

    tkpi = PIL.ImageTk.PhotoImage(hist_img)
    label_image = Tkinter.Label(histogram_window, image=tkpi)
    label_image.place(x=0,y=0,width=hist_img.size[0],height=hist_img.size[1])
    
    hmenubar = Tkinter.Menu(root)

    hmenubar.add_command(label="apply", command=apply)
    hmenubar.add_command(label="cancel", command=cancel)
    histogram_window.config(menu=hmenubar)
    click(mv=1)
    histogram_window.mainloop()

def save_gs_buffer(slot):
    #save primary buffer (img) in greyscale buffer number (slot - int)
    print slot
    global greyscale_buffers
    global img
    if img.mode=="L":
        greyscale_buffers[slot] = img.copy()
        print "# Saved in greyscale slot " + str(slot)
        if logger:
            print "img" + str(slot)  + " = img.copy()"
    else:
        print "# Not a greyscale image"

def load_gs_buffer(slot):
    #load greyscale buffer number (slot - int) into primary buffer (img)
    global greyscale_buffers
    global img
    if greyscale_buffers[slot]!=None:
        img = greyscale_buffers[slot]
        print "# Loaded greyscale slot " + str(slot)
        if logger:
            print "img = img" + str(slot)
        update(root, img)
    else:
        print "# greyscale buffer " + str(slot) + " is empty"

def save_bin_buffer(slot):
    #save primary buffer (img) in binary buffer number (slot - int)
    global bin_buffers
    global img
    if img.mode=="1":
        bin_buffers[slot] = img.copy()
        print "# Saved in greyscale slot " + str(slot)
        if logger:
            print "binimg" + str(slot)  + " = img.copy()"
    else:
        print "# Not a binary image"

def load_bin_buffer(slot):
    #load binary buffer number (slot - int) into primary buffer (img)
    global bin_buffers
    global img
    if bin_buffers[slot]!=None:
        img = bin_buffers[slot]
        print "# Loaded binary slot " + str(slot)
        if logger:
            print "img = binimg" + str(slot)
        update(root, img)
    else:
        print "# binary buffer " + str(slot) + " is empty"

def edgedetect_find_edges():
    global img
    img = img.filter(PIL.ImageFilter.FIND_EDGES)
    update(root,img)

def sobel():
    global img
    tmpdata = scipy.ndimage.sobel(img)
    img.putdata(tmpdata.flatten())
    update(root, img)

def invert_bin():
    global img
    if img.mode=="1":
        img = PIL.ImageChops.invert(img)
        update(root,img)
    else:
        print "Not a binary image"

def logical_operators():
    top = Tkinter.Toplevel(height=300, width=300)
    e = Tkinter.Entry(top)
    e.insert(Tkinter.END, "1")
    e.pack()
    e.focus_set()
    def AND():
        global img
        global bin_buffers
        num = int(e.get())
        if bin_buffers!=None:
            img = PIL.ImageChops.logical_and(img,bin_buffers[num])
            if logger:
                print "img = (img and binimg" + str(num) + ")"
        top.destroy()
        update(root,img)
    def OR():
        global img
        global bin_buffers
        num = int(e.get())
        if bin_buffers!=None:
            img = PIL.ImageChops.logical_or(img,bin_buffers[num])
            if logger:
                print "img = (img or binimg" + str(num) + ")"
        top.destroy()
        update(root,img)
    def XOR():
        global img
        global bin_buffers
        num = int(e.get())
        if bin_buffers!=None:
            tmp1 = PIL.ImageChops.subtract(img, bin_buffers[num])
            tmp2 = PIL.ImageChops.subtract(bin_buffers[num], img)
            img = PIL.ImageChops.logical_or(tmp1, tmp2)
            if logger:
                print "tmp1 = PIL.ImageChops.subtract(img, binimg" + str(num) + ")"
                print "tmp2 = PIL.ImageChops.subtract(binimg" + str(num) + ", img)"
        top.destroy()
        update(root,img)

    l_and= Tkinter.Button(top, text="AND", width=10, command=AND)
    l_or = Tkinter.Button(top, text="OR", width=10, command=OR)
    l_xor = Tkinter.Button(top, text="XOR", width=10, command=XOR)
    l_and.pack()
    l_or.pack()
    l_xor.pack()

    top.mainloop()

def morph(count, top, morph_type):
    #Should change to non-square shape for nicer binary images
    global img
    img2 = PIL.Image.new(mode="1", size = img.size)
    #convertion to "L" was necessary for some reason
    if morph_type=="erosion":
        tmp = scipy.ndimage.binary_erosion(img.convert("L"), 
                structure=np.ones((count, count)))
    elif morph_type=="dilation":
        tmp = scipy.ndimage.binary_dilation(img.convert("L"), 
                structure=np.ones((count, count)))
    elif morph_type=="opening":
        tmp = scipy.ndimage.binary_opening(img.convert("L"), 
                structure=np.ones((count, count)))
    elif morph_type=="closing":
        tmp = scipy.ndimage.binary_closing(img.convert("L"), 
                structure=np.ones((count, count)))
    elif morph_type=="fill_holes":
        tmp = scipy.ndimage.binary_fill_holes(img.convert("L"))
    if logger:
        print 'tmpimg = PIL.Image.new(mode="1", size = img.size)'
        print ('tmp = scipy.ndimage.binary_' +
        morph_type + '(img.convert("L"), structure=np.ones((' +
                str(count) + ", " + str(count) + "))")
        print "tmpimg.putdata(255*tmp.flatten())"
        print 'img = tmpimg.convert("1")'

    img2.putdata(255*tmp.flatten())
    top.destroy()
    img = img2.convert("1")
    update(root, img)


def erosion_option():
    top = Tkinter.Toplevel(height=300, width=300)
    e = Tkinter.Entry(top)
    e.insert(Tkinter.END, "2")
    e.pack()
    e.focus_set()
    apply = Tkinter.Button(top, text="apply", width=10, 
            command=lambda : morph(int(e.get()), top, "erosion"))
    apply.pack()
    top.mainloop()

def dilation_option():
    top = Tkinter.Toplevel(height=300, width=300)
    e = Tkinter.Entry(top)
    e.insert(Tkinter.END, "2")
    e.pack()
    e.focus_set()
    apply = Tkinter.Button(top, text="apply", width=10, 
            command=lambda : morph(int(e.get()), top, "dilation"))
    apply.pack()
    top.mainloop()

def opening_option():
    top = Tkinter.Toplevel(height=300, width=300)
    e = Tkinter.Entry(top)
    e.insert(Tkinter.END, "2")
    e.pack()
    e.focus_set()
    apply = Tkinter.Button(top, text="apply", width=10, 
            command=lambda : morph(int(e.get()), top, "opening"))
    apply.pack()
    top.mainloop()

def closing_option():
    top = Tkinter.Toplevel(height=300, width=300)
    e = Tkinter.Entry(top)
    e.insert(Tkinter.END, "2")
    e.pack()
    e.focus_set()
    apply = Tkinter.Button(top, text="apply", width=10, 
            command=lambda : morph(int(e.get()), top, "closing"))
    apply.pack()
    top.mainloop()

def fill_option():
    #maybe add options later though not strictly necessary
    top = Tkinter.Toplevel(height=300, width=300)
    morph(0,top,"fill_holes")

menubar = Tkinter.Menu(root)

filemenu = Tkinter.Menu(menubar, tearoff=0)
filemenu.add_command(label="Open", command=openfile)
filemenu.add_command(label="Save as", command=savefile)
filemenu.add_command(label="Quit", command=root.destroy)
menubar.add_cascade(label="File", menu=filemenu)

savemenu = Tkinter.Menu(menubar, tearoff=0)
#Would be nice to have this variable but became buggy when attempted
savemenu.add_command(label="Greyscale Slot 1", command=lambda : save_gs_buffer(1))
savemenu.add_command(label="Greyscale Slot 2", command=lambda : save_gs_buffer(2))
savemenu.add_command(label="Greyscale Slot 3", command=lambda : save_gs_buffer(3))
savemenu.add_command(label="Binary Slot 1", command=lambda : save_bin_buffer(1))
savemenu.add_command(label="Binary Slot 2", command=lambda : save_bin_buffer(2))
savemenu.add_command(label="Binary Slot 3", command=lambda : save_bin_buffer(3))
menubar.add_cascade(label="Save buffer", menu=savemenu)

loadmenu = Tkinter.Menu(menubar, tearoff=0)
loadmenu.add_command(label="Original", command=load_original)
loadmenu.add_command(label="Greyscale Slot 1", command=lambda : load_gs_buffer(1))
loadmenu.add_command(label="Greyscale Slot 2", command=lambda : load_gs_buffer(2))
loadmenu.add_command(label="Greyscale Slot 3", command=lambda : load_gs_buffer(3))
loadmenu.add_command(label="Binary Slot 1", command=lambda : load_bin_buffer(1))
loadmenu.add_command(label="Binary Slot 2", command=lambda : load_bin_buffer(2))
loadmenu.add_command(label="Binary Slot 3", command=lambda : load_bin_buffer(3))
menubar.add_cascade(label="Load buffer", menu=loadmenu)

menubar.add_command(label="Threshold", command=threshold)

gsops = Tkinter.Menu(menubar, tearoff=0)
gsops.add_command(label="FIND_EDGES", command=edgedetect_find_edges)
gsops.add_command(label="Sobel", command=sobel)
menubar.add_cascade(label="GS Operations", menu=gsops)

binops = Tkinter.Menu(menubar, tearoff=0)
binops.add_command(label="Logical operations", command = logical_operators)
binops.add_command(label="Invert", command = invert_bin)
binops.add_command(label="Erosion", command = erosion_option)
binops.add_command(label="Dilation", command = dilation_option)
binops.add_command(label="Opening", command = opening_option)
binops.add_command(label="Closing", command = closing_option)
binops.add_command(label="Fill holes", command = fill_option)

menubar.add_cascade(label="Bin Operations", menu=binops)

root.config(menu=menubar)

img = PIL.Image.new('L', (800,500))
update(root, img)

