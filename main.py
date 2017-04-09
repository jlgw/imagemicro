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


logger=True

class imagewindow:
    def __init__(self, master):
        self.master = master
        menubar = Tkinter.Menu(master)
        
        self.img = PIL.Image.new('L', (800,500)) #keep here for now

        filemenu = Tkinter.Menu(menubar, tearoff=0)
        filemenu.add_command(label="Open", command=self.openfile)
        filemenu.add_command(label="Save as", command=self.savefile)
        filemenu.add_command(label="Quit", command=master.destroy)
        menubar.add_cascade(label="File", menu=filemenu)
        
        #constants, move to better place
        self.gs_buffer_count = 4
        self.bin_buffer_count = 4
        self.greyscale_buffers = self.gs_buffer_count*[None]
        self.bin_buffers = self.bin_buffer_count*[None]

        savemenu = Tkinter.Menu(menubar, tearoff=0)
        #Would be nice to have this variable but became buggy when attempted
        savemenu.add_command(label="Greyscale Slot 1", command=lambda : self.save_gs_buffer(1))
        savemenu.add_command(label="Greyscale Slot 2", command=lambda : self.save_gs_buffer(2))
        savemenu.add_command(label="Greyscale Slot 3", command=lambda : self.save_gs_buffer(3))
        savemenu.add_command(label="Binary Slot 1", command=lambda : self.save_bin_buffer(1))
        savemenu.add_command(label="Binary Slot 2", command=lambda : self.save_bin_buffer(2))
        savemenu.add_command(label="Binary Slot 3", command=lambda : self.save_bin_buffer(3))
        menubar.add_cascade(label="Save buffer", menu=savemenu)

        loadmenu = Tkinter.Menu(menubar, tearoff=0)
        loadmenu.add_command(label="Original", command=self.load_original)
        loadmenu.add_command(label="Greyscale Slot 1", command=lambda : self.load_gs_buffer(1))
        loadmenu.add_command(label="Greyscale Slot 2", command=lambda : self.load_gs_buffer(2))
        loadmenu.add_command(label="Greyscale Slot 3", command=lambda : self.load_gs_buffer(3))
        loadmenu.add_command(label="Binary Slot 1", command=lambda : self.load_bin_buffer(1))
        loadmenu.add_command(label="Binary Slot 2", command=lambda : self.load_bin_buffer(2))
        loadmenu.add_command(label="Binary Slot 3", command=lambda : self.load_bin_buffer(3))
        menubar.add_cascade(label="Load buffer", menu=loadmenu)

        menubar.add_command(label="Threshold", command=self.threshold)

        gsops = Tkinter.Menu(menubar, tearoff=0)
        gsops.add_command(label="FIND_EDGES", command=self.edgedetect_find_edges)
        gsops.add_command(label="Sobel", command=self.sobel)
        gsops.add_command(label="Watersheds", command=lambda : self.select_points(self.watersheds))
        menubar.add_cascade(label="GS Operations", menu=gsops)

        binops = Tkinter.Menu(menubar, tearoff=0)
        binops.add_command(label="Logical operations", command = self.logical_operators)
        binops.add_command(label="Invert", command = self.invert_bin)
        binops.add_command(label="Erosion", command = self.erosion_option)
        binops.add_command(label="Dilation", command = self.dilation_option)
        binops.add_command(label="Opening", command = self.opening_option)
        binops.add_command(label="Closing", command = self.closing_option)
        binops.add_command(label="Fill holes", command = self.fill_option)

        menubar.add_cascade(label="Bin Operations", menu=binops)
        master.config(menu=menubar)
        master.protocol("WM_DELETE_WINDOW", quit) #close by x button does the same as quit

    def update(self):
        image1 = PIL.ImageTk.PhotoImage(self.img)
        #updates the current instance
        #there seems to be a recursion issue here
        self.master.geometry('%dx%d' % (self.img.size[0],self.img.size[1]))
        label_image = Tkinter.Label(self.master, image=image1)
        label_image.place(x=0,y=0,width=self.img.size[0],height=self.img.size[1])
        self.master.quit()
        self.master.mainloop()

    def greyscale(self):
        #convert to greyscale - deprecated? we do this by default with openfile()
        self.img = self.img.convert('L')
        self.update()

    def openfile(self):
        filename = tkFileDialog.askopenfilename()
        self.img = PIL.Image.open(filename).convert("L")
        self.orig = self.img.copy()
        print 'img = PIL.Image.open("' + filename + '").convert("L")'
        print 'orig = img.copy()'
        self.update()

    def load_original(self):
        self.img = self.orig
        self.update()

    def savefile(self):
        filename = tkFileDialog.asksaveasfilename()
        self.img.save(filename)
        print 'img.save("' + filename + '")'

    def select_points(self,fn):
        #ugly hack to solve the return issue, try to come up with better solution
        #Permits the user to interactively select points which are then returned as a list
        img_copy = self.img.copy()
        self.img = self.img.convert(mode="RGB")
        point_list = []
        def point_self(all=None):
            r = 2
            if all!=None: # start may be unnecessary
                start=0
            else:
                start=-1
            draw = PIL.ImageDraw.Draw(self.img)
            for i in point_list[start:]:
                draw.ellipse((i[0]-r, i[1]-r, i[0]+r, i[1]+r), fill=(0,255,0,0))
            self.update()

        def click(event):
            point_list.append((event.x, event.y))
            point_self()

        def undo(self):
            if point_list!=[]:
                point_list.pop()
            self.img = img_copy.convert(mode="RGB")
            point_self(1)
            self.update()

        def done():
            root.unbind("<ButtonPress-1>")
            self.img = img_copy.copy()
            select_window.destroy()
            fn(point_list)

        root.bind("<ButtonPress-1>", click)
        select_window = Tkinter.Toplevel()
        
        hmenubar = Tkinter.Menu(root)

        hmenubar.add_command(label="done", command=done)
        hmenubar.add_command(label="undo", command=undo)
        select_window.protocol("WM_DELETE_WINDOW", done)
        select_window.config(menu=hmenubar)
            
    def threshold(self):
        # Opens a visual histogram
        # add thresholding to this
        # would be nice to have better interactivity
        # not sure what option would be best for that
        global bin
        old = self.img.copy()
        hist_img = visual_histogram(self.img)
        if hist_img==None:
            #Only accept greyscale image
            return
        def cancel():
            #reverts to the old image and closes the histogram prompt
            histogram_window.destroy()
            self.img=old
            self.update()
        def click(event=None, mv=None):
            #shows red binary image superimposed on input greyscale image
            global bin
            global mval
            if mv==None:
                mval = event.x # cut-off value
            else:
                # this gives the default mval (cut-off value) as the maximum position in the histogram
                # this isn't ideal but it's better than no default value at all
                mval = self.img.histogram().index(max(self.img.histogram()))
            colors = old.convert(mode="RGBA").split()
            colors[0].paste(old.point(lambda p: (255*(p>mval))))
            colors[1].paste(old.point(lambda p: 0))
            colors[2].paste(old.point(lambda p: 0))
            
            bin = PIL.Image.merge("RGBA", colors)
            bin.putalpha(old.point(lambda p: (255*(p>mval))))
            self.img = PIL.Image.blend(bin, old.convert(mode="RGBA"), 0.5)
            self.update()

        def apply():
            #applies changes and saves the binary image in the primary buffer
            global bin
            self.img = bin.convert(mode='L').point(lambda p: 255*(p>0)).convert(mode='1')
            if logger:
                print "img = img.point(lambda p: 255*(p>" + str(mval) + "))"
            histogram_window.destroy()
            self.update()

        histogram_window = Tkinter.Toplevel()
        histogram_window.bind("<ButtonPress-1>", click)
        histogram_window.geometry('%dx%d' % (hist_img.size[0],hist_img.size[1]))

        tkpi = PIL.ImageTk.PhotoImage(hist_img)
        label_image = Tkinter.Label(histogram_window, image=tkpi)
        label_image.place(x=0,y=0,width=hist_img.size[0],height=hist_img.size[1])
        
        hmenubar = Tkinter.Menu(root)

        hmenubar.add_command(label="apply", command=apply)
        histogram_window.protocol("WM_DELETE_WINDOW", cancel) #close by x button does the same as cancel
        histogram_window.config(menu=hmenubar)
        click(mv=1)

    def save_gs_buffer(self, slot):
        #save primary buffer (img) in greyscale buffer number (slot - int)
        #print slot
        if self.img.mode=="L":
            self.greyscale_buffers[slot] = self.img.copy()
            print "# Saved in greyscale slot " + str(slot)
            if logger:
                print "img" + str(slot)  + " = img.copy()"
        else:
            print "# Not a greyscale image"

    def load_gs_buffer(self, slot):
        #load greyscale buffer number (slot - int) into primary buffer (img)
        if self.greyscale_buffers[slot]!=None:
            self.img = self.greyscale_buffers[slot]
            print "# Loaded greyscale slot " + str(slot)
            if logger:
                print "img = img" + str(slot)
            self.update()
        else:
            print "# greyscale buffer " + str(slot) + " is empty"

    def save_bin_buffer(self, slot):
        #save primary buffer (img) in binary buffer number (slot - int)
        if self.img.mode=="1":
            self.bin_buffers[slot] = self.img.copy()
            print "# Saved in greyscale slot " + str(slot)
            if logger:
                print "binimg" + str(slot)  + " = img.copy()"
        else:
            print "# Not a binary image"

    def load_bin_buffer(self, slot):
        #load binary buffer number (slot - int) into primary buffer (img)
        if self.bin_buffers[slot]!=None:
            self.img = self.bin_buffers[slot]
            print "# Loaded binary slot " + str(slot)
            if logger:
                print "img = binimg" + str(slot)
            self.update()
        else:
            print "# binary buffer " + str(slot) + " is empty"

    def edgedetect_find_edges(self):
        self.img = self.img.filter(PIL.ImageFilter.FIND_EDGES)
        self.update()

    def sobel(self):
        tmpdata = scipy.ndimage.sobel(self.img)
        self.img.putdata(tmpdata.flatten())
        self.update()

    def watersheds(self, point_list):
        markers = np.zeros(self.img.size).astype(np.int16)
        factor = 256/len(point_list) # just for looks and easier use
        for i,a in enumerate(point_list):
            markers[a[0],a[1]] = i*factor
        imin = np.array(self.img).transpose()
        imgdata = scipy.ndimage.measurements.watershed_ift(np.array(imin), markers)
        img2 = PIL.Image.new(mode="L", size=self.img.size)
        img2.putdata(imgdata.transpose().flatten())
        if logger:
            print "points = " + str(point_list)
            print "# ATTENTION REQUIRED - definition on markers missing"
            print "data = scipy.ndimage.measurements.watershed_ift(np.array(" +str(imin) + "), markers)"
        self.img = img2.copy()
        self.update()

    def invert_bin(self):
        if self.img.mode=="1":
            self.img = PIL.ImageChops.invert(self.img)
            self.update()
        else:
            print "Not a binary image"

    def logical_operators(self):
        top = Tkinter.Toplevel(height=300, width=300)
        e = Tkinter.Entry(top)
        e.insert(Tkinter.END, "1")
        e.pack()
        e.focus_set()
        def AND():
            num = int(e.get())
            if self.bin_buffers!=None:
                self.img = PIL.ImageChops.logical_and(self.img,bin_buffers[num])
                if logger:
                    print "img = (img and binimg" + str(num) + ")"
            top.destroy()
            self.update()
        def OR():
            num = int(e.get())
            if self.bin_buffers!=None:
                self.img = PIL.ImageChops.logical_or(self.img,bin_buffers[num])
                if logger:
                    print "img = (img or binimg" + str(num) + ")"
            top.destroy()
            self.update()
        def XOR():
            num = int(e.get())
            if self.bin_buffers!=None:
                tmp1 = PIL.ImageChops.subtract(self.img, self.bin_buffers[num])
                tmp2 = PIL.ImageChops.subtract(bin_buffers[num], self.img)
                self.img = PIL.ImageChops.logical_or(tmp1, tmp2)
                if logger:
                    print "tmp1 = PIL.ImageChops.subtract(img, binimg" + str(num) + ")"
                    print "tmp2 = PIL.ImageChops.subtract(binimg" + str(num) + ", img)"
            top.destroy()
            self.update()

        l_and= Tkinter.Button(top, text="AND", width=10, command=AND)
        l_or = Tkinter.Button(top, text="OR", width=10, command=OR)
        l_xor = Tkinter.Button(top, text="XOR", width=10, command=XOR)
        l_and.pack()
        l_or.pack()
        l_xor.pack()


    def morph(count, top, morph_type):
        #Should change to non-square shape for nicer binary images
        img2 = PIL.Image.new(mode="1", size = self.img.size)
        #convertion to "L" was necessary for some reason
        if morph_type=="erosion":
            tmp = scipy.ndimage.binary_erosion(self.img.convert("L"), 
                    structure=unit_disk(count))
        elif morph_type=="dilation":
            tmp = scipy.ndimage.binary_dilation(self.img.convert("L"), 
                    structure=unit_disk(count))
        elif morph_type=="opening":
            tmp = scipy.ndimage.binary_opening(self.img.convert("L"), 
                    structure=unit_disk(count))
        elif morph_type=="closing":
            tmp = scipy.ndimage.binary_closing(self.img.convert("L"), 
                    structure=unit_disk(count))
        elif morph_type=="fill_holes":
            tmp = scipy.ndimage.binary_fill_holes(self.img.convert("L"))

        if logger:
            print 'tmpimg = PIL.Image.new(mode="1", size = img.size)'
            print ('tmp = scipy.ndimage.binary_' +
            morph_type + '(img.convert("L"), structure=unit_disk(' +
                    str(count) + "))")
            print "tmpimg.putdata(255*tmp.flatten())"
            print 'img = tmpimg.convert("1")'

        img2.putdata(255*tmp.flatten())
        top.destroy()
        self.img = img2.convert("1")
        self.update()


    def erosion_option(self):
        top = Tkinter.Toplevel(height=300, width=300)
        e = Tkinter.Entry(top)
        e.insert(Tkinter.END, "2")
        e.pack()
        e.focus_set()
        apply = Tkinter.Button(top, text="apply", width=10, 
                command=lambda : morph(int(e.get()), top, "erosion"))
        apply.pack()

    def dilation_option(self):
        top = Tkinter.Toplevel(height=300, width=300)
        e = Tkinter.Entry(top)
        e.insert(Tkinter.END, "2")
        e.pack()
        e.focus_set()
        apply = Tkinter.Button(top, text="apply", width=10, 
                command=lambda : morph(int(e.get()), top, "dilation"))
        apply.pack()

    def opening_option(self):
        top = Tkinter.Toplevel(height=300, width=300)
        e = Tkinter.Entry(top)
        e.insert(Tkinter.END, "2")
        e.pack()
        e.focus_set()
        apply = Tkinter.Button(top, text="apply", width=10, 
                command=lambda : morph(int(e.get()), top, "opening"))
        apply.pack()

    def closing_option(self):
        top = Tkinter.Toplevel(height=300, width=300)
        e = Tkinter.Entry(top)
        e.insert(Tkinter.END, "2")
        e.pack()
        e.focus_set()
        apply = Tkinter.Button(top, text="apply", width=10, 
                command=lambda : morph(int(e.get()), top, "closing"))
        apply.pack()

    def fill_option(self):
        #maybe add options later though not strictly necessary
        top = Tkinter.Toplevel(height=300, width=300)
        morph(0,top,"fill_holes")

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

def unit_disk(n):
    #Make this faster later if necessary, probably good enough for now
    h = (n-1)/2.
    mat = np.ones((n,n))
    for i in range(n):
        for j in range(n):
            if (i-h)**2+(j-h)**2>h**2+.5:
                mat[i,j]=0
    return mat


root = Tkinter.Tk()
root.geometry('+%d+%d' % (100,100))

gui = imagewindow(root)
gui.update()
