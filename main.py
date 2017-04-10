import PIL
import PIL.Image
import PIL.ImageTk
import PIL.ImageDraw
import PIL.ImageChops
import PIL.ImageFilter
import Tkinter
import tkFileDialog
import scipy
import scipy.ndimage
import numpy as np

#img is the primary display buffer
#should also split up into multiple files and separate UI from methods


logger=True

class imagewindow:
    def __init__(self, master):
        self.master = master
        menubar = Tkinter.Menu(master)
        
        self.img = PIL.Image.new('L', (800,500)) #keep here for now
        self.gimg = PIL.Image.new('L', (800,500)) #keep here for now
        self.bimg = PIL.Image.new('1', (800,500)) #keep here for now

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

        self.show_binary = True

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
        menubar.add_command(label="Toggle binary", command=self.toggle_binary)

        gsops = Tkinter.Menu(menubar, tearoff=0)
        gsops.add_command(label="FIND_EDGES", command=self.edgedetect_find_edges)
        gsops.add_command(label="Sobel", command=self.sobel)
        wshedmenu = Tkinter.Menu(gsops, tearoff=0)
        wshedmenu.add_command(label="Watershed", command=lambda : self.select_points(self.watershed))
        wshedmenu.add_command(label="Watershed grid", command=lambda : self.grid(self.watershed))
        gsops.add_cascade(label="Watershed", menu=wshedmenu)
        menubar.add_cascade(label="GS Operations", menu=gsops)

        binops = Tkinter.Menu(menubar, tearoff=0)
        binops.add_command(label="Logical operations", command = self.logical_operators)
        binops.add_command(label="Invert", command = self.invert_bin)
        binops.add_command(label="Erosion", command = lambda : self.morph_options("erosion"))
        binops.add_command(label="Dilation", command =lambda : self.morph_options("dilation"))
        binops.add_command(label="Opening", command = lambda : self.morph_options("opening"))
        binops.add_command(label="Closing", command = lambda : self.morph_options("closing"))
        binops.add_command(label="Fill holes", command = self.fill_option)

        menubar.add_cascade(label="Bin Operations", menu=binops)
        master.config(menu=menubar)
        master.protocol("WM_DELETE_WINDOW", quit) #close by x button does the same as quit
        
        image1 = PIL.ImageTk.PhotoImage(self.img)
        self.label_image = Tkinter.Label(self.master, image=image1)
        self.label_image.pack(side = "bottom", fill = "both", expand = "yes")
        master.mainloop()

    def update(self):
        #updates the current instance
        #there seems to be a recursion issue here
        colors = self.bimg.convert(mode="RGBA").split()
        colors[0].paste(self.bimg.point(lambda p: (255*(p>0))))
        colors[1].paste(self.bimg.point(lambda p: 0))
        colors[2].paste(self.bimg.point(lambda p: 0))
        
        bin = PIL.Image.merge("RGBA", colors) #maybe change to something more reasonable
        bin.putalpha(self.bimg)

        self.img = self.gimg.convert("RGBA")
        if self.show_binary:
            self.img.paste(bin, (0,0), bin)
        self.master.geometry('%dx%d' % (self.gimg.size[0],self.gimg.size[1]))
        
        image1 = PIL.ImageTk.PhotoImage(self.img)
        self.label_image.configure(image = image1)
        self.label_image.image = image1

    def greyscale(self):
        #convert to greyscale - deprecated? we do this by default with openfile()
        self.img = self.img.convert('L')
        self.update()

    def toggle_binary(self):
        if self.show_binary==1:
            self.show_binary=0
        else:
            self.show_binary=1
        self.update()

    def openfile(self):
        filename = tkFileDialog.askopenfilename()
        self.img = PIL.Image.open(filename).convert("L") #remove?
        self.gimg = PIL.Image.open(filename).convert("L")
        self.bimg = PIL.Image.new(mode="1", size=self.gimg.size)
        self.orig = self.img.copy()
        print 'img = PIL.Image.open("' + filename + '").convert("L")'
        print 'orig = img.copy()'
        self.update()

    def load_original(self):
        self.gimg = self.orig
        self.update()

    def savefile(self):
        filename = tkFileDialog.asksaveasfilename()
        self.img.save(filename)
        print 'img.save("' + filename + '")'

    def select_points(self,fn):
        #Permits the user to interactively select points which are then returned as a list
        img_copy = self.gimg.copy()
        self.gimg = self.gimg.convert(mode="RGB")
        point_list = []
        def point_self(all=None):
            r = 2
            if all!=None: # start may be unnecessary
                start=0
            else:
                start=-1
            draw = PIL.ImageDraw.Draw(self.gimg)
            for i in point_list[start:]:
                draw.ellipse((i[0]-r, i[1]-r, i[0]+r, i[1]+r), fill=(0,255,0,0))
            self.update()
            del draw

        def click(event):
            point_list.append((event.x, event.y))
            point_self()

        def undo(self):
            if point_list!=[]:
                point_list.pop()
            self.gimg = img_copy.convert(mode="RGB")
            point_self(1)
            self.update()

        def done(event=None):
            root.unbind("<ButtonPress-1>")
            root.unbind("<ButtonPress-3>")
            self.gimg = img_copy.copy()
            select_window.destroy()
            fn(point_list)

        root.bind("<ButtonPress-1>", click)
        root.bind("<ButtonPress-3>", done)
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
        #old = self.gimg.copy()
        old = self.gimg.copy()
        self.bimg = PIL.Image.new('1', self.gimg.size)
        hist_img = visual_histogram(self.gimg)
        
        self.mval = 0 #think about better solution here
        self.uval = 256

        if hist_img==None:
            #Only accept greyscale image
            return
        def cancel(event=None):
            #reverts to the old image and closes the histogram prompt
            histogram_window.destroy()
            self.gimg=old
            self.update()
        def click(event=None, mv=None):
            #shows red binary image superimposed on input greyscale image
            if mv==None:
                self.mval = event.x # cut-off value
                if self.mval > self.uval:
                    self.uval=self.mval
            elif mv==1:
                # this gives the default self.mval (cut-off value) as the maximum value position in the histogram
                # this isn't ideal but it's better than no default value at all
                self.mval = self.gimg.histogram().index(max(self.gimg.histogram()))
            colors = old.convert(mode="RGBA").split()
            colors[0].paste(old.point(lambda p: (255*(p>self.mval)*(p<self.uval))))
            colors[1].paste(old.point(lambda p: 0))
            colors[2].paste(old.point(lambda p: 0))
            
            self.bin = PIL.Image.merge("RGBA", colors) #maybe change to something more reasonable
            self.bin.putalpha(old.point(lambda p: (255*(p>self.mval)*(p<self.uval))))
            self.gimg = PIL.Image.blend(self.bin, old.convert(mode="RGBA"), 0.5)
            himg = visual_histogram(old, self.mval, self.uval)
            image1 = PIL.ImageTk.PhotoImage(himg)
            label_image.configure(image = image1)
            label_image.image = image1
            self.show_binary = True # This is preferred
            self.update()

        def rightclick(event):
            self.uval = event.x
            if self.uval<self.mval:
                self.uval=self.mval
            click(mv=2)

        def apply(event=None):
            #applies changes and saves the binary image in the primary buffer
            self.bimg = self.bin.convert(mode='L').point(lambda p: 255*(p>0)).convert(mode='1')
            self.gimg = old
            if logger:
                print ""
                #change later
                #print "img = img.point(lambda p: 255*(p>" + str(self.mval) + "))"
            histogram_window.destroy()
            self.update()

        histogram_window = Tkinter.Toplevel()
        histogram_window.bind("<ButtonPress-1>", click)
        histogram_window.bind("<ButtonPress-3>", rightclick)
        histogram_window.bind("<Return>", apply)
        histogram_window.bind("<Escape>", cancel)
        histogram_window.protocol("WM_DELETE_WINDOW", cancel) #close by x button does the same as cancel
        histogram_window.geometry('%dx%d' % (hist_img.size[0],hist_img.size[1]))

        tkpi = PIL.ImageTk.PhotoImage(hist_img)
        label_image = Tkinter.Label(histogram_window, image=tkpi)
        label_image.place(x=0,y=0,width=hist_img.size[0],height=hist_img.size[1])
        
        hmenubar = Tkinter.Menu(root)

        hmenubar.add_command(label="apply", command=apply)
        histogram_window.config(menu=hmenubar)
        click(mv=1)
        histogram_window.mainloop()

    def save_gs_buffer(self, slot):
        #save primary buffer (img) in greyscale buffer number (slot - int)
        #print slot
        if self.gimg.mode=="L":
            self.greyscale_buffers[slot] = self.gimg.copy()
            print "# Saved in greyscale slot " + str(slot)
            if logger:
                print "gimg" + str(slot)  + " = gimg.copy()"
        else:
            print "# Not a greyscale image"

    def load_gs_buffer(self, slot):
        #load greyscale buffer number (slot - int) into primary buffer (img)
        if self.greyscale_buffers[slot]!=None:
            self.gimg = self.greyscale_buffers[slot]
            print "# Loaded greyscale slot " + str(slot)
            if logger:
                print "gimg = gimg" + str(slot)
            self.update()
        else:
            print "# greyscale buffer " + str(slot) + " is empty"

    def save_bin_buffer(self, slot):
        #save primary buffer (img) in binary buffer number (slot - int)
        if self.bimg.mode=="1":
            self.bin_buffers[slot] = self.bimg.copy()
            print "# Saved in greyscale slot " + str(slot)
            if logger:
                print "binimg" + str(slot)  + " = img.copy()"
        else:
            print "# Not a binary image"

    def load_bin_buffer(self, slot):
        #load binary buffer number (slot - int) into primary buffer (img)
        if self.bin_buffers[slot]!=None:
            self.bimg = self.bin_buffers[slot]
            print "# Loaded binary slot " + str(slot)
            if logger:
                print "img = binimg" + str(slot)
            self.update()
        else:
            print "# binary buffer " + str(slot) + " is empty"

    def edgedetect_find_edges(self):
        self.gimg = self.gimg.filter(PIL.ImageFilter.FIND_EDGES)
        self.update()

    def sobel(self):
        tmpdata = scipy.ndimage.sobel(self.gimg)
        self.gimg.putdata(tmpdata.flatten())
        self.update()

    def grid(self, fn):
        # runs the function with n by n grid, n supplied by user
        top = Tkinter.Toplevel(height=300, width=300)
        e = Tkinter.Entry(top)
        e.insert(Tkinter.END, "5")
        e.pack()
        e.focus_set()
        def construct_run(n): #this can probably be simplified
            xpts = [int((0.5+i)*self.img.size[0]/n) for i in range(n)]
            ypts = [int((0.5+i)*self.img.size[1]/n) for i in range(n)]
            points = n*n*[0]
            for i,a in enumerate(xpts):
                for j,b in enumerate(ypts):
                    points[i*n+j] = (a,b)
            fn(points)
            top.destroy()
        
        apply = Tkinter.Button(top, text="apply", width=10, 
                command=lambda : construct_run(int(e.get())))
        apply.pack()
        
        
    def watershed(self, point_list):
        markers = np.zeros(self.gimg.size).astype(np.int16)
        factor = 256/len(point_list) # just for looks and easier use
        for i,a in enumerate(point_list):
            markers[a[0],a[1]] = i*factor
        imin = np.array(self.gimg).transpose()
        imgdata = scipy.ndimage.measurements.watershed_ift(np.array(imin), markers)
        img2 = PIL.Image.new(mode="L", size=self.gimg.size)
        img2.putdata(imgdata.transpose().flatten())
        if logger:
            print "points = " + str(point_list)
            print "# ATTENTION REQUIRED - definition on markers missing"
            print "data = scipy.ndimage.measurements.watershed_ift(np.array(" +str(imin) + "), markers)"
        self.gimg = img2.copy()
        self.update()

    def invert_bin(self):
        if self.bimg.mode=="1":
            self.bimg = PIL.ImageChops.invert(self.bimg)
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
                self.bimg = PIL.ImageChops.logical_and(self.img,self.bin_buffers[num])
                if logger:
                    print "img = (img and binimg" + str(num) + ")"
            top.destroy()
            self.update()
        def OR():
            num = int(e.get())
            if self.bin_buffers!=None:
                self.bimg = PIL.ImageChops.logical_or(self.img,self.bin_buffers[num])
                if logger:
                    print "img = (img or binimg" + str(num) + ")"
            top.destroy()
            self.update()
        def XOR():
            num = int(e.get())
            if self.bin_buffers!=None:
                tmp1 = PIL.ImageChops.subtract(self.bimg, self.bin_buffers[num])
                tmp2 = PIL.ImageChops.subtract(self.bin_buffers[num], self.bimg)
                self.bimg = PIL.ImageChops.logical_or(tmp1, tmp2)
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


    def morph(self, count, top, morph_type):
        #Should change to non-square shape for nicer binary images
        img2 = PIL.Image.new(mode="1", size = self.bimg.size)
        #convertion to "L" was necessary for some reason
        if morph_type=="erosion":
            tmp = scipy.ndimage.binary_erosion(self.bimg.convert("L"), 
                    structure=unit_disk(count))
        elif morph_type=="dilation":
            tmp = scipy.ndimage.binary_dilation(self.bimg.convert("L"), 
                    structure=unit_disk(count))
        elif morph_type=="opening":
            tmp = scipy.ndimage.binary_opening(self.bimg.convert("L"), 
                    structure=unit_disk(count))
        elif morph_type=="closing":
            tmp = scipy.ndimage.binary_closing(self.bimg.convert("L"), 
                    structure=unit_disk(count))
        elif morph_type=="fill_holes":
            tmp = scipy.ndimage.binary_fill_holes(self.bimg.convert("L"))

        if logger:
            print 'tmpimg = PIL.Image.new(mode="1", size = img.size)'
            print ('tmp = scipy.ndimage.binary_' +
            morph_type + '(img.convert("L"), structure=unit_disk(' +
                    str(count) + "))")
            print "tmpimg.putdata(255*tmp.flatten())"
            print 'img = tmpimg.convert("1")'

        img2.putdata(255*tmp.flatten())
        top.destroy()
        self.bimg = img2.convert("1")
        self.update()


    def morph_options(self, morph_type):
        top = Tkinter.Toplevel(height=300, width=300)
        e = Tkinter.Entry(top)
        e.insert(Tkinter.END, "2")
        e.pack()
        e.focus_set()
        apply = Tkinter.Button(top, text="apply", width=10, 
                command=lambda : self.morph(int(e.get()), top, morph_type))
        apply.pack()

    def fill_option(self):
        #maybe add options later though not strictly necessary
        top = Tkinter.Toplevel(height=300, width=300)
        self.morph(0,top,"fill_holes")

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


root = Tkinter.Tk()
root.geometry('+%d+%d' % (100,100))

gui = imagewindow(root)
