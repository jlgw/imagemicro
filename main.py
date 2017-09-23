import os.path
import sys
import PIL
import PIL.Image
import PIL.ImageMath
import PIL.ImageTk
import PIL.ImageDraw
import PIL.ImageChops
import PIL.ImageFilter
import Tkinter
import tkFileDialog
import scipy
import scipy.fftpack
import scipy.ndimage
import numpy as np
import inspect
import matplotlib.pyplot as plt

from dialog_window import dialog_window
from tools import visual_histogram, unit_disk, watershed_pts, grid, \
        PIL_filters, filter, shade_by_size, linbin, logbin, \
        highest_2nd_derivative

#should also split up into multiple files and separate UI from methods

logger=True

class Imagewindow:
    def __init__(self, master, filename=None):
        self.master = master
        
        self.img  = PIL.Image.new('L', (800,500)) #keep here for now
        self.gimg = PIL.Image.new('L', (800,500))
        self.bimg = PIL.Image.new('1', (800,500))
        
        #constants, move to better place
        self.gs_buffer_count  = 3
        self.bin_buffer_count = 4
        self.greyscale_buffers = self.gs_buffer_count*[None]
        self.bin_buffers = self.bin_buffer_count*[None]

        self.show_binary  = True
        self.resize_ratio = 1

        self.max_undos = 8
        self.undolist  = []

        self.conversion_ratio = 0

        menubar = Tkinter.Menu(master)

        filemenu = Tkinter.Menu(menubar, tearoff=0)
        filemenu.add_command(label="Open", command=self.openfile)
        filemenu.add_command(label="Save as", command=self.savefile)
        filemenu.add_command(label="Quit", command=master.destroy)
        menubar.add_cascade(label="File", menu=filemenu)

        editmenu = Tkinter.Menu(menubar, tearoff=0)
        editmenu.add_command(label="Undo", command=self.undo)
        menubar.add_cascade(label="Edit", menu=editmenu)
        
        savemenu = Tkinter.Menu(menubar, tearoff=0)
        loadmenu = Tkinter.Menu(menubar, tearoff=0)
        loadmenu.add_command(label="Original greyscale", 
                command=self.load_original)

        for i in range(self.gs_buffer_count):
            savemenu.add_command(label="Greyscale Slot" + str(i), 
                    command=lambda i=i: self.save_gs_buffer(i))
            loadmenu.add_command(label="Greyscale Slot" + str(i), 
                    command=lambda i=i: self.load_gs_buffer(i))
        
        for i in range(self.bin_buffer_count):
            savemenu.add_command(label="Binary Slot" + str(i), 
                    command=lambda i=i: self.save_bin_buffer(i))
            loadmenu.add_command(label="Binary Slot" + str(i), 
                    command=lambda i=i: self.load_bin_buffer(i))

        menubar.add_cascade(label="Save buffer", menu=savemenu)
        menubar.add_cascade(label="Load buffer", menu=loadmenu)

        menubar.add_command(label="Threshold", 
                command=self.threshold)
        menubar.add_command(label="Toggle binary", 
                command=self.toggle_binary)

        gsops = Tkinter.Menu(menubar, tearoff=0)
        gsops.add_command(label="Apply binary mask",
                command=self.apply_binary_mask)
        gsops.add_command(label="Arithmetic",
                command=self.greyscale_arithmetic)
        gsops.add_command(label="FIND_EDGES", 
                command=self.edgedetect_find_edges)
        gsops.add_command(label="Sobel",
                command=self.sobel)
        
        filtermenu = Tkinter.Menu(gsops, tearoff=0)
        for f in PIL_filters:
            filtermenu.add_command(label=f, 
                    command=lambda f=f: self.filter_function(f))
        gsops.add_cascade(label="PIL Filters", menu=filtermenu)

        wshedmenu = Tkinter.Menu(gsops, tearoff=0)
        wshedmenu.add_command(label="Manual", 
                command=lambda : self.select_points(self.watershed))
        wshedmenu.add_command(label="Grid", 
                command=lambda : self.run_grid(self.watershed))
        gsops.add_cascade(label="Watershed", menu=wshedmenu)
        gsops.add_command(label="LUT Transform", command=self.lut)
        fftmenu = Tkinter.Menu(gsops, tearoff=0)
        fftmenu.add_command(label="Transform", command=self.fft)
        fftmenu.add_command(label="Inverse transform", command=self.ifft)
        gsops.add_cascade(label="FFT", menu=fftmenu)
        menubar.add_cascade(label="GS Operations", menu=gsops)

        binops = Tkinter.Menu(menubar, tearoff=0)
        binops.add_command(label="Logical operations", command = self.logical_operators)
        binops.add_command(label="Invert", 
                command = self.invert_bin)
        pickmenu = Tkinter.Menu(binops, tearoff=0)
        pickmenu.add_command(label="Reject",
                command = lambda : self.picker("reject"))
        pickmenu.add_command(label="Keep",
                command = lambda : self.picker("keep"))
        binops.add_cascade(label="Pick elements", menu=pickmenu)
        morphmenu = Tkinter.Menu(binops, tearoff=0)
        morphmenu.add_command(label="Erosion", 
                command = lambda : self.morph_options("erosion"))
        morphmenu.add_command(label="Dilation", 
                command = lambda : self.morph_options("dilation"))
        morphmenu.add_command(label="Opening", 
                command = lambda : self.morph_options("opening"))
        morphmenu.add_command(label="Closing", 
                command = lambda : self.morph_options("closing"))
        morphmenu.add_command(label="Fill holes", 
                command = self.fill_option)
        binops.add_cascade(label = "Morphology", menu=morphmenu)
        binops.add_command(label="Largest brightest", 
                command = self.shade)
        drawmenu = Tkinter.Menu(binops, tearoff=0)
        drawmenu.add_command(label="Draw center circle", command=self.draw_center_circle)
        binops.add_cascade(label = "Draw", menu=drawmenu)
        menubar.add_cascade(label="Bin Operations", menu=binops)

        anmenu = Tkinter.Menu(menubar, tearoff=0)
        anmenu.add_command(label="Log bins", 
                command = lambda : self.log_divide_bins())
        anmenu.add_command(label="Lin bins", 
                command = lambda : self.lin_divide_bins())
        anmenu.add_command(label="Set conversion ratio", 
                command = lambda : self.set_conversion_ratio())
        menubar.add_cascade(label="Analysis", menu=anmenu)

        master.config(menu=menubar)

        master.protocol("WM_DELETE_WINDOW", quit) #close by x button does the same as quit
        master.bind("<Control-z>", self.undo)
        master.bind("<plus>", self.zoom_in)
        master.bind("<minus>", self.zoom_out)
        
        image1 = PIL.ImageTk.PhotoImage(self.img)
        self.label_image = Tkinter.Label(self.master, image=image1)
        self.label_image.pack(side = "bottom", fill = "both", expand = "yes")
        
        if logger:
            print("from tools import *")
        if filename!=None:
            self.openfile(filename)
        master.mainloop()

    def update(self, saveundo=True):

        if saveundo:
            self.undolist.append((self.gimg.copy(),self.bimg.copy()))
            if len(self.undolist)>self.max_undos:
                self.undolist.pop(0)

        colors = self.bimg.convert(mode="RGBA").split()
        colors[0].paste(self.bimg.point(lambda p: (255*(p>0))))
        colors[1].paste(self.bimg.point(lambda p: 0))
        colors[2].paste(self.bimg.point(lambda p: 0))
        
        bin = PIL.Image.merge("RGBA", colors) 
        bin.putalpha(self.bimg)

        self.img = self.gimg.convert("RGBA")
        if self.show_binary:
            self.img.paste(bin, (0,0), bin)
        
        if self.resize_ratio!=1:
            self.img = self.img.resize((int(self.img.width*self.resize_ratio), 
                int(self.img.height*self.resize_ratio)))
        self.master.geometry('%dx%d' % (self.img.size[0],self.img.size[1]))
        
        image1 = PIL.ImageTk.PhotoImage(self.img)
        self.label_image.configure(image = image1)
        self.label_image.image = image1
        if(self.conversion_ratio!=0):
            cr = " ("+str(self.conversion_ratio) + " um per pixel)"
        else:
            cr = ""
        self.master.wm_title(self.filename + " (" + 
                str(self.resize_ratio*100) + "%)" + cr)
    
    def zoom_in(self, event=None):
        self.resize_ratio += 0.1
        self.update(False)

    def zoom_out(self, event=None):
        self.resize_ratio -= 0.1
        self.update(False)

    def undo(self, event=None):
        if len(self.undolist)>1:
            self.undolist.pop()
            self.gimg, self.bimg = self.undolist[-1]
            self.update(False)
            if logger:
                self.logging("# Undo")
        else:
            if logger:
                self.logging("# Max undos reached")

    def toggle_binary(self):
        if self.show_binary==1:
            self.show_binary=0
        else:
            self.show_binary=1
        self.update()

    def openfile(self, filename=None):
        if filename==None:
            filename = tkFileDialog.askopenfilename()
        if type(filename)==tuple or filename=="":
            return
        self.gimg = PIL.Image.open(filename).convert("L")
        self.bimg = PIL.Image.new(mode="1", size=self.gimg.size)
        self.orig = self.gimg.copy()
        fhead, ftail = os.path.split(filename)
        self.filename = ftail
        self.update()
        if logger:
            self.logging('gimg = PIL.Image.open("' + filename + '").convert("L")')
            self.logging('bimg = PIL.Image.new(mode="1", size=gimg.size)')
            self.logging('orig = gimg.copy()')

    def load_original(self):
        self.gimg = self.orig.copy()
        if logger:
            self.logging("gimg = orig.copy()")
        self.update()

    def savefile(self):
        filename = tkFileDialog.asksaveasfilename()
        if type(filename)==tuple or filename=="":
            return
        self.gimg.save(filename)
        self.logging('gimg.save("' + filename + '")')

    def select_points(self,fn):
        #Permits the user to interactively select points
        #function fn is then run with list of points as argument
        #right-click - done
        img_copy = self.gimg.copy()
        self.gimg = self.gimg.convert(mode="RGB")
        #print id(img_copy)
        #print id(self.gimg)
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
            self.update(False)
            del draw

        def click(event):
            point_list.append((int(event.x/self.resize_ratio), 
                int(event.y/self.resize_ratio)))
            point_self()

        def undo_pt(event=None):
            if point_list!=[]:
                point_list.pop()
            self.gimg = img_copy.convert(mode="RGB")
            point_self(1)
            self.update()

        def done(event=None, cancel=False):
            self.master.unbind("<ButtonPress-1>")
            self.master.unbind("<ButtonPress-3>")
            self.gimg = img_copy.copy()
            select_window.destroy()
            if cancel==False:
                fn(point_list)

        cancel = lambda event=None: done(event, cancel=True)
        self.master.bind("<ButtonPress-1>", click)
        self.master.bind("<ButtonPress-3>", done)

        select_window = dialog_window(self.master)
        select_window.setapply(done)
        select_window.setcancel(cancel)
        undobtn = Tkinter.Button(select_window, text="Undo", width=10, 
                command=undo_pt)
        undobtn.grid(row=2, column=0)
            
    def threshold(self):
        # Opens a visual histogram
        # add thresholding to this
        # would be nice to have better interactivity
        # not sure what option would be best for that
        gold = self.gimg.copy()
        bold = self.bimg.copy()
        self.bimg = PIL.Image.new('1', self.gimg.size)
        hist_img  = visual_histogram(self.gimg)
        
        self.mval = 0 
        self.uval = 256

        if hist_img==None:
            #Only accept greyscale image
            return
        def cancel(event=None):
            #reverts to the gold image and closes the histogram prompt
            histogram_window.destroy()
            self.gimg = gold
            self.bimg = bold
            self.update()

        def click(event=None, mv=None):
            if event!=None:
                if event.widget!=label_image:
                    return
            #shows red binary image superimposed on input greyscale image
            if mv==None:
                self.mval = event.x # cut-off value
                if self.mval > self.uval:
                    self.uval=self.mval
            elif mv==1:
                # this gives the default self.mval (cut-off value) as the maximum value position in the histogram
                # this isn't ideal but it's better than no default value at all
                #self.mval = self.gimg.histogram().index(max(self.gimg.histogram()))
                # alternate solution: highest 2nd derivative
                # if it works create self.logging for this
                # split into two lines where this function is called if mval is
                # unchanged
                self.mval = highest_2nd_derivative(self.gimg)
                
            colors = gold.convert(mode="RGBA").split()
            colors[0].paste(gold.point(
                lambda p: (255*(p>self.mval)*(p<self.uval))))
            colors[1].paste(gold.point(lambda p: 0))
            colors[2].paste(gold.point(lambda p: 0))
            
            self.bin = PIL.Image.merge("RGBA", colors) #change name? bin is bad
            self.bin.putalpha(gold.point(
                lambda p: (255*(p>self.mval)*(p<self.uval))))
            self.gimg = PIL.Image.blend(self.bin, gold.convert(mode="RGBA"), 0.5)
            himg = visual_histogram(gold, self.mval, self.uval)
            image1 = PIL.ImageTk.PhotoImage(himg)
            label_image.configure(image = image1)
            label_image.image = image1
            # People generally want the hidden binary to show when they use a
            # feature that operates on the binary image, even if it were hidden
            # before, so self.show_binary is set to True on all calls to such a
            # function
            self.show_binary = True 
            self.update(False)

        def rightclick(event):
            if event.widget!=label_image:
                return
            self.uval = event.x
            if self.uval<self.mval:
                self.uval=self.mval
            click(mv=2)

        def apply(event=None):
            #applies changes and saves the binary image in the primary buffer
            self.bimg = self.bin.convert(mode='L').point(lambda p: 255*(p>0)).convert(mode='1')
            self.gimg = gold
            if logger:
                self.logging("bimg = gimg.point(lambda p: 255*(p>" + str(self.mval) + ")).convert(mode='1')")
            histogram_window.destroy()
            self.update()

        histogram_window = dialog_window(self.master)
        histogram_window.wm_title("Threshold")
        histogram_window.bind("<ButtonPress-1>", click)
        histogram_window.bind("<ButtonPress-3>", rightclick)
        histogram_window.setapply(apply)
        histogram_window.setcancel(cancel)
        #histogram_window.geometry('%dx%d' % (hist_img.size[0],hist_img.size[1]))
        histogram_window.geometry('%dx%d' % (hist_img.size[0],hist_img.size[1]+30))

        tkpi = PIL.ImageTk.PhotoImage(hist_img)
        label_image = Tkinter.Label(histogram_window, image=tkpi)
        label_image.place(x=0,y=30,width=hist_img.size[0],height=hist_img.size[1])
        
        #hmenubar = Tkinter.Menu(self.master)

        #optionframe= Tkinter.Frame(histogram_window)
        #optionframe.pack(side=Tkinter.BOTTOM)
        #apply = Tkinter.Button(optionframe, text="apply", width=10, command=apply)
        #cancel = Tkinter.Button(optionframe, text="cancel", width=10, command=cancel)
        #apply.pack(in_ = optionframe, side=Tkinter.LEFT)
        #cancel.pack(in_ = optionframe, side=Tkinter.LEFT)
        #apply = Tkinter.Button(histogram_window, text="apply", width=10, command=apply)
        #cancel = Tkinter.Button(histogram_window, text="cancel", width=10, command=cancel)
        #apply.grid(row=0, column=0, sticky=Tkinter.W)
        #cancel.grid(row=0, column=1, sticky=Tkinter.W)

        #hmenubar.add_command(label="apply", command=apply)
        #hmenubar.add_command(label="cancel", command=cancel)
        #histogram_window.config(menu=hmenubar)
        click(mv=1)
        histogram_window.mainloop()

    def save_gs_buffer(self, slot):
        #save primary greyscale buffer (gimg) in greyscale buffer number (slot - int)
        #print slot
        if self.gimg.mode=="L":
            self.greyscale_buffers[slot] = self.gimg.copy()
            print "# Saved in greyscale slot " + str(slot)
            if logger:
                print "gimg" + str(slot)  + " = gimg.copy()"
        else:
            print "# Not a greyscale image"

    def load_gs_buffer(self, slot):
        #load greyscale buffer number (slot - int) into primary greyscale buffer (gimg)
        if self.greyscale_buffers[slot]!=None:
            self.gimg = self.greyscale_buffers[slot]
            print "# Loaded greyscale slot " + str(slot)
            if logger:
                self.logging("gimg = gimg" + str(slot))
            self.update()
        else:
            print "# greyscale buffer " + str(slot) + " is empty"

    def save_bin_buffer(self, slot):
        #save primary binary buffer (bimg) in binary buffer number (slot - int)
        if self.bimg.mode=="1":
            self.bin_buffers[slot] = self.bimg.copy()
            print "# Saved in greyscale slot " + str(slot)
            if logger:
                self.logging("binimg" + str(slot)  + " = bimg.copy()")
        else:
            print "# Not a binary image"

    def load_bin_buffer(self, slot):
        #load binary buffer number (slot - int) into primary binary buffer (bimg)
        if self.bin_buffers[slot]!=None:
            self.bimg = self.bin_buffers[slot]
            print "# Loaded binary slot " + str(slot)
            if logger:
                self.logging("bimg = binimg" + str(slot))
            self.update()
        else:
            print "# binary buffer " + str(slot) + " is empty"

    def edgedetect_find_edges(self):
        self.filter_function("FIND_EDGES")

    def sobel(self):
        tmpdata = scipy.ndimage.sobel(self.gimg, 0)
        tmpimg = PIL.Image.new(mode="L", size=self.gimg.size)
        tmpimg.putdata(tmpdata.flatten())
        self.gimg = tmpimg
        self.update()
    
    def filter_function(self, filtername):
        self.gimg = self.gimg.filter(filter(filtername))
        if logger:
            self.logging("gimg = self.gimg.filter(PIL.Filters." + filtername + ")")

        self.update()

    def run_grid(self, fn):
        # runs the function with n by n grid, n supplied by user
        top = dialog_window(self.master)
        e = Tkinter.Entry(top)
        e.insert(Tkinter.END, "5")
        e.grid(row=2, column=0)
        e.focus_set()
        def run(n): #this can probably be simplified
            fn(grid(self.gimg, n))
            top.destroy()
        applyfn = lambda : run(int(e.get()))
        top.setapply(applyfn)
        #apply = Tkinter.Button(top, text="apply", width=10, command=applyfn)
        #apply.pack()
        
    def watershed(self, point_list):
        if logger:
            self.logging("points = " + str(point_list))
            self.logging("data = watershed_pts(points)")
        self.gimg = watershed_pts(self.gimg, point_list)
        self.update()

    def lut(self):
        top = dialog_window(self.master)
        old = self.gimg.copy()
        self.fn=None
        def transform(method):
            if method=="sqrt":
                lutfn = lambda x: np.sqrt(x)*16
            if method=="log":
                lutfn = lambda x: np.log(x)*46
            if method=="square":
                lutfn = lambda x: x**2/255
            self.fn = method
            self.gimg = old.point(lutfn)
            self.update(False)

        def cancel(event=None):
            self.gimg = old
            top.destroy()
            self.update(False)

        def apply(event=None):
            if self.fn!=None:
                if self.fn=="sqrt":
                    self.logging("gimg = gimg.point(lambda x: np.sqrt(x)*16)")
                if self.fn=="log":
                    self.logging("gimg = gimg.point(lambda x: np.log(x)*46)")
                if self.fn=="square":
                    self.logging("gimg = gimg.point(lambda x: x**2/255)")
            top.destroy()
            self.update()
            
        #menubar = Tkinter.Menu(self.master)
        #menubar.add_command(label="apply", command=apply)
        #menubar.add_command(label="cancel", command=cancel)
        #top.config(menu=menubar)
        sqrtbtn = Tkinter.Button(top, text="sqrt", width=10, 
                command=lambda : transform("sqrt"))
        logbtn = Tkinter.Button(top, text="log", width=10, 
                command=lambda : transform("log"))
        squarebtn = Tkinter.Button(top, text="square", width=10, 
                command=lambda : transform("square"))
        #sqrtbtn.pack()
        #logbtn.pack()
        #squarebtn.pack()
        sqrtbtn.grid(row=2, column=0)
        logbtn.grid(row=3, column=0)
        squarebtn.grid(row=4, column=0)

        top.setcancel(cancel)
        top.setapply(apply)
            

    def invert_bin(self):
        self.show_binary=1 
        if self.bimg.mode=="1":
            self.bimg = PIL.ImageChops.invert(self.bimg)
            self.update()
        else:
            print "#Not a binary image"

    def picker(self, action):
        gold = self.gimg.copy()
        bold = self.bimg.copy()
        mao, n = scipy.ndimage.measurements.label(np.array(
            self.bimg.convert("L")))

        if action=="keep":
            self.picked = PIL.Image.new(mode="1", size=self.img.size)
        elif action=="reject":
            self.picked = self.bimg.copy()
            self.bimg = self.bimg.point(lambda x: 0)

        self.pick = PIL.Image.new(mode="1", size=self.gimg.size)
        self.gimg = self.gimg.convert("RGBA")
        def pick(event=None):
            ma = np.ndarray.copy(mao)
            tmp = self.bimg.copy()
            if event!=None:
                x = int(event.x/self.resize_ratio)
                y = int(event.y/self.resize_ratio)
                if ma[y,x]!=0:
                    ma[ma!=ma[y,x]] = 0
                    ma[ma==ma[y,x]] = 1
                else:
                    return
            else:
                ma[ma!=-1] = 0
            self.pick.putdata(ma.flatten())
            self.picked = PIL.ImageChops.logical_xor(self.picked, self.pick)
            #self.picked.convert("L").point(lambda x: 255*x).show()
            self.bimg = PIL.ImageChops.logical_xor(self.bimg, self.pick)

            colors = PIL.Image.new(mode="RGBA", size=self.gimg.size).split()
            colors[1].paste(self.picked.convert("L").point(
                lambda x: 255*x))
            binmask = PIL.Image.merge("RGBA", colors)
            binmask.putalpha(colors[1])
            #binmask.show() 
            #self.picked.convert("L").point(lambda x: 255*x).show()
            self.gimg.paste(binmask, (0,0), self.picked)
            self.update(False)


        def apply(event=None):
            self.bimg = self.picked
            self.master.unbind("<ButtonPress-1>")
            self.gimg = gold
            top.destroy()
            self.update()
        
        def cancel(event=None):
            self.bimg = bold #These names are confusing
            self.master.unbind("<ButtonPress-1>")
            self.gimg = gold
            top.destroy()
            self.update()

        self.master.bind("<ButtonPress-1>", pick)
        top = dialog_window(self.master)
        top.setcancel(cancel)
        top.setapply(apply)
        pick()
        

    def logical_operators(self):
        # Performs the chosen binary operation on the selected binary 
        # buffers (def 1 and 2) and places the result in primary buffer (bimg)
        self.show_binary=1 
        top = dialog_window(self.master)
        e1 = Tkinter.Entry(top)
        e1.insert(Tkinter.END, "0")
        #e1.pack()
        e1.grid(row=2, column=0)
        e1.focus_set()
        e2 = Tkinter.Entry(top)
        e2.insert(Tkinter.END, "1")
        #e2.pack()
        e2.grid(row=2, column=1)

        def logic(op):
            n,m = (int(e1.get()),int(e2.get()))
            if self.bin_buffers[n]!=None and self.bin_buffers[m]!=None:
                if op=="and":
                    self.bimg = PIL.ImageChops.logical_and(self.bin_buffers[n], 
                            self.bin_buffers[m])
                elif op=="or":
                    self.bimg = PIL.ImageChops.logical_or(self.bin_buffers[n], 
                            self.bin_buffers[m])
                elif op=="xor":
                    self.bimg = PIL.ImageChops.logical_xor(self.bin_buffers[n], 
                            self.bin_buffers[m])
                if logger:
                    self.logging(("bimg = PIL.ImageChops.logical_" + op + "(binimg" +
                            str(n) + ", binimg" + str(m) + ")"))
                top.destroy()
                self.update()
            elif self.bin_buffers[n]==None:
                print "# binary buffer " + str(n) + " is empty"
            elif self.bin_buffers[m]==None:
                print "# binary buffer " + str(m) + " is empty"
                

        l_and= Tkinter.Button(top, text="AND", width=10, 
                command=lambda : logic("and"))
        l_or = Tkinter.Button(top, text="OR", width=10, 
                command=lambda : logic("or"))
        l_xor = Tkinter.Button(top, text="XOR", width=10, 
                command=lambda : logic("xor"))
        #l_and.pack()
        #l_or.pack()
        #l_xor.pack()
        l_and.grid(row=3, column=0)
        l_or.grid(row=4, column=0)
        l_xor.grid(row=5, column=0)


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
            self.logging('tmp = PIL.Image.new(mode="1", size = img.size)')
            self.logging(('tmp = scipy.ndimage.binary_' +
            morph_type + '(img.convert("L"), structure=unit_disk(' +
                    str(count) + "))"))
            self.logging("tmp.putdata(255*tmp.flatten())")
            self.logging('bimg = tmp.convert("1")')

        img2.putdata(255*tmp.flatten())
        top.destroy()
        self.bimg = img2.convert("1")
        self.update()

    def morph_options(self, morph_type):
        top = dialog_window(self.master)
        e = Tkinter.Entry(top)
        e.insert(Tkinter.END, "2")
        #e.pack()
        e.grid(row=2, column=0)
        e.focus_set()
        applyfn = lambda : self.morph(int(e.get()), top, morph_type)
        top.setapply(applyfn)
        #apply = Tkinter.Button(top, text="apply", width=10, command= applyfn)
        #apply.pack()

    def fill_option(self):
        #maybe add options later though not strictly necessary
        top = dialog_window(self.master)
        self.morph(0,top,"fill_holes")

    def shade(self):
        self.gimg = shade_by_size(self.bimg)
        self.toggle_binary()
        self.update()
        if logger:
            self.logging("gimg = shade_by_size(bimg)")

    def lin_divide_bins(self):
        # Plot or output? Maybe both? Decide later
        # Should there be a log?
        count, bins = linbin(self.bimg)
        plt.bar(range(len(bins)), count, width=1)
        plt.ylabel("Count")
        if(self.conversion_ratio):
            plt.xlabel("Size (microns^2)")
            plt.xticks([z for z in range(0,len(bins))],
                    ["<"+'%s' % float('%.2g' % (b*self.conversion_ratio**2))
                        for b in bins])
        else:
            plt.xlabel("Size (pixels)")
            plt.xticks([z for z in range(0,len(bins))],
                ["<"+str(int(b)) for b in bins])
        plt.show()

    def log_divide_bins(self):
        # Plot or output? Maybe both? Decide later
        # Should there be a log?
        count, bins = linbin(self.bimg)
        plt.bar(range(len(bins)), count, width=1)

        plt.ylabel("Count")
        if(self.conversion_ratio):
            plt.xlabel("Size (microns^2)")
            plt.xticks([z for z in range(0,len(bins))],
                    ["<"+'%s' % float('%.2g' % (b*self.conversion_ratio**2)) 
                        for b in bins])
        else:
            plt.xlabel("Size (pixels)")
            plt.xticks([z for z in range(0,len(bins))],
                ["<"+str(int(b)) for b in bins])
        plt.show()

    def set_conversion_ratio(self):
        def apply(event=None):
            self.conversion_ratio = float(e.get())
            self.update()
            conv_window.destroy()
        conv_window = dialog_window()
        l = Tkinter.Label(conv_window, text="Microns per pixel: ")
        e = Tkinter.Entry(conv_window)
        e.insert(Tkinter.END, str(self.conversion_ratio))
        l.grid(row=2, column=0)
        e.grid(row=2, column=1)
        conv_window.setapply(apply)

    def greyscale_arithmetic(self):
        # TODO: Add exception for missing gimg
        #       Add logger
        #       
        g_window = dialog_window()
        fn_string = Tkinter.StringVar()
        fn_string.set("a")
        abuff_label = Tkinter.Label(g_window, text="input a (gimg)")
        abuff_entry= Tkinter.Entry(g_window)
        abuff_entry.insert(0, "0")
        bbuff_label = Tkinter.Label(g_window, text="input b (gimg)")
        bbuff_entry= Tkinter.Entry(g_window)
        bbuff_entry.insert(0, "1")
        fn_entry = Tkinter.Entry(g_window, textvariable=fn_string)
        old = self.gimg.copy()

        def apply(event=None):
            self.gimg = (PIL.ImageMath.eval("convert(" + fn_entry.get() + ",'L')", 
                    a=self.greyscale_buffers[int(abuff_entry.get())], 
                    b=self.greyscale_buffers[int(bbuff_entry.get())]))
            if logger:
                self.logging("gimg = PIL.ImageMath.eval('convert("+fn_entry.get()+',"L"'+")',"+ 
                        "a=gimg" + abuff_entry.get()+", b=gimg" + bbuff_entry.get() + ")")
            self.update()
            g_window.destroy()

        def cancel(event=None):
            self.gimg = old
            g_window.destroy()
            self.update(False)

        def update(event=None):
            self.gimg = (PIL.ImageMath.eval("convert(" + fn_entry.get() + ",'L')", 
                    a=self.greyscale_buffers[int(abuff_entry.get())], 
                    b=self.greyscale_buffers[int(bbuff_entry.get())]))
            self.update(False)
        def minval(event=None):
            fn_string.set("min(a,b)")
            update()

        def maxval(event=None):
            fn_string.set("max(a,b)")
            update()

        def add(event=None):
            fn_string.set("a+b")
            update()

        def subtract(event=None):
            fn_string.set("a-b")
            update()

        updatebutton    = Tkinter.Button(g_window, text = "Update", width = 10,
                command = update)
        minbutton       = Tkinter.Button(g_window, text = "Min", width = 10,
                command = minval)
        maxbutton       = Tkinter.Button(g_window, text = "Max", width = 10,
                command = maxval)
        addbutton       = Tkinter.Button(g_window, text = "Add", width = 10,
                command = add)
        subtractbutton  = Tkinter.Button(g_window, text = "Subtract", width = 10,
                command = subtract)
        fn_entry.grid(row=2, column=0)
        updatebutton.grid(row=2, column=1)
        minbutton.grid(row=3, column=0)
        maxbutton.grid(row=4, column=0)
        addbutton.grid(row=5, column=0)
        subtractbutton.grid(row=6, column=0)
        abuff_label.grid(row=3, column=1)
        abuff_entry.grid(row=4, column=1)
        bbuff_label.grid(row=5, column=1)
        bbuff_entry.grid(row=6, column=1)
        g_window.setapply(apply)
        g_window.setcancel(cancel)

    def apply_binary_mask(self):
        self.gimg = PIL.ImageMath.eval("convert(min(a,b),'L')", 
                a=self.gimg, b=self.bimg.convert('L'))
        if logger:
            self.logging('PIL.ImageMath.eval("convert(min(a,b),'+"'L')"+'"'+", a=self.gimg, b=self.bimg.convert('L'))")
        self.update()

    def fft(self):
        #Add logging
        #if self.phase!=None:
        #print "# Phase data being overwritten"
        tmpimg = PIL.Image.new(mode="L", size=self.gimg.size)
        fftvals = np.fft.fftshift(np.fft.fft2(self.gimg))
        freq = np.abs(fftvals)
        self.phase = np.angle(fftvals)
        freq = np.log(freq+np.e)
        self.fftmax = np.max(freq)
        freq = freq*255./self.fftmax
        tmpimg.putdata(freq.flatten())
        self.gimg = tmpimg
        self.update()
        #self.logging("fft = np.fft.fftshift(np.fft.fft2(gimg))")
        #self.logging("freq = np.abs(fftvals)")
        #self.logging("phase = np.angle(fftvals)")


    def ifft(self):
        #Add logging
        if self.phase==None:
            print "# No phase data available"
            return
        tmpimg = PIL.Image.new(mode="L", size=self.gimg.size)
        fftvals = np.array(self.gimg.getdata()).reshape((self.gimg.size[1],self.gimg.size[0]))

        #fftvals = (1-2*(fftvals>=128))*(np.exp((np.abs((1-2*(fftvals>=128))*fftvals+256*(fftvals>=128)))*self.fftmax/126.)-np.e)
        fftvals = ((np.exp(fftvals*self.fftmax/255))-np.e)
        ifftvals = np.fft.ifft2(np.fft.ifftshift(fftvals*np.exp(1j*self.phase)))
        tmpimg.putdata(ifftvals.flatten())
        self.gimg = tmpimg
        self.phase = None
        self.update()
    
    def draw_center_circle(self):
        circ_dialog = dialog_window()
        orig = self.bimg.copy()
        size_entry = Tkinter.Entry(circ_dialog)
        size_entry.insert(Tkinter.END, "1")
        def apply(event=None):
            draw = PIL.ImageDraw.Draw(self.bimg)
            r = int(size_entry.get())
            draw.ellipse([self.bimg.size[0]/2-r, self.bimg.size[1]/2-r, 
                self.bimg.size[0]/2+r, self.bimg.size[1]/2+r], fill=255)
            self.update()
            circ_dialog.destroy()
            if logger:
                self.logging("draw.ellipse(["+str(self.bimg.size[0]/2-r)+", "+ str(self.bimg.size[1]/2-r)+", "+
                    str(self.bimg.size[0]/2+r)+", "+str(self.bimg.size[1]/2+r)+"], fill=255)")
        def cancel(event=None):
            self.bimg = orig
            circ_dialog.destroy()
            self.update()
        size_entry.grid(row=2,column=0)
        circ_dialog.setapply(apply)
        circ_dialog.setcancel(cancel)
            
    def logging(self, logstr, func=None):
        print logstr

        
root = Tkinter.Tk()
root.geometry('+%d+%d' % (100,100))
filename=None
if len(sys.argv)>1:
    filename=str(sys.argv[1])
gui = Imagewindow(root, filename)
