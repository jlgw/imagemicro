import Tkinter

class dialog_window(Tkinter.Toplevel):
    def __init__(self, topwindow=None):
        Tkinter.Toplevel.__init__(self,topwindow)
        self.wm_transient(topwindow) #Keep dialog buttons on top
        self.bind("<Return>", self.apply)
        self.bind("<Escape>", self.cancel)
        self.protocol("WM_DELETE_WINDOW", self.cancel) #close by x button does the same as cancel
        self.wm_title("")
        self.applybtn = Tkinter.Button(self, text="apply", width=10, 
                command=self.apply)
        self.cancelbtn = Tkinter.Button(self, text="cancel", width=10, 
                command=self.cancel)
        self.sep = Tkinter.Frame(self, height=5, bg="grey")
        
        self.applybtn.grid(row=0, column=0, sticky=Tkinter.W)
        self.cancelbtn.grid(row=0, column=1, sticky=Tkinter.W)
        self.sep.grid(row=1, columnspan=2, sticky="ew")

    def cancel(self, event=None):
        self.destroy()
    def apply(self, event=None):
        pass

    #fn can't have any arguments, usually not needed but 
    #may need to be changed for more complex stuff
    def setcancel(self, fn):
        self.cancel= lambda x=None: fn()
        # Maybe these duplicate calls can be replaced by a second call to
        # __init__?
        self.bind("<Escape>", self.cancel)
        self.protocol("WM_DELETE_WINDOW", self.cancel)
        self.cancelbtn.configure(command=self.cancel)

    def setapply(self, fn):
        self.apply = lambda x=None: fn()
        self.bind("<Return>", self.apply)
        self.applybtn.configure(command=self.apply)
