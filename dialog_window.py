import Tkinter

class dialog_window(Tkinter.Toplevel):
    def __init__(self):
        Tkinter.Toplevel.__init__(self)
        self.bind("<Return>", self.apply)
        self.bind("<Escape>", self.cancel)
        self.protocol("WM_DELETE_WINDOW", self.cancel) #close by x button does the same as cancel

    def cancel(self, event=None):
        self.destroy()
    def apply(self, event=None):
        pass

    #fn can't have any arguments, usually not needed but 
    #may need to be changed for more complex stuff
    def setcancel(self, fn):
        self.cancel= lambda x: fn()
        self.bind("<Escape>", self.cancel)
        self.protocol("WM_DELETE_WINDOW", self.cancel) #close by x button does the same as cancel

    def setapply(self, fn):
        self.apply = lambda x: fn()
        self.bind("<Return>", self.apply)
