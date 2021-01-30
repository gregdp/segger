# ------------------------------------------------------------------------------
# Display an image in a window.
#
def showImage(image, title = '',
              new_window = False, resize_window = False, fit_to_window = False):

    if new_window:
        v = ImageViewer()
        resize_window = True
    else:
        from chimera import dialogs
        v = dialogs.find(ImageViewer.name, create=False)
        if v is None:
            v = dialogs.find(ImageViewer.name, create=True)
            resize_window = True
    v.set_image(image, resize_window, fit_to_window)
    if title:
        v.set_title(title)
    v.enter()

# ------------------------------------------------------------------------------
#
from chimera.baseDialog import ModelessDialog
class ImageViewer(ModelessDialog):
    title = "Image Viewer"
    name = "image viewer"
    buttons = ("Fit", "Full", "Window", "Save", "Close")
    help = 'UsersGuide/midas/segment.html#sliceimage-window'

    def fillInUI(self, parent):

        self.image = None
        self.shown_size = None
        self.image_tag = None

        import Tkinter

        parent.rowconfigure(0, weight=1)
        parent.columnconfigure(0, weight=1)

        self.canvas = c = Tkinter.Canvas(parent, highlightthickness = 0)
        c.grid(row = 0, column = 0, sticky = 'news')

        sv = Tkinter.Scrollbar(parent, orient=Tkinter.VERTICAL, command=c.yview)
        sv.grid(row = 0, column = 1, sticky = 'ns')
        sh = Tkinter.Scrollbar(parent, orient=Tkinter.HORIZONTAL,
                               command=c.xview)
        sh.grid(row = 1, column = 0, sticky = 'ew')

        c.config(yscrollcommand=sv.set)
        c.config(xscrollcommand=sh.set)

    def set_title(self, title):

        self._toplevel.title(title)

    def set_image(self, image, resize_window = False, fit_to_window = False):

        self.image = image
        if fit_to_window and not resize_window:
            self.Fit()
        else:
            self.set_canvas_image(image, resize_window)

    def set_canvas_image(self, image, resize_window = False):

        from PIL.ImageTk import PhotoImage
        pi = PhotoImage(image)
        c = self.canvas
        c.save_image = pi	# Tk Label keeps no count.
        t = c.create_image(0,0,anchor="nw",image=pi)
        if not self.image_tag is None:
            c.delete(self.image_tag)
        self.image_tag = t
        self.shown_size = w,h = image.size
        c.config(scrollregion=(0,0,w,h))

        if resize_window:
            self.Window()

    def Window(self):

        if self.shown_size is None:
            return

        # Resize window to fit image.
        w,h = self.shown_size
        self.canvas.config(width = w, height = h)
        self.uiMaster().winfo_toplevel().geometry('')

    def Fit(self):

        if self.image is None:
            return

        c = self.canvas
        w,h = c.winfo_width(), c.winfo_height()

        # Preserve aspect ratio.
        wi,hi = self.image.size
        if h*wi > w*hi:
            h = int((float(w)/wi)*hi)
        else:
            w = int((float(h)/hi)*wi)

        if w <= 0 or h <= 0:
            return

        from PIL import Image
        ri = self.image.resize((w,h), Image.ANTIALIAS)
        self.set_canvas_image(ri)

    def Full(self):

        if self.image is None:
            return

        self.set_canvas_image(self.image)

    def Save(self):

        if self.image is None:
            return

	def save(okay, dialog, self = self):
            if not okay:
                return
            paths_and_types = dialog.getPathsAndTypes()
            if len(paths_and_types) == 0:
                return
            path, format = paths_and_types[0]
            self.image.save(path, format)

	from OpenSave import SaveModeless
	SaveModeless(title = 'Save Image',
		     filters = [('JPEG', '*.jpg', '.jpg'),
				('PNG', '*.png', '.png'),
				('TIFF', '*.tif', '.tif')],
		     command = save )

# ------------------------------------------------------------------------------
#
from chimera import dialogs
dialogs.register(ImageViewer.name, ImageViewer, replace = True)
