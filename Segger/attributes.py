# -----------------------------------------------------------------------------
# Allow associating number, string and image attributes with segmentation
# regions, saving them segmentation files, browsing them in a dialog, and
# searching for regions with specific attribute values.
#

import regions
from chimera.baseDialog import ModelessDialog

import numpy
int_types = (int, numpy.int8, numpy.int16, numpy.int32, numpy.int64,
             numpy.uint8, numpy.uint16, numpy.uint32, numpy.uint64)
float_types = (float, numpy.float32, numpy.float64)

class Attribute_Dialog(ModelessDialog):

  title = 'Segmentation Region Attributes'
  name = 'segmentation region attributes'
  buttons = ('Update', 'Options', 'Close',)
  help = 'ContributedSoftware/segger/segment.html#attributes'
  
  def fillInUI(self, parent):

    self.sel_handler = None
    self.last_selected = set()
    
    import Tkinter
    from CGLtk import Hybrid

    top = parent.winfo_toplevel()
    menubar = Tkinter.Menu(top, type="menubar", tearoff=False)
    top.config(menu=menubar)

    file_menu_entries = (('Export...', self.export_cb),)
    fmenu = Hybrid.cascade_menu(menubar, 'File', file_menu_entries)

    self.columnMenu = Tkinter.Menu(menubar)
    menubar.add_cascade(label="Columns", menu=self.columnMenu)

    from chimera.tkgui import aquaMenuBar
    aquaMenuBar(menubar, parent, row = 0)

    row = 1


    columns = (('region', lambda r: r.rid),
               ('grid points', lambda r: r.point_count()),
               ('grouped', lambda r: len(r.children())),
               ('has surface', lambda r: 1 if r.has_surface() else 0),
               ('contacts', lambda r: len(r.contacting_regions())),
               ('edge distance', lambda r: r.edge_distance()),
               ('bounds', lambda r: tuple(r.bounds()[0])+tuple(r.bounds()[1]),
                '%d,%d,%d,%d,%d,%d', False),
               )
    self.computed_attributes = dict([c[:2] for c in columns])
    self.keys = set([c[0] for c in columns])
    defaults = dict([(c[0],True) for c in columns])
    prefs = {}
    mi = (self.columnMenu, prefs, defaults, True)
    from CGLtk.Table import SortableTable
    self.attribute_table = t = SortableTable(parent, menuInfo=mi)
    for cspec in columns:
      name, func = cspec[:2]
      format = cspec[2] if len(cspec) >= 3 else format_item
      display = cspec[3] if len(cspec) >= 4 else True
      t.addColumn(name, func, format = format, display = display,
                  font='TkFixedFont', refresh = False)
    self.add_columns()  # Additional region attributes.
    t.setData([])
    t.launch(browseCmd=self.table_line_selected)

    t.grid(row=row, column=0, sticky="nsew")
    parent.rowconfigure(row, weight=1)
    parent.columnconfigure(0, weight=1)
    row += 1

    naf = Tkinter.Frame(parent)
    naf.grid(row=row, column=0, sticky = 'new')
    naf.columnconfigure(1, weight=1)
    row += 1

    ke = Hybrid.Entry(naf, 'Set attribute ', 15)
    ke.frame.grid(row = 0, column = 0, sticky = 'w')
    self.key_name = ke.variable
    ke.entry.bind('<KeyPress-Return>', self.new_value_cb)
    
    ve = Hybrid.Entry(naf, ' to value ', 15)
    ve.frame.grid(row = 0, column = 1, sticky = 'ew')
    self.value = ve.variable
    ve.entry.bind('<KeyPress-Return>', self.new_value_cb)

    ol = Tkinter.Label(naf, text = ' or ')
    ol.grid(row = 0, column = 2, sticky = 'e')

    sb = Tkinter.Button(naf, text = 'snapshot', command = self.snapshot_cb)
    sb.grid(row = 0, column = 3, sticky = 'e')

    fl = Hybrid.Checkbutton_Entries(parent, False, 'Filter list ', (50, 'grid_points > 1000'))
    fl.frame.grid(row = row, column = 0, sticky = 'new')
    row += 1
    self.use_filter, self.filter_text = fl.variables
    self.use_filter.add_callback(self.filter_cb)
    e = fl.entries[0]
    fl.frame.columnconfigure(1, weight=1)
    e.grid(sticky = 'ew')
    e.bind('<KeyPress-Return>', self.filter_cb)

    op = Hybrid.Popup_Panel(parent, resize_dialog = False)
    opf = op.frame
    opf.grid(row = row, column = 0, sticky = 'news')
    opf.grid_remove()
    opf.columnconfigure(0, weight=1)
    self.optionsPanel = op.panel_shown_variable
    row += 1
    orow = 0

    cb = op.make_close_button(opf)
    cb.grid(row = orow, column = 0, sticky = 'e')

    l = Tkinter.Label(opf, text='Options', font = 'TkCaptionFont')
    l.grid(column=0, row=orow, sticky='w', pady=5)
    orow += 1

    ss = Hybrid.Checkbutton(opf, 'Show selected region surface', False)
    ss.button.grid(row = orow, column = 0, sticky = 'w')
    orow += 1
    self.show_selected = ss.variable

    msg = Tkinter.Label(parent, width = 60, anchor = 'w', justify = 'left')
    msg.grid(column=0, row=row, sticky='ew')
    self.msg = msg

    self.Update()

  def map(self, event = None):

    from chimera import triggers
    self.sel_handler = triggers.addHandler('selection changed',
                                           self.selection_cb, None)

  def unmap(self, event = None):

    from chimera import triggers
    triggers.deleteHandler('selection changed', self.sel_handler)
    self.sel_handler = None

  def status(self, message):

    text = message.rstrip('\n')
    self.msg.configure(text = text)
    self.msg.update_idletasks()

  def Options(self):

    self.optionsPanel.set(not self.optionsPanel.get())

  def export_cb(self):

    def save ( okay, dialog, self=self ):
      if okay:
        paths = dialog.getPaths ( )
        if paths:
          self.write_comma_separated_values(paths[0])

    from segment_dialog import current_segmentation
    seg = current_segmentation()
    if seg is None:
      idir = ifile = None
    elif hasattr(seg, 'path'):
      import os.path
      idir, ifile = os.path.split(seg.path)
      ifile = os.path.splitext(ifile)[0] + '.csv'
    else:
      import os.path
      idir = None
      ifile = os.path.splitext(seg.name)[0] + '.csv'
      
    from OpenSave import SaveModeless
    SaveModeless ( title = 'Export Comma Separated Values',
                 initialdir = idir, initialfile = ifile, command = save )

  def write_comma_separated_values(self, path):

    f = open(path, 'w')
    t = self.attribute_table
    col = [c for c in t.columns if c.display]
    f.write(csv_line([c.title for c in col]))
    for r in t._sortedData():
      fields = [c.displayValue(r) for c in col]
      # Don't export images or other non-string values.
      fields = [v if isinstance(v,basestring) else '' for v in fields]
      f.write(csv_line(fields))
    f.close()

  def Update(self):

    from segment_dialog import current_segmentation
    s = current_segmentation()
    if s is None:
      return

    rlist = self.filter_regions(s.regions)
    if rlist is None:
      return
      
    self.status('Table has %d regions from %s' % (len(rlist), s.name))

    self.add_columns()
    self.attribute_table.setData(rlist)

  def add_columns(self):
    
    from segment_dialog import current_segmentation
    s = current_segmentation()
    if s is None:
      return

    k = set(self.computed_attributes.keys())
    for r in s.regions:
      k.update(r.attributes().keys())

    knew = list(k.difference(self.keys))
    kgone = self.keys.difference(k)
    if knew or kgone:
      t = self.attribute_table
      mreg = t.middleRow()
      knew.sort()
      for key in knew:
        self.add_column(key, refresh = False)
      for key in kgone:
        t.removeColumn(key, refresh = False)
      self.keys = k
      t.refresh(rebuild = True)
      t.showRow(mreg)   # TODO: Adding column loses scroll position.

  def add_column(self, key, refresh = True):

    if key in self.keys:
      return

    def get_val(region, key = key):
      v = region.get_attribute(key)

      from PIL.Image import Image
      if isinstance(v, Image):
        v = self.image_thumbnail(v, region, key)
      return v

    t = self.attribute_table
    t.addColumn(key, get_val, format = format_item, font = 'TkFixedFont',
                refresh = refresh)
    self.keys.add(key)

  def image_thumbnail(self, image, region, key, size = (32,32)):

    if hasattr(image, 'thumbnail_table_image'):
      return image.thumbnail_table_image

    ti = image.copy()
    from PIL.Image import ANTIALIAS
    ti.thumbnail(size, ANTIALIAS)
    from PIL.ImageTk import PhotoImage
    image.thumbnail_table_image = pi = PhotoImage(ti)
    def show_image(event, r=region, k=key, i=image, self=self):
      from imageviewer import showImage
      shift_mask = 1
      showImage(i, title = '%s region %d' % (k, r.rid),
                new_window = (event.state & shift_mask))
      self.table_line_selected([r])
    pi.button_press_callback = show_image
    return pi

  def selection_cb(self, trigger, unused, selection):

    rlist = regions.SelectedRegions()

    t = self.attribute_table
    if rlist:
      t.select(rlist)

    rset = set(rlist)
    if rlist and rset != self.last_selected:
      t.showRow(rlist[0])
    self.last_selected = rset

  def table_line_selected(self, rsel):

    rsel = [r for r in rsel if not r.segmentation.__destroyed__]
    regions.select_regions(rsel)
    self.last_selected = set(rsel)

    self.show_attribute_value(rsel)

    if rsel and self.show_selected.get():
      regions.show_only_regions(rsel)

  def show_attribute_value(self, rsel):

    key = self.key_name.get().strip()
    if len(key) == 0:
      return

    rlist = [r for r in rsel if r.has_attribute(key)]
    if len(rlist) == 0:
      return

    val = rlist[0].get_attribute(key)
    if isinstance(val, (basestring, int, float)):
      self.value.set(str(val))

  def new_value_cb(self, event = None):

    # Try interpreting value as int or float, otherwise string
    val = self.value.get().strip()
    try:
      val = int(val)
    except ValueError:
      try:
        val = float(val)
      except ValueError:
        pass

    self.set_attribute(val)

  def set_attribute(self, val):

    key = self.key_name.get().strip()
    if len(key) == 0:
      return

    t = self.attribute_table
    rlist = t.selected()
    if len(rlist) == 0:
      self.status('No list rows selected')
      return

    if rlist[0].is_reserved_name(key):
      self.status('%s is a reserved attribute name' % key)
      return

    if isinstance(val, basestring) and len(val) == 0:
      for r in rlist:
        r.remove_attribute(key)
    else:
      for r in rlist:
        r.set_attribute(key, val)

    if not key in self.keys:
      mreg = t.middleRow()
      self.add_column(key)
      t.showRow(mreg)   # TODO: Adding column loses scroller position.
    else:
      t.refresh(selection = rlist)

  def snapshot_cb(self):

    import chimera
    width,height = chimera.viewer.windowSize
    image = chimera.viewer.pilImages(width, height, supersample = 3)[0]
    self.set_attribute(image)
    
  def filter_cb(self, event=None):

    self.Update()

  def filter_regions(self, regions):

    if not self.use_filter.get():
      return list(regions)

    ftext = self.filter_text.get().strip()
    if len(ftext) == 0:
      return list(regions)

    rlist = []

    # Make globals dictionary for evaluating filter that has all attributes
    # set to None.
    n2n = dict([(a, pyvar_name(a)) for a in self.keys])
    glb = dict([(n2n[a],None) for a in self.keys])
    glb['__builtins__'] = __builtins__
    glb.update(__builtins__)
    
    try:
      for r in regions:
        loc = dict([(n2n[a],v) for a,v in r.attributes().items()])
        loc.update([(n2n[a],f(r)) for a,f in self.computed_attributes.items()])
        if eval(ftext, glb, loc):
          rlist.append(r)
    except Exception as e:
      self.status('Error in filter expression: %s' % e)
      return None

    return rlist

  def get_region_attribute(self, r, name):

    ca = self.computed_attributes
    if name in ca:
      return ca[name](r)
    elif r.has_attribute(name):
      return r.get_attribute(name)
    raise NoAttributeValue, (r,name)

  def has_region_attribute(self, r, name):

    return name in self.computed_attributes or r.has_attribute(name)

# -----------------------------------------------------------------------------
#
def format_item(value):

  from PIL.ImageTk import PhotoImage
  if isinstance(value, basestring):
    return value
  elif isinstance(value, int_types):
    return '%d' % value
  elif isinstance(value, float_types):
    return '%.5g' % value
  elif value is None:
    return ''
  elif isinstance(value, PhotoImage):
    return value
  return str(value)
  
# -----------------------------------------------------------------------------
# Convert name to legal Python variable, eliminating spaces and syntax
# characters like '.'
# 
def pyvar_name(name):

  n = ''.join([(c if c.isalnum() else '_') for c in name])
  if len(n) == 0:
    n = 'blank'
  elif n[0].isdigit():
    n = '_' + n
  return n

# -----------------------------------------------------------------------------
#
def csv_line(fields):

  return ','.join([csv_escape(f) for f in fields]) + '\n'

# -----------------------------------------------------------------------------
#
def csv_escape(f):

  fq = f.replace('"', '""')
  if ',' in fq or '\n' in fq:
    fq = '"' + fq + '"'
  return fq
  
# -----------------------------------------------------------------------------
#
def region_attributes_dialog(create = False):

  from chimera import dialogs
  return dialogs.find ( Attribute_Dialog.name, create=create )


# -----------------------------------------------------------------------------
#
def show_region_attributes_dialog ():

    from chimera import dialogs
    d = region_attributes_dialog(create = True)
    # Avoid transient dialog resizing when created and mapped for first time.
#    d.toplevel_widget.update_idletasks ()
    d.enter()

# -----------------------------------------------------------------------------
#
from chimera import dialogs
dialogs.register (Attribute_Dialog.name, Attribute_Dialog, replace = True)

