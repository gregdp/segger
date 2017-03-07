# -----------------------------------------------------------------------------
# Read and write segmentation data in hdf5 format.
#
# Example layout:
#
# format = "segger"
# format_version = 1
#
# name = "somedata segmentation"
# mask = <3-d array of region indices>
# region_ids = <array or region ids numbers, length N>
# region_colors = <array of rgba, N by 4>
# ref_points = <array of region reference points, N by 3>
# parent_ids = <array each regions parent (0 = no parent), length N>
# smoothing_levels = <array of float, length N>
#
# map_path = "/Users/smith/somedata.mrc"
# map_size = (512, 512, 200)
# map_level = 1.245
# ijk_to_xyz_transform = <3x4 matrix>
#
# Region attributes are each written in a separate group named to match
# the attribute name with a type int, float, string appened to the name.
# An array named "attributes" contains names of these group nodes.
#
# attributes = <array of node names>, e.g. ["curvature float", ...]
#
# /curvature float
#   attribute_name = "curvature"
#   ids = <array of region indices, length M>
#   values = <array of values, length M>
#
# skeleton = 'string encoding chimera marker file'
#
# The file is saved with the Python PyTables modules which includes
# additional attributes "VERSION", "CLASS", "TITLE", "PYTABLES_FORMAT_VERSION".
#
# Tests with alternate data storage with every region being a separate HDF
# node and every contact being a separate HDF node gave extremely slow
# read/write speed.
#
def write_segmentation(seg, path = None):

  if path is None:
    show_save_dialog(seg)
    return

  import tables
  h5file = tables.openFile(path, mode = 'w')

  try:

    root = h5file.root
    a = root._v_attrs
    a.format = 'segger'
    a.format_version = 2
    a.name = seg.name

    m = seg.mask
    atom = tables.Atom.from_dtype(m.dtype)
    filters = tables.Filters(complevel=5, complib='zlib')
    ma = h5file.createCArray(root, 'mask', atom, m.shape, filters = filters)
    ma[:] = m
    
    print " - updating region colors..."
    seg.region_colors ()

    from numpy import array, int32, float32
    rlist = seg.id_to_region.values()
    rlist.sort(lambda r1,r2: cmp(r1.rid,r2.rid))

    rids = array([r.rid for r in rlist], int32)
    h5file.createArray(root, 'region_ids', rids)

    rcolors = array([r.color for r in rlist], float32)
    h5file.createArray(root, 'region_colors', rcolors)

    refpts = array([r.max_point for r in rlist], float32)
    h5file.createArray(root, 'ref_points', refpts)

    slev = array([r.smoothing_level for r in rlist], float32)
    h5file.createArray(root, 'smoothing_levels', slev)

    pids = array([(r.preg.rid if r.preg else 0) for r in rlist], int32)
    h5file.createArray(root, 'parent_ids', pids)

    map = seg.volume_data()
    if map:
      d = map.data
      a.map_path = d.path
      print " - map path: ", d.path
      a.map_size = array(d.size, int32)

    if not seg.map_level is None:
      a.map_level = seg.map_level

    t = seg.point_transform()
    if t:
      from numpy import array, float32
      a.ijk_to_xyz_transform = array(t, float32)

    write_attributes(h5file, seg)

    if seg.adj_graph:
      write_skeleton(h5file, seg.adj_graph)

  finally:

    h5file.close()

  seg.path = path

# -----------------------------------------------------------------------------
#
def show_save_dialog(seg, saved_cb = None):

  def save ( okay, dialog, saved_cb = saved_cb ):
    if okay:
      paths = dialog.getPaths ( )
      if paths:
        write_segmentation ( seg, paths[0] )
        if saved_cb:
          saved_cb(seg)

  if hasattr(seg, 'path'):
    import os.path
    idir, ifile = os.path.split(seg.path)
  else:
    idir = None
    ifile = seg.name

  from OpenSave import SaveModeless
  SaveModeless ( title = 'Save Segmentation %s' % seg.name,
                 filters = [('Segmentation', '*.seg', '.seg')],
                 initialdir = idir, initialfile = ifile, command = save )

# -----------------------------------------------------------------------------
#
def write_attributes(h5file, seg):

  aa = {}
  for r in seg.all_regions():
    for a,v in r.attributes().items():
      ta = (a, attribute_value_type(v))
      if ta in aa:
        aa[ta].append((r.rid, v))
      else:
        aa[ta] = [(r.rid, v)]

  if len(aa) == 0:
    return              # HDF5 doesn't handle 0 length arrays.

  gnames = []
  from numpy import array, uint32, int32, float64
  for (a,t), vals in aa.items():
    gname = a.replace('/','_') + ' ' + t
    gnames.append(gname)
    g = h5file.createGroup("/", gname, 'region attribute')
    g._v_attrs.attribute_name = a
    rid = array([i for i,v in vals], uint32)
    h5file.createArray(g, 'ids', rid, 'region id numbers')
    if t == 'int':
      va = array([v for i,v in vals], int32)
    elif t == 'float':
      va = array([v for i,v in vals], float64)
    elif t == 'string':
      va = [v for i,v in vals]
    elif t == 'image':
      va = [image_to_string(v) for i,v in vals]
    h5file.createArray(g, 'values', va, 'attribute values')
    if t == 'image':
      g._v_attrs.value_type = 'PNG image'

  h5file.createArray(h5file.root, 'attributes', gnames)


# -----------------------------------------------------------------------------
#
def read_attributes(h5file, seg):

  r = h5file.root
  if not hasattr(r, 'attributes'):
    return

  id2r = seg.id_to_region
  for gname in r.attributes:
    g = getattr(r, gname)
    a = g._v_attrs.attribute_name
    ids = g.ids
    values = g.values
    img = (hasattr(g._v_attrs, 'value_type') and
           g._v_attrs.value_type == 'PNG image')
    for id,v in zip(ids,values):
      if id in id2r:
        if img:
          v = string_to_image(v)
        id2r[id].set_attribute(a, v)

# -----------------------------------------------------------------------------
#
import numpy
int_types = (int, numpy.int8, numpy.int16, numpy.int32, numpy.int64,
             numpy.uint8, numpy.uint16, numpy.uint32, numpy.uint64)
float_types = (float, numpy.float32, numpy.float64)
def attribute_value_type(v):

  from PIL.Image import Image
  if isinstance(v, int_types):
    return 'int'
  elif isinstance(v, float_types):
    return 'float'
  elif isinstance(v, basestring):
    return 'string'
  elif isinstance(v, Image):
    return 'image'

  raise TypeError, "Can't save value type %s" % str(type(v))

# -----------------------------------------------------------------------------
#
def image_to_string(image):

  from cStringIO import StringIO
  s = StringIO()
  image.save(s, 'PNG')
  return s.getvalue()

# -----------------------------------------------------------------------------
#
def string_to_image(string):

  from cStringIO import StringIO
  f = StringIO(string)
  from PIL import Image
  i = Image.open(f)
  return i

# -----------------------------------------------------------------------------
#
def write_skeleton(h5file, mset):

  import StringIO
  s = StringIO.StringIO()
  mset.save_as_xml(s)
  import numpy
  skel = numpy.char.array(str(s.getvalue()), itemsize = 1)
  s.close()
  h5file.createArray(h5file.root, 'skeleton', skel)

# -----------------------------------------------------------------------------
#
def read_segmentation(path, open = True, task = None):

  import tables
  f = tables.openFile(path)

  try :

    r = f.root
    a = r._v_attrs
    for n in ('format', 'format_version', 'name'):
      if not n in a:
        raise ValueError, 'Segmentation file does not have "%s" attribute' % n
    if a.format != 'segger':
      raise ValueError, 'Segmentation file format is not "segger"'
    if a.format_version != 2:
      raise ValueError, 'Segmentation file format is not 2'

    import os.path
    fname = os.path.basename(path)

    from regions import Segmentation
    s = Segmentation(fname, open = open)

    if 'map_path' in a:
      s.map_path = a.map_path
      print " - map path: " + s.map_path
    if 'map_level' in a:
      s.map_level = a.map_level
    if 'map_name' in a:
      s.map_name = a.map_name
      print " - map name: " + s.map_name

    if 'ijk_to_xyz_transform' in a:
      s.ijk_to_xyz_transform = a.ijk_to_xyz_transform

    s.mask = r.mask.read()
    #print "mask:"pl
    #print s.mask
    rids = r.region_ids.read()
    #print "rids:"
    #print rids
    rcolors = r.region_colors.read()
    refpts = r.ref_points.read()
    slevels = r.smoothing_levels.read() if hasattr(r, 'smoothing_levels') else None
    pids = r.parent_ids.read()
    #print "pids:"
    #print pids

    create_regions(s, rids, rcolors, refpts, slevels, pids, task)

    print " - created regions"

    read_attributes(f, s)

    read_skeleton(f, s)

    read_patches (f, s)

  finally:

    f.close()

  s.path = path

  print " - done reading seg file: " + path
  return s


# -----------------------------------------------------------------------------
#
def create_regions(s, rids, rcolors, refpts, slevels, pids, task):

  if task:
    task.updateStatus('Making ID table')
  id_to_index = dict([(id,i) for i,id in enumerate(rids)])

  if task:
    task.updateStatus('Collecting child region IDs')
  id_to_child_ids = {}
  n = len(rids)
  for i in range(n):
    pid = pids[i]
    if pid > 0:
      if pid in id_to_child_ids:
        id_to_child_ids[pid].append(rids[i])
      else:
        id_to_child_ids[pid] = [rids[i]]

  if task:
    task.updateStatus('Ordering IDs')
  from regions import Region
  ids = depth_order(rids, id_to_child_ids, set())
  rlist = []
  for c,rid in enumerate(ids):
    if rid in id_to_child_ids:
      children = [s.id_to_region[cid] for cid in id_to_child_ids[rid]]
    else:
      children = []
    i = id_to_index[rid]
    r = Region(s, rid, refpts[i], children)
    # TODO: Get wrappy error setting surface piece color to numpy array.
    r.color = tuple(rcolors[i])
    if not slevels is None:
      r.smoothing_level = slevels[i]
    rlist.append(r)
    if task and c % 1000 == 0:
      task.updateStatus('Created %d of %d regions' % (c,n))

  if not slevels is None:
    s.smoothing_level = max(slevels)

  return rlist

# -----------------------------------------------------------------------------
#
def depth_order(rids, id_to_child_ids, used):

  idlist = []
  for rid in rids:
    if not rid in used:
      used.add(rid)
      if rid in id_to_child_ids:
        cids = id_to_child_ids[rid]
        idlist.extend(depth_order(cids, id_to_child_ids, used))
      idlist.append(rid)
  return idlist

# -----------------------------------------------------------------------------
#
def show_open_dialog(dir, callback):

  def open ( okay, dialog ):
    if okay:
      paths_types = dialog.getPathsAndTypes ( )
      if paths_types:
        callback ( paths_types )

  from OpenSave import OpenModeless
  OpenModeless ( title = 'Open Segmentation',
                 initialdir = dir,
                 filters = [('Segmentation', ['*.seg']),
                            ('Old regions file', ['*_regions'])],
                 defaultFilter = 'Segmentation',
                 command = open )

# -----------------------------------------------------------------------------
#
def read_skeleton(f, s):

  a = f.root
  if not 'skeleton' in a :
    return

  sks = a.skeleton.read().tostring()
  import StringIO
  xml = StringIO.StringIO(sks)
  from VolumePath import markerset
  marker_sets = markerset.load_marker_set_xml(xml, model_id = s.id)
  skel = marker_sets[0]
  skel.show_model(True)

  # Map markers to regions
  id2r = dict([(r.rid, r) for r in s.all_regions()])
  for m in skel.markers():
    rid = int(m.extra_attributes['region_id'])
    m.region = id2r.get(rid, None)
    if m.region is None:
      print 'missing skeleton region', rid

  s.adj_graph = skel


# -----------------------------------------------------------------------------
#
def read_patches(f, s):

    a = f.root
    if not 'patches' in a :
        return

    print " - reading patches:"

    patches = list( a.patches )
    print patches

    import chimera
    import Mesh
    reload ( Mesh )
    mesh = None

    for rg in patches:

        rgPath = rg._v_pathname
        print " - path: ", rgPath
        rid = rgPath [ len("/patches/"): ]
        print "  - patch for region ", rid, " - type: ", rg._v_attrs["type"]

        if not 'verts' in rg or not 'tris' in rg :
            print "  - tris or verts not found"
            continue

        #print " - region ids:"
        #print s.id_to_region

        try :
            reg = s.id_to_region[int(rid)]
        except :
            print " - did not find region for id"
            continue


        verts = rg.verts.read()
        tris = rg.tris.read()
        print "   - %d verts, %d tris" % (len(verts), len(tris))
        #print verts
        #print tris

        if mesh == None :
            mesh = Mesh.MeshFromVertsTris (verts, tris, color=reg.color, m=mesh)
            mesh.name = "patches"
            chimera.openModels.add([mesh])
        else :
            mesh = Mesh.MeshFromVertsTris (verts, tris, color=reg.color, m=mesh)
