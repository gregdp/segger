# ------------------------------------------------------------------------------
# segment copygroups #2 #3
# segment unbin #1 #5
#

# -----------------------------------------------------------------------------
#
def segment_command(cmdname, args):

    from Commands import perform_operation, string_arg, volume_arg
    from Commands import bool_arg, float_arg, float3_arg, int_arg, int3_arg

    ops = {'copygroups': (copy_groups,
                          (('from_seg', segmentation_arg),
                           ('to_seg', segmentation_arg)),
                          (),
                          ()),
           'unbin': (unbin_mask,
                     (('segmentation', segmentation_arg),
                      ('volume', volume_arg)),
                     (),
                     ()),
           'directioncolor': (color_by_direction,
                              (('segmentation', segmentation_arg),),
                              (),
                              (('pattern', string_arg),
                               ('spherekey', bool_arg))),
           'sliceimage': (slice_view,
                          (('rlist', regions_arg),
                           ('volume', volume_arg)),
                          (),
                          (('traceSpacing', float_arg),
                           ('traceTipLength', float_arg),
                           ('unbendSize', float_arg),
                           ('unbendYAxis', float3_arg),
                           ('unbendGridSpacing', float_arg),
                           ('sliceSpacing', float3_arg),
                           ('xyTrim', float_arg),
                           ('panelAspect', float_arg),
                           ('imageSpacing', int_arg),
                           ('showImage', bool_arg))),
           'exportmask': (export_mask,
                          (('segmentation', segmentation_arg),),
                          (),
                          (('savePath', string_arg),
                           ('format', string_arg),
                           ('sequentialIds', bool_arg),
                           ('binSize', int3_arg))),
           }

    perform_operation(cmdname, args, ops)


# -----------------------------------------------------------------------------
#
def copy_groups(from_seg, to_seg):

  # Find transform from to_seg mask indices to from_seg mask indices.
  from Matrix import multiply_matrices, invert_matrix, xform_matrix
  tf = multiply_matrices(invert_matrix(from_seg.point_transform()),
                         invert_matrix(xform_matrix(from_seg.openState.xform)),
                         xform_matrix(to_seg.openState.xform),
                         to_seg.point_transform())

  # Find from_seg regions containing to_seg region maximum points.
  fmask = from_seg.mask
  rlist = list(to_seg.regions)
  ijk = [r.max_point for r in rlist]
  from _interpolate import interpolate_volume_data
  rnums, outside = interpolate_volume_data(ijk, tf, fmask, 'nearest')

  # Find to_seg regions that have max_point in same from_seg region.
  fr2tr = {}
  id2r = from_seg.id_to_region
  for tr, frnum in zip(rlist, rnums):
    if frnum > 0:
      fr = id2r[frnum].top_parent()
      if fr in fr2tr:
        fr2tr[fr].append(tr)
      else:
        fr2tr[fr] = [tr]
  groups = [g for g in fr2tr.values() if len(g) > 1]

  # Form groups.
  for g in groups:
    r = to_seg.join_regions(g)
    for gr in g:
        gr.set_color(r.color)

  print 'Made %d groups for %s matching %s' % (len(groups), to_seg.name, from_seg.name)
  

# -----------------------------------------------------------------------------
#
def unbin_mask(segmentation, volume):

    from regions import bin_size, Segmentation

    ssize = segmentation.grid_size()
    vsize = volume.data.size
    bsize = bin_size(ssize, vsize)
    if bsize is None:
        from Commands import CommandError
        raise CommandError, 'Map size %d,%d,%d is not compatible with segmentation size %d,%d,%d' % (tuple(ssize) + tuple(vsize))

    name = segmentation.name + ' unbin'
    seg = Segmentation(name, volume)

    # Copy mask
    b2, b1, b0 = bsize
    s2, s1, s0 = [s*b for s, b in zip(ssize, bsize)]
    for o0 in range(b0):
        for o1 in range(b1):
            for o2 in range(b2):
                seg.mask[o0:s0:b0,o1:s1:b1,o2:s2:b2] = segmentation.mask

    # Copy regions
    copy_regions(segmentation.all_regions(), bsize, seg)
    
    return s

# -----------------------------------------------------------------------------
#
def copy_regions(regions, bin_size, seg):

    from regions import Region
    rlist = []
    for r in regions:
        if r.rid in seg.id_to_region:
            rc = seg.id_to_region[r.rid]
        else:
            cc = copy_regions(r.children(), bin_size, seg)
            max_point = tuple([b*i for b,i in zip(bin_size, r.max_point)])
            rc = Region(seg, r.rid, max_point, cc)
        rlist.append(rc)
    return rlist


# -----------------------------------------------------------------------------
# Color segmented regions according to the direction of their principle axis.
# This is to visualize how much neighbor bacteria align with each other.
#
def color_by_direction(segmentation, pattern = 'circle',
                       spherekey = False):

    regions = segmentation.regions
    from Measure.spine import points_long_axis
    for region in regions:
        axis = points_long_axis(region.map_points())
        color = direction_color(axis, pattern)
        region.set_color(color)

    if spherekey:
        radius = 0.2 * segmentation.bsphere()[1].radius
        color_sphere(radius, pattern)
    else:
        remove_sphere_key()

# -----------------------------------------------------------------------------
#
def direction_color(axis, pattern = 'circle', opacity = 1):

  import Matrix
  x,y,z = Matrix.normalize_vector(axis)
  if pattern in ('circle', 'circle111'):
      if pattern == 'circle':
          tf = Matrix.vector_rotation_transform((0,0,1),(1,1,1))
          x,y,z = Matrix.apply_matrix(tf, (x,y,z))
      r, g, b = (x*x+2*y*z+1)/2, (y*y+2*x*z+1)/2, (z*z+2*x*y+1)/2
      from math import pow
      r, g, b = [pow(k,0.25) for k in (r,g,b)]
  elif pattern == 'rgb':
      r, g, b = abs(x), abs(y), abs(z)
  elif pattern == 'rgb2':
      r, g, b = x*x, y*y, z*z
  elif pattern == 'cmy':
      r, g, b = y*y+z*z, x*x+z*z, x*x+y*y
      from math import sqrt
      r, g, b = [sqrt(k) for k in (r,g,b)]
  elif pattern == 'cmy2':
      r, g, b = y*y+z*z, x*x+z*z, x*x+y*y
  elif pattern == 'cmz':
      r, g, b = (x*x-2*y*z+1)/2, (y*y-2*x*z+1)/2, (z*z-2*x*y+1)/2
      from math import sqrt
      r, g, b = [sqrt(k) for k in (r,g,b)]
  elif pattern == 'rgb111':
      r, g, b = (x*x+y*z+.5)/1.5, (y*y+x*z+.5)/1.5, (z*z+x*y+.5)/1.5
  else:
      from Commands import CommandError
      raise CommandError('Unknown color pattern %s ("circle", "circle111", "rgb", "rgb2", "cmy", "cmy2", "cmz", "rgb111")'
                         % pattern)

#  r, g, b = [pow(k,0.33) for k in (r,g,b)]
#  r, g, b = [pow(k,0.5) for k in (r,g,b)]
#  n = r+g+b
#  n = max((r,g,b))
#  n = 1
#  r /= n
#  g /= n
#  b /= n

  color = (r,g,b,opacity)
  return color

# -----------------------------------------------------------------------------
# Color sphere so opposite points have the same color, colors vary continuously
# and no two points (except opposite one another) have the same color.
#
# I think this is impossible because it would require mapping the projective
# plane P2 into 3-dimensions which is not possible.
#
# Purpose of this is to color long bacteria according to the orientation of
# their long axes.
#
# Assume sphere surface model already open.
#
def color_sphere(radius = 1, pattern = 'circle'):

  from chimera import openModels as om
  from _surface import SurfaceModel
  cs = [m for m in om.list(modelTypes = [SurfaceModel])
        if m.name == 'Direction colors']
  if cs:
    p = cs[0].surfacePieces[0]
  else:
    from Shape.shapecmd import sphere_shape
    p = sphere_shape(radius = radius, modelName = 'Direction colors')
    
  v,t = p.geometry
  #vmask = [1 if z >= 0 else 0 for x,y,z in v]
  #p.triangleAndEdgeMask = None
  #p.setTriangleMaskFromVertexMask(vmask)
  from numpy import empty, float32
  c = empty((p.vertexCount,4), float32)
  for i, axis in enumerate(v):
      c[i,:] = direction_color(axis, pattern)
  p.vertexColors = c
  p.useLighting = False

#  f = s.addPiece(c[:,:3], t, (0.5,0.5,0.5,1))
#  f.displayStyle = f.Mesh
#  f.setTriangleMaskFromVertexMask(vmask)

  return p
        
# -----------------------------------------------------------------------------
#
def remove_sphere_key():

  from chimera import openModels as om
  from _surface import SurfaceModel
  om.close([m for m in om.list(modelTypes = [SurfaceModel])
            if m.name == 'Direction colors'])
        
# -----------------------------------------------------------------------------
#
def slice_view(**kw):

    # Translate keyword names.
    n2n = {'traceSpacing': 'trace_spacing',
           'traceTipLength': 'trace_tip_length',
           'unbendSize': 'unbend_size',
           'unbendYAxis': 'unbend_yaxis',
           'unbendGridSpacing': 'unbend_grid_spacing',
           'sliceSpacing': 'slice_spacing',
           'xyTrim': 'xy_trim',
           'panelAspect': 'panel_aspect',
           'imageSpacing': 'image_spacing',
           'showImage': 'show_image',
           }
    kw = dict([(n2n[k],v) if k in n2n else (k,v) for k,v in kw.items()])

    from chimera import tasks, CancelOperation
    task = tasks.Task('Slice images', modal = True)
    kw['task'] = task
    import orthoview
    try:
        orthoview.make_orthoslice_images(**kw)
    except CancelOperation:
        pass
    finally:
        task.finished()

# -----------------------------------------------------------------------------
#
def export_mask(segmentation, savePath = None, format = 'mrc',
                binSize = (1,1,1), sequentialIds = True):

    m = segmentation.mask
    origin = segmentation.grid_origin()
    step = segmentation.grid_step()

    # Include only id numbers of top-level region groups.
    from numpy import zeros, empty
    parent = zeros((segmentation.max_region_id+1,), dtype = m.dtype)
    for r in segmentation.all_regions():
        parent[r.rid] = r.top_parent().rid
    if sequentialIds:
        used_ids = list(set(parent))
        used_ids.sort()
        seq_id = zeros((used_ids[-1]+1,), dtype = m.dtype)
        for i,id in enumerate(used_ids):
            seq_id[id] = i
        parent = seq_id[parent]
    array = parent[m]

    # Expand array to match unbinned density map size.
    if tuple(binSize) != (1,1,1):
        step = [float(s)/b for s,b in zip(step,binSize)]
        shape = [a*b for a,b in zip(array.shape,binSize[::-1])]
        aub = empty(shape, array.dtype)
        b2, b1, b0 = binSize
        for o0 in range(b0):
            for o1 in range(b1):
                for o2 in range(b2):
                    aub[o0::b0,o1::b1,o2::b2] = array
        array = aub
                    
    from VolumeData import Array_Grid_Data
    g = Array_Grid_Data(array, origin, step)
    g.name = segmentation.name + ' region ids'

    if savePath is None:
        # Open mask map as a volume.
        import VolumeViewer
        v = VolumeViewer.volume_from_grid_data(g)
        v.openState.xform = segmentation.openState.xform
    else:
        # Write map file.
        from VolumeData import save_grid_data
        save_grid_data(g, savePath, format)

    return g
        
# -----------------------------------------------------------------------------
#
def segmentation_arg(s):

    from chimera.specifier import evalSpec
    sel = evalSpec(s)
    mlist = sel.models()
    from regions import Segmentation
    slist = [s for s in mlist if isinstance(s, Segmentation)]
    if len(slist) == 0:
        raise CommandError, 'No segmentation specified'
    elif len(slist) > 1:
        raise CommandError, 'Multiple segmentations specified'
    return slist[0]
        
# -----------------------------------------------------------------------------
#
def regions_arg(s):

    from chimera.specifier import evalSpec
    sel = evalSpec(s)
    import Surface
    plist = Surface.selected_surface_pieces(sel)
    from Segger import Region
    rlist = [p.region for p in plist
             if hasattr(p,'region') and isinstance(p.region, Region)]
    if len(rlist) == 0:
        raise CommandError, 'No segmentation regions specified'
    return rlist
