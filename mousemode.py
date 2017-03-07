# -----------------------------------------------------------------------------
# Add mouse mode to click on density map and group watershed regions connected
# to the clicked on region by density above the display contour level.
#

class Group_Connected_Mouse_Mode:

  def __init__(self):

    self.bound_button = None
    self.mode_name = 'group connected'
    self.region = None
    self.level = None
    self.last_level = None      # Last used grouping level.
    self.last_y = 0
    self.dragging = False
    self.connected = set()
    self.ungrouped_color = (0.8, 0.8, 0.8, 1.0)
    self.background_show_all = True
    
    from chimera import mousemodes
    callbacks = (self.mouse_down_cb, self.mouse_drag_cb, self.mouse_up_cb)
    icon = None
    mousemodes.addFunction(self.mode_name, callbacks, icon)

  def bind_mouse_button(self, button, modifiers):

    self.unbind_mouse_button()
    from chimera import mousemodes
    mousemodes.setButtonFunction(button, modifiers, self.mode_name)
    self.bound_button = (button, modifiers)
    
  def unbind_mouse_button(self):

    if self.bound_button:
      button, modifiers = self.bound_button
      from chimera import mousemodes
      def_mode = mousemodes.getDefault(button, modifiers)
      if def_mode:
        mousemodes.setButtonFunction(button, modifiers, def_mode)
      self.bound_button = None
  
  def mouse_down_cb(self, viewer, event):

    s = self.segmentation()
    if s is None:
      return

    if self.region:
      # Mouse up when segmentation first created is lost on Mac.
      self.mouse_up_cb(viewer, event)

    r = clicked_mask_region(s, event.x, event.y)
    if r is None:
      # Clicked on background.
      all = self.background_show_all
      shift_mask = 1
      if event.state & shift_mask:
        if all:
          unhide_regions(s.all_regions())
        else:
          unhide_regions([r for r in s.all_regions() if not r.in_group()])
      else:
        if all:
          show_groups(s.regions)
        else:
          g = set([r for r in s.regions if r.is_group()])
          from regions import childless_regions, boundary_groups
          g.update(boundary_groups(childless_regions(g), s.region_contacts()))
          show_groups(g)
      self.background_show_all = not all
      return
    else:
      self.background_show_all = False

#    print 'clicked on region', r.rid

    t = r.top_parent()
    from regions import random_color
    self.color = random_color(self.ungrouped_color) if t is r else t.color
    self.level = self.density_level(r)
    lev = self.level if t is r else None
    self.region = r
    self.last_y = event.y
    self.dragging = False
    hide_regions(s.all_regions())
    self.connected = show_connected_regions(r, lev, self.color,
                                            set(), self.ungrouped_color)
    if self.level is None:
        print 'mouse down: no density map, needed for contact densities'
    
  def mouse_drag_cb(self, viewer, event):

    r = self.region
    if r is None:
      return

    dy = event.y - self.last_y
    if abs(dy) < 3 and not self.dragging:
      return          # Ignore small initial motion.
    
    lev = self.level
    if lev is None:
        return

    self.last_y = event.y
    self.dragging = True

    efactor_pixels = 300.0

    shift_mask = 1
    if event.state & shift_mask:
      efactor_pixels *= 10

    import math
    lev *= math.exp(dy/efactor_pixels)

    self.last_y = event.y
    self.level = lev
    self.connected = show_connected_regions(r, lev, self.color, self.connected,
                                            self.ungrouped_color)
    
  def mouse_up_cb(self, viewer, event):

    print 'got mouse up'
    r = self.region
    c = self.connected
    if r and len(c) >= 1:
      s = r.segmentation
      rg = s.join_regions(c, color = self.color)
      self.last_level = rg.connect_level = self.level
      from regions import boundary_groups
      unhide_regions(boundary_groups(c, s.region_contacts()))

    self.region = None
    self.level = None
    self.color = None
    self.last_y = None
    self.dragging = False
    self.connected = set()

  def segmentation(self):

    import segment_dialog as sd
    s = sd.current_segmentation()
    if s is None and sd.segmentation_map():
      d = sd.volume_segmentation_dialog()
      s = d.Segment(group = False)
      for r in s.regions:
        r.set_color(self.ungrouped_color)
      v = s.volume_data()
      if v:
        v.display = False
    if s is None:
      print 'mouse down: no current segmentation'
      return None
    return s

  def density_level(self, region):

    p = region.top_parent()
    if hasattr(p, 'connect_level'):
      return p.connect_level

    if p.has_children():
      rcons = p.segmentation.region_contacts()
      d = contact_density(p.childless_regions(), rcons)
      if not d is None:
        return d

    if not self.last_level is None:
      return self.last_level

    # Get density level to connect at.
    m = region.segmentation.volume_data()
    if m and m.surface_levels:
      return min(m.surface_levels)

    return None

    
# ---------------------------------------------------------------------------
#
def clicked_mask_region(segmentation, pointer_x, pointer_y):

  s = segmentation
  if not s.display:
    return None

  if s.mask is None:
    print 'mouse down: no current segmentation mask'
    return None

  mk,mj,mi = s.mask.shape
  box = ((0,0,0),(mi,mj,mk))
  from Matrix import xform_matrix, multiply_matrices
  ijk_to_eye = multiply_matrices(xform_matrix(s.openState.xform),
                                 s.point_transform())
  from VolumeViewer.slice import box_intercepts, array_slice_values
  ijk_in, ijk_out = box_intercepts(pointer_x, pointer_y, ijk_to_eye, box, s)
  if ijk_in is None or ijk_out is None:
    print 'mouse down: no intercept with mask'
    return None

  rnums = array_slice_values(s.mask, ijk_in, ijk_out, method = 'nearest')[:,1]
  r = first_shown_region(s, rnums)
  if r is None:
    r = clicked_volume_region(segmentation, pointer_x, pointer_y)
    if r is None:
      print 'mouse down: no intercept with region'
      return None

  return r
    
# ---------------------------------------------------------------------------
#
def clicked_volume_region(segmentation, pointer_x, pointer_y):

  r = None
  s = segmentation
  m = s.mask
  v = s.volume_data()
  if v and v.shown() and v.single_plane() and not m is None:
    box = v.region[:2]
    from Matrix import xform_matrix, multiply_matrices
    ijk_to_eye = multiply_matrices(xform_matrix(s.openState.xform),
                                   s.point_transform())
    from VolumeViewer.slice import box_intercepts
    ijk_in, ijk_out = box_intercepts(pointer_x, pointer_y, ijk_to_eye, box, v)
    if not ijk_in is None:
      i,j,k = [int(round(h)) for h in ijk_in]
      ksz,jsz,isz = m.shape
      if i >= 0 and i < isz and j >= 0 and j < jsz and k >= 0 and k < ksz:
        rid = s.mask[k,j,i]
        r = s.id_to_region.get(rid, None)
  return r
        
# ---------------------------------------------------------------------------
#
def first_shown_region(segmentation, rnums):

  id2r = segmentation.id_to_region
  for rnum in rnums:
    if rnum > 0:
      r = id2r[rnum]
      if r.visible():
        return r
  return None
    
# ---------------------------------------------------------------------------
#
def clicked_density_region(segmentation, pointer_x, pointer_y):

  s = segmentation

  if s.mask is None:
    print 'mouse down: no current segmentation mask'
    return None
  
  v = s.volume_data()
  if v is None:
    print 'mouse down: no current segmentation map'
    return None

  v_xyz = density_maximum(v, pointer_x, pointer_y, s.map_level)
  if v_xyz is None:
    return

  i,j,k = [int(round(i)) for i in v.data.xyz_to_ijk(v_xyz)]
  rnum = s.mask[k,j,i]
  if rnum == 0:
    print 'mouse down: outside contour'
    return None

  r = s.id_to_region[rnum]
  return r
    
# ---------------------------------------------------------------------------
#
def density_maximum(v, pointer_x, pointer_y, map_level):

  from VolumeViewer.slice import volume_segment
  xyz_in, xyz_out = volume_segment(v, pointer_x, pointer_y)

  if xyz_in is None or xyz_out is None:
    print 'mouse down: no intersection with volume'
    return None

  t_max, v_max = first_maximum_above_threshold(v, xyz_in, xyz_out,
                                               map_level)
  if t_max is None:
    print 'mouse down: no maximum under mouse'
    return None

  from VolumePath.gui import linear_combination
  v_xyz = linear_combination(1-t_max, xyz_in, t_max, xyz_out)
  return v_xyz
    
# ---------------------------------------------------------------------------
#
def first_maximum_above_threshold(v, xyz_in, xyz_out, thresh):

    from VolumePath.gui import slice_data_values
    trace = slice_data_values(v, xyz_in, xyz_out)

    n = len(trace)
    for k in range(n):
      t, v = trace[k]
      if v >= thresh:
        if ((k-1 < 0 or trace[k-1][1] < v) and
            (k+1 >= n or trace[k+1][1] <= v)):
          return t, v
    return None, None
    
# ---------------------------------------------------------------------------
#
def show_connected_regions(region, level, color, rset, ungrouped_color):

  s = region.segmentation

  # Remove grouping
  top = region.top_parent()
  rsibling = top.childless_regions()
  parents = [p for p in top.all_regions() if p.has_children()]
  if parents:
    s.remove_regions(parents, update_surfaces = False)
  if unhide_regions(rsibling) == 0:
    for r in rsibling:
      r.make_surface()

  # Find connected regions
  rcons = s.region_contacts()
  if level is None:
    cset = set(rsibling)
  else:
    cset = set([region])
    bndry = set([region])
    dmax = None
    while bndry:
      r = bndry.pop()
      if r in rcons:
        for cr, c in rcons[r].items():
          if cr.preg is None and not cr in cset:
            if c.maximum_density >= level:
              cset.add(cr)
              bndry.add(cr)
            if dmax is None or c.maximum_density > dmax:
              dmax = c.maximum_density

  # Show only connected regions and boundary.
  from regions import boundary_regions
  dset = cset | boundary_regions(cset, rcons)   # display
  hset = rset | boundary_regions(rset, rcons)   # extra
  for r in dset | hset:
    sp = r.surface_piece
    if sp:
      sp.display = (r in dset)

  if (len(cset) == 1 and not level is None and not dmax is None and
      level > 1.2 * dmax):
    cset.clear()

  # Update colors.
  for r in cset:
    r.set_color(color)
  for r in rset - cset:
    r.set_color(ungrouped_color)        # Uncolor regions no longer in group
  
#  print 'connected %d regions to region %d at level %.3g' % (len(cset), region.rid, level)
  
  return cset

# -----------------------------------------------------------------------------
#
def contact_density(regions, rcons, scale = 1.000001):

  demax = None           # external contacts
  dimin = None           # internal contacts
  rset = set(regions)
  for r in rset:
    if r in rcons:
      for rc, c in rcons[r].items():
        if not c.maximum_density is None:
          if rc in rset:
            if dimin is None or c.maximum_density < dimin:
              dimin = c.maximum_density
          elif not rc.in_group():
            if demax is None or c.maximum_density > demax:
              demax = c.maximum_density
  if not demax is None:
    d = demax * scale
  elif not dimin is None:
    d = dimin / scale
  else:
    d = None
  return d

# -----------------------------------------------------------------------------
#
def hide_regions(regions):

  count = 0
  for r in regions:
    sp = r.surface_piece
    if sp:
      sp.display = False
      count += 1
  return count

# -----------------------------------------------------------------------------
#
def unhide_regions(regions):

  count = 0
  for r in regions:
    sp = r.surface_piece
    if sp:
      sp.display = True
      count += 1
  return count

# -----------------------------------------------------------------------------
#
def show_groups(regions):

  for r in regions:
    if r.has_surface():
      r.show_surface()
    else:
      r.make_surface()
    hide_regions(r.all_children())
