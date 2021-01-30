# ------------------------------------------------------------------------------
# Produce tiled grids of cross-sections of density of segmented regions along
# 3 orthogonal axes of straightened density.
#
def make_orthoslice_images(
    rlist,                       # Segmentation regions.
    volume,                      # Volume model.
    trace_spacing = None,        # Physical units.
    trace_tip_length = None,     # Physical units.
    unbend_size = 1.5,           # Factor multipied by region diameter.
    unbend_yaxis = (0,0,1),
    unbend_grid_spacing = 1,     # Factor multiplied by minimum voxel size.
    slice_spacing = None,        # 3-tuple, physical units.
    xy_trim = 0.3,               # Factor multipled by unbent grid width/height.
    panel_aspect = 0.5,          # Minimum aspect ratio for tiled layout.
    image_spacing = 20,          # Pixels between 3 sets of slices.
    show_image = True,
    task = None,
    ):

    from math import ceil
    from Measure.spine import trace_spine, measure_diameter
    from VolumeFilter.unbend import atom_path, unbend_volume
    from VolumeFilter.tile import tile_planes
    from chimera import openModels

    ubgs = unbend_grid_spacing * min(volume.data.step)

    axes = ('x','y','z')
    pstep = (10,10,10) if slice_spacing is None else [int(ceil(s/ubgs))
                                                      for s in slice_spacing]
    orders = ('ulh','urv','ulhr')

    for ri,r in enumerate(rlist):

        if task:
            task.updateStatus('region %d (%d of %d)' % (r.rid,ri+1,len(rlist)))

        # Trace center-line.
        mset = trace_spine(r, trace_spacing, trace_tip_length) 

        # Unbend volume.
        p = atom_path([m.atom for m in mset.markers()])
        dmax, dmin = measure_diameter(r, mset)
        if dmax is None:
            print 'Region %d has no diameter' % r.rid
            mset.close()
            continue
        xsize = ysize = unbend_size*dmax
        ubv = unbend_volume(volume, p, unbend_yaxis, xsize, ysize, ubgs,
                            open = False)
#        ubv.set_representation('solid')
#        ubv.set_parameters(show_outline_box = True)

        # Create tiled cross-section volumes along 3 axes.
        etrim = (xy_trim*ubv.data.size[0],xy_trim*ubv.data.size[1],0)
        etrim = [int(ceil(t)) for t in etrim]
        rowcol = rows_and_columns(ubv.data.size, pstep, etrim, panel_aspect)
        tparams = zip(axes,pstep,etrim,rowcol,orders)
        tv = [tile_planes(ubv, axis, step, trim, rows, cols, order, open=False)
              for axis,step,trim,(rows,cols),order in tparams]
        if [v for v in tv if v is None]:
            print 'Region %d has no sections for some axes' % r.rid
            mset.close()
            openModels.close([v for v in tv if v] + [ubv])
            continue

        # Make images for each set of slices and combine them side-by-side.
        images = (volume_image(tv[0]),
                  volume_image(tv[1], xflip=True),
                  volume_image(tv[2]))
        image = montage_image(images, pad = image_spacing)
        image.title = 'Region %d' % r.rid
        if show_image:
            from Segger.imageviewer import showImage
            showImage(image, title = image.title, fit_to_window = True)

        # Set image region attribute
        r.set_attribute('slices', image)

        # Close center line.
        mset.close()

# ------------------------------------------------------------------------------
# Figure out optimal row/column tiling for each set of sections so that the
# composite 3 panels has a convenient aspect ratio around 1.5.
# yz slices
#
def rows_and_columns(vsize, step, trim, min_aspect):

    from math import ceil, floor
    tc = [max(1,(s-t+st-1)/st - (t+st-1)/st) for s,st,t in zip(vsize,step,trim)]
    tc0 = tc[0]
    c = 1
    while c <= tc0:
        r = int(ceil((tc0-1+c)/c))
        aspect = float(c*vsize[1])/(r*vsize[2])
        if r <= 1 or aspect >= min_aspect:
            break
        while int(ceil((tc0-1+c)/c)) == r:
            c += 1
    rz = max(1,int(floor(float(r*vsize[2])/vsize[1])))
    rowcol = ((None,c),(c,None),(rz,None))
    return rowcol

# ------------------------------------------------------------------------------
# Look down x, y, z axes.  For x and y use z as vertical.
#
def position_orthoviews(tv):

    # Put all 3 planes in xy plane.
    vx, vy, vz = tv
    from chimera import Xform
    vx.openState.localXform(Xform.xRotation(-90))
    vx.openState.localXform(Xform.zRotation(-90))
    vy.openState.localXform(Xform.xRotation(-90))
    vy.openState.localXform(Xform.zRotation(180))
    # Horz/vert axes vx:(y,z), vy:(x,z), vz:(x,y)

    sizes = []
    for v in tv:
        (x0,y0,z0), (x1,y1,z1) = v.xyz_bounds(step = 1, subregion = 'all')
        sizes.append((x1-x0,y1-y0,z1-z0))

    # Align lower left corners
    vy.openState.localXform(Xform.translation(-sizes[1][0]))

    # Center vertically and spread out horizontally.
    vx.openState.localXform(Xform.translation(0,0,-0.5*sizes[0][2]))
    offset = 1.05 * sizes[0][1]
    vy.openState.localXform(Xform.translation(offset,0,-0.5*sizes[1][2]))
    offset += .05 * sizes[0][1] + sizes[1][0]
    vz.openState.localXform(Xform.translation(offset,-0.5*sizes[2][1],0))

    # Hide all other models.
    from chimera import openModels
    mset = set(tv + [v.solid_model() for v in tv])
    for m in openModels.list():
        if not m in mset:
            m.display = False

    # Standard orientation, fit in window.
    from Midas import reset
    reset()
    from chimera import viewer
    viewer.viewAll()
    viewer.scaleFactor *= 0.95
    viewer.depthCue = False

    # Capture image
#        width,height = viewer.windowSize
#        image = viewer.pilImages(width, height, supersample = 3)[0]

# ------------------------------------------------------------------------------
# Creates a PIL image from a volume plane with 1 pixel exactly matching one
# voxel, and using the solid rendering transfer function.
#
# The volume must be rendering as single plane in solid rendering style.
#
def volume_image(v, xflip = False, yflip = False):

    if v.solid is None:
        v.update_solid(v.rendering_options)
    colors = v.solid.color_values()
    if colors.shape[0] == 1:
        c = colors[0,:,:,:]
    elif colors.shape[1] == 1:
        c = colors[:,0,:,:]
    elif colors.shape[2] == 1:
        c = colors[:,:,0,:]
    else:
        raise ValueError('Color array not a plane %s' % str(colors.shape))
    size = (c.shape[1], c.shape[0])
    mode = {1:'L', 2:'LA', 3:'BGR', 4:'BGRA'}[c.shape[2]]
    if mode.endswith('A'):
        # Don't include alpha channel in image.
        mode = mode[:-1]
        c = c[:,:,0] if c.shape[2] == 2 else c[:,:,:-1]
    if xflip:
        c = c[:,::-1]
    if not yflip:
        c = c[::-1,:]
    from PIL import Image
    im = Image.fromarray(c, mode)
    return im

# ------------------------------------------------------------------------------
#
def montage_image(images, pad = 0):

    w = sum([im.size[0] for im in images]) + pad * (len(images)-1)
    h = max([im.size[1] for im in images])
    from PIL import Image
    mi = Image.new(images[0].mode, (w,h))
    x = y = 0
    for im in images:
        mi.paste(im, (x,y))
        x += im.size[0] + pad
    return mi
    
# ------------------------------------------------------------------------------
#
def test():

    from Segger import SelectedRegions
    rlist = SelectedRegions()

    from VolumeViewer import active_volume
    volume = active_volume()

    make_orthoslice_images(rlist, volume)

#test()
