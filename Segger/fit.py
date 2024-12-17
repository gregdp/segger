

import chimera
import VolumeViewer
import VolumeData
import numpy
import FitMap
import Matrix
import time
import os
import sys
import shutil


def CopyMolX ( mol, xf ) :

    nmol = chimera.Molecule()
    nmol.name = mol.name

    aMap = dict()
    from random import random as rand
    clr = ( rand(), rand(), rand() )

    for res in mol.residues :
        #nres = nmol.newResidue (res.type, chimera.MolResId(res.id.chainId, res.id.position))
        nres = nmol.newResidue (res.type, chimera.MolResId(res.id.chainId, res.id.position))
        # print "New res: %s %d" % (nres.id.chainId, nres.id.position)
        for at in res.atoms :
            nat = nmol.newAtom (at.name, chimera.Element(at.element.number))
            aMap[at] = nat
            nres.addAtom( nat )
            nat.setCoord ( xf.apply(at.coord()) )
            nat.altLoc = at.altLoc
            nat.occupancy = at.occupancy
            nat.bfactor = at.bfactor
            #if res.isProt or res.isNA :
            #    nat.display = False
            #else :
            #    nat.display = True
            #    nat.radius=1.46
            nat.display = True
            nat.color = chimera.MaterialColor( clr[0], clr[1], clr[2], 1.0 )
            nat.drawMode = nat.EndCap

        nres.isHelix = res.isHelix
        nres.isHet = res.isHet
        nres.isSheet = res.isSheet
        nres.isStrand = res.isStrand
        nres.ribbonDisplay = True
        nres.ribbonDrawMode = 2
        nres.ribbonColor = chimera.MaterialColor( clr[0], clr[1], clr[2], 1.0 );

    for bond in mol.bonds :
        nb = nmol.newBond ( aMap[bond.atoms[0]], aMap[bond.atoms[1]] )
        nb.display = nb.Smart

    return nmol



def CopyAfMol ( mol, xf, bCutoff=55.0 ) :

    nmol = chimera.Molecule()
    nmol.name = mol.name

    aMap = dict()
    from random import random as rand
    clr = ( rand(), rand(), rand() )

    for res in mol.residues :

        #sumb = 0.0
        #for at in res.atoms :
        #    sumb += at.bfactor
        #avgb = sumb / float(len(res.atoms))
        #if avgb < bCutoff :
        #    continue

        if res.atomsMap['CA'][0].bfactor < bCutoff :
            continue

        #nres = nmol.newResidue (res.type, chimera.MolResId(res.id.chainId, res.id.position))
        nres = nmol.newResidue (res.type, chimera.MolResId(res.id.chainId, res.id.position))
        # print "New res: %s %d" % (nres.id.chainId, nres.id.position)
        for at in res.atoms :
            nat = nmol.newAtom (at.name, chimera.Element(at.element.number))
            aMap[at] = nat
            nres.addAtom( nat )
            nat.setCoord ( xf.apply(at.coord()) )
            nat.altLoc = at.altLoc
            nat.occupancy = at.occupancy
            nat.bfactor = at.bfactor
            #if res.isProt or res.isNA :
            #    nat.display = False
            #else :
            #    nat.display = True
            #    nat.radius=1.46
            nat.display = True
            nat.color = chimera.MaterialColor( clr[0], clr[1], clr[2], 1.0 )
            nat.drawMode = nat.EndCap

        nres.isHelix = res.isHelix
        nres.isHet = res.isHet
        nres.isSheet = res.isSheet
        nres.isStrand = res.isStrand
        nres.ribbonDisplay = True
        nres.ribbonDrawMode = 2
        nres.ribbonColor = chimera.MaterialColor( clr[0], clr[1], clr[2], 1.0 );

    for bond in mol.bonds :
        if bond.atoms[0] in aMap and bond.atoms[1] in aMap :
            nb = nmol.newBond ( aMap[bond.atoms[0]], aMap[bond.atoms[1]] )
            nb.display = nb.Smart

    return nmol



def fit_points_g (fdata, threshold = 0.3) :

    mat = fdata.full_matrix()

    import _volume
    points = _volume.high_indices(mat, threshold)
    fpoints = points.astype(numpy.single)
    fpoint_weights = mat[points[:,2],points[:,1],points[:,0]]

    nz = numpy.nonzero( fpoint_weights )[0]
    if len(nz) < len (fpoint_weights) :
        fpoints = numpy.take( fpoints, nz, axis=0 )
        fpoint_weights = numpy.take(fpoint_weights, nz, axis=0)

    from _contour import affine_transform_vertices as transform_vertices
    transform_vertices ( fpoints, fdata.ijk_to_xyz_transform )

    if 0 : print "FitPoints from %s with threshold %.4f, %d nonzero" % (
        fmap.name, threshold, len(nz) )

    return fpoints, fpoint_weights


def resizemap ( dmap, step=6.0, fout=None ) :

    print "step %.2f -> %.2f" % (dmap.data.step[0], step)

    nstep = (step, step, step)

    nn = dmap.data.size

    dim = numpy.array(dmap.data.size) * dmap.data.step
    nn1 = int ( numpy.ceil ( dim[0] / nstep[0] ) )
    nn2 = int ( numpy.ceil ( dim[1] / nstep[1] ) )
    nn3 = int ( numpy.ceil ( dim[2] / nstep[2] ) )
    nnn = [nn1, nn2, nn3]
    print nn, " -> ", nnn

    O = dmap.data.origin

    nmat = numpy.zeros ( (nn1,nn2,nn3), numpy.float32 )
    ndata = VolumeData.Array_Grid_Data ( nmat, O, nstep, dmap.data.cell_angles )

    from VolumeData import grid_indices
    npoints = grid_indices ( (nn1, nn2, nn3), numpy.single)  # i,j,k indices

    from _contour import affine_transform_vertices as transform_vertices
    transform_vertices ( npoints, ndata.ijk_to_xyz_transform )

    # todo - don't interpolate

    dvals = dmap.interpolated_values ( npoints, dmap.openState.xform )
    #dvals = numpy.where ( dvals > threshold, dvals, numpy.zeros_like(dvals) )
    #nze = numpy.nonzero ( dvals )

    nmat = dvals.reshape( (nn3,nn2,nn1) )

    ndata = VolumeData.Array_Grid_Data ( nmat, O, nstep, dmap.data.cell_angles )

    if fout == None :
        try : nv = VolumeViewer.volume.add_data_set ( ndata, None )
        except : nv = VolumeViewer.volume.volume_from_grid_data ( ndata )
        dmap_base = os.path.splitext(dmap.name)[0]
        dmap_path = os.path.splitext (dmap.data.path)[0]
        nv.name = dmap_base + "_step_%.2f" % step
        nv.openState.xform = dmap.openState.xform
        return nv

    else :

        from VolumeData import save_grid_data
        #d = self.grid_data()
        format = save_grid_data(ndata, fout, None, {}, False)
        #print " - saved data"


def SpherePts ( ctr, rad, N ) :

    thetas, phis = [], []
    from math import acos, sin, cos, sqrt, pi
    for k in range ( 1, N+1 ) :
        h = -1.0 + ( 2.0*float(k-1)/float(N-1) )
        phis.append ( acos(h) )
        thetas.append ( 0 if k == 1 or k == N else
                        (thetas[k-2] + 3.6/sqrt(N*(1.0-h**2.0))) % (2*pi) )

    pts = [None] * N
    for i, theta, phi in zip ( range(N), thetas, phis ):
        v = chimera.Vector (sin(phi)*cos(theta), sin(phi)*sin(theta), cos(phi))
        #if numpy.abs ( v.length - 1.0 ) > 1e-3 :
        #    print "x"
        pt = ctr + v * rad
        pts[i] = pt

    return pts





def GetMapGridPts ( dmap, bRad, thr=None, step=10.0, numMax=None, showPts=False ) :

    #dmapr = resizemap ( dmap, step )
    #pts, vals = fit_points_g ( dmapr.data, dmap.surface_levels[0] )
    #print " - %d points above %.2f" % (len(pts), dmap.surface_levels[0])
    #chimera.openModels.close ( [dmapr] )
    #return

    if thr == None :
        #thr = dmap.surface_levels[0]
        #print " - using surface level as threshold: %.2f" % thr
        M = dmap.data.full_matrix()
        thr = numpy.average(M)+numpy.std(M)*2.0
        print " - using threshold at 2-sigma for map points: %.2f" % thr


    points, mapVals = fit_points_g ( dmap.data, thr )
    print " - %d map points above thr %.2f" % ( len(points), thr )

    mpts = []

    #from _contour import affine_transform_vertices
    #points1 = numpy.copy ( points )
    #affine_transform_vertices ( points, Matrix.xform_matrix( dmap.openState.xform ) )

    bound = bRad
    lx,ly,lz = numpy.min ( points, axis=0 ) - (bRad, bRad, bRad)
    hx,hy,hz = numpy.max ( points, axis=0 ) + (bRad, bRad, bRad)

    print " - bounds: %.2f, %.2f, %.2f -> %.2f, %.2f, %.2f" % (lx,ly,lz, hx,hy,hz)

    from numpy import ceil
    nx, ny, nz = int(ceil((hx-lx)/step)), int(ceil((hy-ly)/step)), int(ceil((hz-lz)/step))
    dx, dy, dz = (hx-lx)/float(nx), (hy-ly)/float(ny), (hz-lz)/float(nz)

    print nx, ny, nz

    for i in range ( nx + 1 ) :
        for j in range ( ny + 1 ) :
            for k in range ( nz + 1 ) :
                pt = [ lx + float(i) * dx, ly + float(j) * dy, lz + float(k) * dz ]
                mpts.append ( pt )


    if showPts :
        ptsMol = chimera.Molecule()
        ptsMol.name = "map points step %.2f" % step
        ptsMol.isRealMolecule = False
        chimera.openModels.add ( [ptsMol], noprefs = True )

        for pi, p in enumerate(mpts) :
            res = ptsMol.newResidue('marker', chimera.MolResId("%d"%pi, pi) )
            a = ptsMol.newAtom('', chimera.elements.H)
            res.addAtom(a)
            a.setCoord ( chimera.Point(*p) )  # ( chimera.Point(*xyz) )
            a.radius = 2.0
            a.drawMode = chimera.Atom.Sphere
            a.color = chimera.MaterialColor ( 1.0, 0.0, 0.0, 1.0 )
            a.surfaceCategory = 'markers'

    return mpts


def GetMapPtsGrid ( dmap, thr=None, step=20.0 ) :

    mpts = []
    dim = dmap.data.step * numpy.array(dmap.data.size)
    nz = numpy.floor ( dim[2] / step )
    ny = numpy.floor ( dim[1] / step )
    nx = numpy.floor ( dim[0] / step )
    dz = dim[2] / (nz - 1.0)
    dy = dim[1] / (ny - 1.0)
    dx = dim[0] / (nx - 1.0)
    print " - grid d %f %f %f" % (dx, dy, dz)
    print " - grid n %f %f %f - %d points" % (nx, ny, nz, nx*ny*nz)
    if nx*ny*nz > 1e6 :
        print " - a lot of points %d - try larger step" % (nx,ny,nz)
        stop
    for xi in range ( int(nx) ) :
        for yi in range ( int(ny) ) :
            for zi in range ( int(nz) ) :
                x = dmap.data.origin[0] + xi * dx
                y = dmap.data.origin[1] + yi * dy
                z = dmap.data.origin[2] + zi * dz
                mpts.append ( [x,y,z] )

    print " - grid %d points" % len(mpts)

    if 1 :
        ptsMol = chimera.Molecule()
        ptsMol.name = "map grid points step %.2f" % step
        ptsMol.isRealMolecule = False
        chimera.openModels.add ( [ptsMol], noprefs = True )

        for pi, p in enumerate(mpts) :
            res = ptsMol.newResidue('marker', chimera.MolResId("%d"%pi, pi) )
            a = ptsMol.newAtom('', chimera.elements.H)
            res.addAtom(a)
            a.setCoord ( chimera.Point(*p) )  # ( chimera.Point(*xyz) )
            a.radius = max ( step / 20.0, 2.0 )
            a.drawMode = chimera.Atom.Sphere
            a.color = chimera.MaterialColor ( 1.0, 0.0, 0.0, 1.0 )
            a.surfaceCategory = 'markers'

    return mpts


def GetMapPts ( dmap, thr=None, step=20.0, numMax=None ) :

    #dmapr = resizemap ( dmap, step )
    #pts, vals = fit_points_g ( dmapr.data, dmap.surface_levels[0] )
    #print " - %d points above %.2f" % (len(pts), dmap.surface_levels[0])
    #chimera.openModels.close ( [dmapr] )
    #return

    if thr == None :
        thr = dmap.surface_levels[0]
        print " - using surface level as threshold: %.2f" % thr

    mapPoints, mapVals = fit_points_g ( dmap.data, thr )
    print " - %d map points above thr %.2f" % ( len(mapPoints), thr )

    opts, i = [], 0
    for pt, val in zip ( mapPoints, mapVals ) :
        p = lambda: None
        p.pt, p.val, p.pi, p.ignore = pt, val, i, False
        i += 1
        opts.append ( p )

    #for p in opts[0:10] :
    #    print p.val, p.pt, p.pi

    #print " - making tree..."
    start = time.time()

    from CGLutil.AdaptiveTree import AdaptiveTree
    treeStep = 5.0
    mpTree = AdaptiveTree ( mapPoints.tolist(), opts, treeStep)

    #print " - made tree (step %.1f) in %.1f sec" % ( treeStep, time.time() - start )

    #print " - sorting by value..."
    spts = sorted ( opts, key=lambda x: x.val, reverse = True)
    #for p in spts[0:10] :
    #    print p.val, p.pt, p.pi

    mpts = []

    start = time.time()

    if numMax :
        mpts = spts[0:numMax]

    else :
        #print " - reducing by step %.2f" % step
        step2 = step*step

        for i, p in enumerate(spts) :
            if p.ignore == True :
                continue
            # take this point
            mpts.append ( p.pt )
            # mark other points within _step_ to be ignored
            pta = numpy.array ( p.pt )
            nearPts = mpTree.searchTree ( p.pt, step )
            for npt in nearPts :
                d = (pta - npt.pt)
                D = numpy.sum(d * d)
                if D < step2 :
                    npt.ignore = True
            #if i % 10000 == 0 :
            #    print i,

    print " - reduced with step %.2f to %d - took %.1f sec" % (step, len(mpts), time.time()-start)

    if 1 :
        ptsMol = chimera.Molecule()
        ptsMol.name = "map points step %.2f" % step
        ptsMol.isRealMolecule = False
        chimera.openModels.add ( [ptsMol], noprefs = True )

        for pi, p in enumerate(mpts) :
            res = ptsMol.newResidue('marker', chimera.MolResId("%d"%pi, pi) )
            a = ptsMol.newAtom('', chimera.elements.H)
            res.addAtom(a)
            a.setCoord ( chimera.Point(*p) )  # ( chimera.Point(*xyz) )
            a.radius = 2.0
            a.drawMode = chimera.Atom.Sphere
            a.color = chimera.MaterialColor ( 1.0, 0.0, 0.0, 1.0 )
            a.surfaceCategory = 'markers'

    return mpts



def add_gaussians(grid, xyz, weights, sdev, cutoff_range, transforms = []):

    from numpy import zeros, float32, empty
    sdevs = zeros((len(xyz),3), float32)
    for a in (0,1,2):
        sdevs[:,a] = sdev / grid.step[a]

    import Matrix as M
    if len(transforms) == 0:
        transforms = [M.identity_matrix()]
    from _gaussian import sum_of_gaussians
    ijk = empty(xyz.shape, float32)
    matrix = grid.matrix()
    for tf in transforms:
        ijk[:] = xyz
        M.transform_points(ijk, M.multiply_matrices(grid.xyz_to_ijk_transform, tf))
        sum_of_gaussians(ijk, weights, sdevs, cutoff_range, matrix)

    from math import pow, pi
    normalization = pow(2*pi,-1.5)*pow(sdev,-3)
    matrix *= normalization

def point_bounds(xyz, transforms = []):

    from _multiscale import bounding_box
    if transforms :
        from numpy import empty, float32
        xyz0 = empty((len(transforms),3), float32)
        xyz1 = empty((len(transforms),3), float32)
        txyz = empty(xyz.shape, float32)
        import Matrix as M
        for i, tf in enumerate(transforms) :
            txyz[:] = xyz
            M.transform_points(txyz, tf)
            xyz0[i,:], xyz1[i,:] = bounding_box(txyz)
        xyz_min, xyz_max = xyz0.min(axis = 0), xyz1.max(axis = 0)
    else:
        xyz_min, xyz_max = bounding_box(xyz)

    return xyz_min, xyz_max


def bounding_grid(xyz, step, pad, transforms):

    xyz_min, xyz_max = point_bounds(xyz, transforms)
    origin = [x-pad for x in xyz_min]
    from math import ceil
    shape = [int(ceil((xyz_max[a] - xyz_min[a] + 2*pad) / step)) for a in (2,1,0)]
    from numpy import zeros, float32
    matrix = zeros(shape, float32)
    from VolumeData import Array_Grid_Data
    grid = Array_Grid_Data(matrix, origin, (step,step,step))
    return grid


def MyMolMapX2 ( atoms, resolution, step=1.0, xf=None ) :

    from math import sqrt, pi
    pad = 3*resolution
    cutoff_range = 5 # in standard deviations
    sigma_factor = 1/(pi*sqrt(2)) # standard deviation / resolution
    sdev = resolution * sigma_factor

    from _multiscale import get_atom_coordinates
    xyz = get_atom_coordinates(atoms, transformed = False)
    anum = [a.element.number for a in atoms]
    grid = bounding_grid(xyz, step, pad, [])

    add_gaussians(grid, xyz, anum, sdev, cutoff_range, [])
    return grid



def GetMolPts ( mol, metric, resolution ) :

    from qscores import SetBBAts
    SetBBAts ( mol )
    atoms = None
    if 1 :
        atoms = [at for at in mol.atoms if not at.element.name == "H"]
        print " - taking all atoms only for model map"
    else :
        atoms = [at for at in mol.atoms if at.isBB == True]
        print " - taking BB atoms only for model map"

    from _multiscale import get_atom_coordinates
    points = get_atom_coordinates ( atoms, transformed = False )
    com = numpy.sum(points, axis=0) / len(points)
    comv = chimera.Vector ( *com )
    print " - center of %d atoms:" % len(atoms), com
    weights = numpy.ones ( len(points), numpy.float32 )
    print " - %d atom/pts" % len(points)

    molMap = None
    if metric == "cc" or metric == "ccm" :

        gstep = min ( resolution/3.0, 2.0 )
        #gstep = 2.0
        molData = MyMolMapX2 ( atoms, resolution, gstep, chimera.Xform.identity() )
        points, weights = fit_points_g ( molData, 0.1 )
        print " - mol to map - %d pts" % len(points)

        molMap = VolumeViewer.volume.volume_from_grid_data ( molData )
        molMap.openState.xform = mol.openState.xform
        molMap.name = mol.name + "_molMap_res%.2fA.mrc" % resolution


    return points, weights, comv, molMap




def FitPoints ( points, weights, points_, weights_, xf, inMap, doTranslate = True, doRotate = True ) :

    xfm = Matrix.xform_matrix ( xf )
    xyz_to_ijk_tf = inMap.data.xyz_to_ijk_transform
    xyz_to_ijk_tf = Matrix.multiply_matrices ( inMap.data.xyz_to_ijk_transform, xfm )
    darray = inMap.data.full_matrix()

    #map_values, outside = VolumeData.interpolate_volume_data(fpoints, xyz_to_ijk_tf, darray)
    #olap0, cc0, other = FitMap.overlap_and_correlation ( fpoint_weights, map_values )
    #print cc0,

    move_tf, stats = FitMap.locate_maximum(points, weights,
                                    darray, xyz_to_ijk_tf, max_steps = 1000,
                                    ijk_step_size_min = 0.01, ijk_step_size_max = 0.5,
                                    optimize_translation = doTranslate,
                                    optimize_rotation = doRotate,
                                    metric = 'sum product', request_stop_cb = None)

    xf = Matrix.chimera_xform ( move_tf )
    #print stats
    #score = stats['average map value']
    #score = stats['correlation about mean']
    #score = stats['correlation']

    xyz_to_ijk_tf = Matrix.multiply_matrices ( xyz_to_ijk_tf, Matrix.xform_matrix ( xf ) )
    map_values, outside = VolumeData.interpolate_volume_data(points_, xyz_to_ijk_tf, darray)
    olap, cc, ccm = FitMap.overlap_and_correlation ( weights_, map_values )
    avgd = numpy.average ( map_values )
    scores = { "avgd":avgd, "cc":cc, "ccm":ccm, "olap":olap, "z-score":-1 }

    #if not isinstance(score, (int, float)) :
    #    print " - score ? (%s)" % metric, score
    #    score = 0.0

    #ApplyXf ( ress, xf )
    #print " -> ", cc1

    return scores, xf


def addFit ( fits, newXf, newScores ) :

    newXfi = newXf.inverse()
    foundSimilar = False
    for fiti, score_xf in enumerate ( fits ) :
        scores, xf = score_xf
        xfc = xf.__copy__()
        xfc.multiply ( newXfi )
        d = xfc.getTranslation().length
        axis, angle = xfc.getRotation ()
        if d <= 3.0 and angle <= 5.0 :
            foundSimilar = True
            #if newScores["avgd"] > scores["avgd"] :
            #    fits[fiti] = [newScores, newXf]
            break

    if not foundSimilar :
        fits.append ( [newScores, newXf] )


def addFit_ ( posMap, newXf, newScores ) :

    if posMap == None :
        posMap = {}

    pos = newXf.getTranslation()
    i, j, k = int(numpy.floor(pos.x/3.0)), int(numpy.floor(pos.y/3.0)), int(numpy.floor(pos.z/3.0))

    xfs = []
    for ii in [i-1,i,i+1] :
        try :
            posMapJ = posMap[ii]
        except :
            continue

        for jj in [j-1, j, j+1] :
            try :
                posMapK = posMapJ[jj]
            except :
                continue

            for kk in [k-1,k,k+1] :
                try :
                    mapXfs = posMapK[kk]
                except :
                    continue
                xfs.extend ( mapXfs )

    foundSimilar = False
    newXfi = newXf.inverse()
    for scores, xf in xfs :
        xfc = xf.__copy__()
        xfc.multiply ( newXfi )
        d = xfc.getTranslation().length
        axis, angle = xfc.getRotation ()
        if d <= 3.0 and angle <= 5.0 :
            foundSimilar = True
            break

    if not foundSimilar :
        try :
            posMapJ = posMap[i]
        except :
            posMapJ = {}
            posMap[i] = posMapJ

        try :
            posMapK = posMapJ[j]
        except :
            posMapK = {}
            posMapJ[j] = posMapK

        try :
            posMapXfs = posMapK[k]
        except :
            posMapXfs = []
            posMapK[k] = posMapXfs

        posMapXfs.append ( [newScores, newXf] )

    return posMap



def addFit__ ( posMap, newXf, newScores ) :

    if posMap == None :
        posMap = {}

    pos = newXf.getTranslation()
    i, j, k = int(numpy.floor(pos.x/3.0)), int(numpy.floor(pos.y/3.0)), int(numpy.floor(pos.z/3.0))

    xfs = []
    for ii in [i-1,i,i+1] :
        if ii in posMap :
            posMapJ = posMap[ii]
            for jj in [j-1, j, j+1] :
                if jj in posMapJ :
                    posMapK = posMapJ[jj]
                    for kk in [k-1,k,k+1] :
                        if kk in posMapK :
                            xfs.extend ( posMapK[kk] )

    foundSimilar = False
    newXfi = newXf.inverse()
    for scores, xf in xfs :
        xfc = xf.__copy__()
        xfc.multiply ( newXfi )
        d = xfc.getTranslation().length
        axis, angle = xfc.getRotation ()
        if d <= 3.0 and angle <= 5.0 :
            foundSimilar = True
            break

    if not foundSimilar :
        posMapJ = None
        if i in posMap :
            posMapJ = posMap[i]
        else :
            posMapJ = {}
            posMap[i] = posMapJ

        posMapK = None
        if j in posMapJ :
            posMapK = posMapJ[j]
        else :
            posMapK = {}
            posMapJ[j] = posMapK

        posMapXfs = None
        if k in posMapK :
            posMapXfs = posMapK[k]
        else :
            posMapXfs = []
            posMapK[k] = posMapXfs

        posMapXfs.append ( [newScores, newXf] )

    return posMap


def fitPts ( mapPts, dmap, molPoints, molWeights, molCenter, molPoints_, molWeights_, nrot, metric, statf=None, numRandFit=None ) :

    dmapXf = dmap.openState.xform # .inverse()
    dmapXfi = dmapXf.inverse()

    fits = []
    posMap = None
    start = time.time()

    totSteps = len(mapPts) * nrot * nrot
    print " - %d map points, %d total steps" % (len(mapPts), totSteps)
    sumT, sumN = 0.0, 0.0

    spts = SpherePts ( chimera.Vector(0,0,0), 1.0, nrot )

    doPts = mapPts
    if numRandFit != None :
        if 0 :
            doPts = mapPts[0:numRandFit]
            print " - took first %d/%d map points" % (len(doPts), len(mapPts))
        else :
            from random import sample
            doPts = sample (mapPts, numRandFit)
            print " - sampled %d/%d map points" % (len(doPts), len(mapPts))

    #rot_xfs = []
    rot_xfs = None

    for p in doPts :

        #print p
        pt = dmapXf.apply ( chimera.Point(*p) )
        #xf.multiply ( chimera.Xform.translation ( pt.toVector() ) )

        for spt in spts :

            for i in range ( nrot ) :

                stepStart = time.time()

                ang = float(i) * 360.0 / float(nrot)
                xfc = chimera.Xform.translation ( -1.0 * molCenter  )
                xfc.premultiply ( chimera.Xform.rotation ( spt, ang ) )
                xfc.premultiply ( chimera.Xform.rotation ( *dmapXf.getRotation() ) )
                xfc.premultiply ( chimera.Xform.translation(chimera.Vector(*pt)) )

                xff = xfc.__copy__()
                xff.premultiply ( dmapXfi )

                if rot_xfs != None :
                    rot_xfs.append ( xfc.__copy__() )

                scores, xff = FitPoints ( molPoints, molWeights, molPoints_, molWeights_, xff, dmap )
                xfc.multiply ( xff )

                #fits.append ( [avgd, xfc] )
                #addFit ( fits, xfc, scores )
                #posMap = addFit_ ( posMap, xfc, scores )
                posMap = addFit__ ( posMap, xfc, scores )

                stepT = time.time() - stepStart
                sumT += stepT
                sumN += 1.0
                if int(sumN) % 100 == 0 :
                    stepsLeft = totSteps - int(sumN)
                    timeLeft = stepsLeft * (sumT/sumN)
                    minLeft = numpy.floor ( timeLeft / 60.0 )
                    secLeft = timeLeft - (minLeft*60.0)
                    #print "%.1f%%_%.0f:%.0f_%d" % ( 100*sumN/float(totSteps), minLeft, secLeft, len(fits) ),

                    if 0 :
                        print "%.1f%% %.0f:%.0f" % ( 100*sumN/float(totSteps), minLeft, secLeft )
                        sys.stdout.flush()
                if statf != None :
                    if int(sumN) % 10 == 0 :
                        stepsLeft = totSteps - int(sumN)
                        timeLeft = stepsLeft * (sumT/sumN)
                        minLeft = numpy.floor ( timeLeft / 60.0 )
                        secLeft = timeLeft - (minLeft*60.0)
                        fp = open ( statf, "w" )
                        fp.write ( "%.1f%% --- eta %.0f:%.0f\n" % ( 100*sumN/float(totSteps), minLeft, secLeft ) )
                        fp.close()

                #break
            #break

        if rot_xfs != None :
            print " - stopping fitPts after 1 point..."
            break

    dur = time.time() - start

    if posMap != None :
        fits = []
        for i, posMapJ in posMap.iteritems() :
            for j, posMapK in posMapJ.iteritems() :
                for k, posMapXfs in posMapK.iteritems() :
                    fits.extend ( posMapXfs )

    #print ""
    min = numpy.floor ( dur / 60.0 )
    sec = dur - (min*60.0)
    print " - took %.0f min %.0f sec" % ( min, sec)

    #if posMap != None : print " - used pos map!"

    print " - %d unique fits, top 10 by %s:" % (len(fits), metric)
    #fits.sort ( reverse=True )
    fits = sorted(fits, key=lambda x: x[0][metric], reverse=True)

    scs = []
    for scores, xf in fits[0:10] :
        #print scores[metric], xf.getTranslation(), xf.getRotation()[0], xf.getRotation()[0]
        print scores[metric]
        scs.append ( scores[metric] )

    z = (scs[0] - numpy.average(scs[1:])) / numpy.std(scs[1:])
    print " - z-score: %.1f" % z

    if rot_xfs != None :
        return fits, rot_xfs
    return fits



def fitPtsZ ( mapPts, dmap, molPoints, molWeights, molCenter, nrot, metric, statf=None, numRandFit=None ) :

    dmapXf = dmap.openState.xform # .inverse()
    dmapXfi = dmapXf.inverse()

    fits = []
    posMap = None
    start = time.time()

    totSteps = len(mapPts) * nrot * nrot
    print " - %d total steps" % totSteps
    sumT, sumN = 0.0, 0.0

    spts = SpherePts ( chimera.Vector(0,0,0), 1.0, nrot )

    doPts = mapPts
    if numRandFit != None :
        if 0 :
            doPts = mapPts[0:numRandFit]
            print " - took first %d/%d map points" % (len(doPts), len(mapPts))
        else :
            from random import sample
            doPts = sample (mapPts, numRandFit)
            print " - sampled %d/%d map points" % (len(doPts), len(mapPts))

    for p in doPts :

        #print p
        pt = dmapXf.apply ( chimera.Point(*p.pt) )
        #xf.multiply ( chimera.Xform.translation ( pt.toVector() ) )

        for spt in spts :

            posRotMap = None

            for i in range ( nrot ) :

                stepStart = time.time()

                ang = float(i) * 360.0 / float(nrot)
                xfc = chimera.Xform.translation ( -1.0 * molCenter  )
                xfc.premultiply ( chimera.Xform.rotation ( spt, ang ) )
                xfc.premultiply ( chimera.Xform.rotation ( *dmapXf.getRotation() ) )
                xfc.premultiply ( chimera.Xform.translation(chimera.Vector(*pt)) )

                xff = xfc.__copy__()
                xff.premultiply ( dmapXfi )
                scores, xff = FitPoints ( molPoints, molWeights, xff, dmap, metric )
                xfc.multiply ( xff )

                #fits.append ( [avgd, xfc] )
                #addFit ( fits, xfc, scores )
                #posMap = addFit_ ( posMap, xfc, scores )
                #posMap = addFit__ ( posMap, xfc, scores )
                posRotMap = addFit__ ( posRotMap, xfc, scores )

                stepT = time.time() - stepStart
                sumT += stepT
                sumN += 1.0
                if int(sumN) % 100 == 0 :
                    stepsLeft = totSteps - int(sumN)
                    timeLeft = stepsLeft * (sumT/sumN)
                    minLeft = numpy.floor ( timeLeft / 60.0 )
                    secLeft = timeLeft - (minLeft*60.0)
                    #print "%.1f%%_%.0f:%.0f_%d" % ( 100*sumN/float(totSteps), minLeft, secLeft, len(fits) ),
                    print "%.1f%% %.0f:%.0f " % ( 100*sumN/float(totSteps), minLeft, secLeft ),
                if statf != None :
                    if int(sumN) % 10 == 0 :
                        stepsLeft = totSteps - int(sumN)
                        timeLeft = stepsLeft * (sumT/sumN)
                        minLeft = numpy.floor ( timeLeft / 60.0 )
                        secLeft = timeLeft - (minLeft*60.0)
                        fp = open ( statf, "w" )
                        fp.write ( "%.1f%% --- eta %.0f:%.0f\n" % ( 100*sumN/float(totSteps), minLeft, secLeft ) )
                        fp.close()

                #break

            fits = []
            for i, posMapJ in posRotMap.iteritems() :
                for j, posMapK in posMapJ.iteritems() :
                    for k, posMapXfs in posMapK.iteritems() :
                        fits.extend ( posMapXfs )

            fits = sorted(fits, key=lambda x: x[0]['ccm'], reverse=True)
            scores = [scores['ccm'] for scores,xf in fits ]

            #print "scores:"
            #for sc in scores :
            #    print sc
            #print ""
            fit0 = fits[0]
            try :
                avg, std = numpy.average ( scores[1:] ), numpy.std ( scores[1:] )
                fit0[0]['z-score'] = (fit0[0]['ccm'] - avg) / std
            except :
                fit0[0]['z-score'] = 0.0

            #print "fit 0 score:", fit0[0]['ccm'], avg, std, fit0[0]['z-score']

            scores, xf = fit0
            posMap = addFit__ ( posMap, xf, scores )

            #break

        #break

    dur = time.time() - start

    if posMap != None :
        fits = []
        for i, posMapJ in posMap.iteritems() :
            for j, posMapK in posMapJ.iteritems() :
                for k, posMapXfs in posMapK.iteritems() :
                    fits.extend ( posMapXfs )

    print ""
    min = numpy.floor ( dur / 60.0 )
    sec = dur - (min*60.0)
    print " - took %.0f min %.0f sec" % ( min, sec)

    if posMap != None :
        print " - used pos map!"

    print " - %d unique fits, top 10 by %s:" % (len(fits), metric)
    #fits.sort ( reverse=True )
    fits = sorted(fits, key=lambda x: x[0]['z-score'], reverse=True)

    scs = []
    for scores, xf in fits[0:10] :
        #print scores[metric], xf.getTranslation(), xf.getRotation()[0], xf.getRotation()[0]
        print scores['ccm'], scores['z-score']
        scs.append ( scores['ccm'] )

    z = (scs[0] - numpy.average(scs[1:])) / numpy.std(scs[1:])
    print " - overall ccm z-score: %.1f" % z

    return fits


def MolBoundR ( mol ) :
    from _multiscale import get_atom_coordinates
    atoms = [at for at in mol.atoms if not at.element.name == "H"]
    points = get_atom_coordinates ( atoms, transformed = False )
    com = numpy.sum(points, axis=0) / len(points)
    comv = chimera.Vector ( *com )
    print " - center of %d atoms in mol:" % len(atoms), com
    C = chimera.Vector ( com[0], com[1], com[2] )
    comv = numpy.ones_like ( points ) * com
    points = points - comv
    bRad = numpy.sqrt ( numpy.max ( numpy.sum ( numpy.multiply ( points, points ), axis=1 ) ) )
    print " - bounding radius for mol: %.2f" % bRad
    return bRad


def fit ( dmap, mol, thr=None, step=8.0, nrot=8, metric="ccm", resolution=8.0, numRandFit=None, numOut=4, bCutoff = 80.0 ) :

    print "Fit %s -> %s" % (mol.name, dmap.name)

    if not hasattr ( dmap, 'surface_levels' ) :
        print " - fitting to a molecule, converting it to a map first"
        dmol = dmap

        from qscores import SetBBAts
        SetBBAts ( dmol )
        #atoms = [at for at in mol.atoms if not at.element.name == "H"]
        atoms = [at for at in dmol.atoms if at.isBB == True]

        dmapRes, dmapStep = 7.0, 2.0
        molData = MyMolMapX2 ( atoms, dmapRes, dmapStep, chimera.Xform.identity() )
        dmap = VolumeViewer.volume.volume_from_grid_data ( molData )
        dmap.openState.xform = dmol.openState.xform
        dmap.name = dmol.name + "_molMap_res%.2fA.mrc" % dmapRes
        thr = None

    print " - in map: %s, thr %.2f" % (dmap.name, dmap.surface_levels[0])
    print " - mol: %s, %d atoms" % ( mol.name, len(mol.atoms) )
    print " - step: %.2f" % step
    print " - resolution: %.2f" % resolution
    print " - num rotations: %.2f x %.2f = %.2f, angle %.2f" % (nrot, nrot, nrot*nrot, 360/float(nrot))
    if thr != None :
        print " - threshold: %.2f" % thr
    else :
        M = dmap.data.full_matrix()
        thr = numpy.average(M)+numpy.std(M)*2.0
        print " - threshold at 2-sigma: %.2f" % thr

    if bCutoff != None :
        bCutoff = 80.0
        nmol = CopyAfMol ( mol, chimera.Xform.identity(), bCutoff=bCutoff )
        nmol.name = os.path.splitext ( mol.name )[0] + ("_af%.0f.pdb" % bCutoff)
        chimera.openModels.add ( [nmol] )
        print " - af mol: %s, %d atoms" % ( nmol.name, len(nmol.atoms) )
        #mfile = os.path.join(fitsPath, nmol.name )
        #print " - saving mol %s" % mfile
        #chimera.PDBio().writePDBfile ( [nmol], mfile )
        mol = nmol


    # molecule for fitting - use lower resolution for less points and avoiding local minima
    molPoints, molWeights, molCenter, molMap = GetMolPts ( mol, metric, 7.0 )
    print " - fitting with %d points with metric %s" % (len(molPoints), metric)

    # molecule for scoring - use actual resolution specified
    molPoints_, molWeights_, molCenter_, molMap_ = GetMolPts ( mol, metric, resolution )
    print " - scoring with %d points with metric %s" % (len(molPoints_), metric)

    mapPts = GetMapPts ( dmap, thr, step, numMax=None )
    #bRad = MolBoundR ( mol )
    #mapPts = GetMapGridPts ( dmap, bRad=bRad/2.0, thr=None, step=step, numMax=None, showPts=True )

    #chimera.openModels.close ( [molMap, molMap_] )
    #return

    mapPts = mapPts [ len(mapPts)/3 : len(mapPts)/3 + 3 ]

    fits, rxfs = fitPts ( mapPts, dmap, molPoints, molWeights, molCenter, molPoints_, molWeights_, nrot, metric, statf=None ), None
    #fits, rxfs = fitPts ( mapPts, dmap, molPoints, molWeights, molCenter, molPoints_, molWeights_, nrot, metric, statf=None )

    if rxfs != None :
        print ( " - showing %d rotations..." % len(rxfs) )
        for rxf in rxfs :
            cmol = CopyMolX ( mol, chimera.Xform() )
            chimera.openModels.add ( [cmol] )
            cmol.openState.xform = rxf

    scores, xf = fits[0]
    print " - putting %s to top fit, %s=%.3f" % (mol.name, metric, scores[metric])
    mol.openState.xform = xf
    if molMap :
        molMap.openState.xform = xf
        molMap_.openState.xform = xf



def fitm ( dmap, toMap, dmapThr=None, toMapThr=None, step=8.0, useGrid=False, nrot=8, metric="ccm", numRandFit=None, numOut=4 ) :

    print "Fit map %s -> %s" % (dmap.name, toMap.name)

    if not hasattr ( dmap, 'surface_levels' ) :
        print " - fitting to a molecule, converting it to a map first"
        dmol = dmap

        from qscores import SetBBAts
        SetBBAts ( dmol )
        #atoms = [at for at in mol.atoms if not at.element.name == "H"]
        atoms = [at for at in dmol.atoms if at.isBB == True]

        dmapRes, dmapStep = 7.0, 2.0
        molData = MyMolMapX2 ( atoms, dmapRes, dmapStep, chimera.Xform.identity() )
        dmap = VolumeViewer.volume.volume_from_grid_data ( molData )
        dmap.openState.xform = dmol.openState.xform
        dmap.name = dmol.name + "_molMap_res%.2fA.mrc" % dmapRes
        thr = None

    print " - to map: %s, thr %.2f" % (toMap.name, toMap.surface_levels[0])
    print " - fit map: %s" % ( dmap.name )
    print " - step: %.2f" % step
    print " - num rotations: %.2f x %.2f = %.2f, angle %.2f" % (nrot, nrot, nrot*nrot, 360/float(nrot))

    if dmapThr != None :
        print " - %s threshold: %.2f" % (dmap.name, dmapThr)
    else :
        M = dmap.data.full_matrix()
        dmapThr = numpy.average(M)+numpy.std(M)*2.0
        print " - %s threshold at 2-sigma: %.2f" % (dmap.name, dmapThr)
        chimera.runCommand ( "vol #%d level %f" % (dmap.id, dmapThr) )

    if toMapThr != None :
        print " - %s threshold: %.2f" % (toMap.name, toMapThr)
    else :
        M = toMap.data.full_matrix()
        toMapThr = numpy.average(M)+numpy.std(M)*2.0
        print " - %s threshold at 2-sigma: %.2f" % (toMap.name, toMapThr)
        chimera.runCommand ( "vol #%d level %f" % (toMap.id, toMapThr) )

    fitPoints, fitWeights = fit_points_g ( dmap.data, dmapThr )
    print " - %s - %d fit pts" % (dmap.name, len(fitPoints))
    if len(fitPoints) > 1e4 :
        print " - too many fit points - try a higher threshold or vop bin the map"
        return

    com = numpy.sum(fitPoints, axis=0) / len(fitPoints)
    comv = chimera.Vector ( *com )
    print " - center of %s: %f,%f,%f" % ( dmap.name, com[0], com[1], com[2] )

    mapPts = None
    if not useGrid :
        mapPts = GetMapPts ( toMap, toMapThr, step, numMax=None )
        mapPts = mapPts [ len(mapPts)/2 : len(mapPts)/2 + 2 ]
    else :
        mapPts = GetMapPtsGrid ( toMap, thr=toMapThr, step=step )
        mapPts = mapPts [ len(mapPts)/2 : len(mapPts)/2 + 2 ]

    fits, rxfs = fitPts ( mapPts, toMap, fitPoints, fitWeights, comv, fitPoints, fitWeights, nrot, metric, statf=None ), None
    if rxfs != None :
        print ( " - showing %d rotations..." % len(rxfs) )
        for rxf in rxfs :
            cmol = CopyMolX ( mol, chimera.Xform() )
            chimera.openModels.add ( [cmol] )
            cmol.openState.xform = rxf

    scores, xf = fits[0]
    print " - putting %s to top fit, %s=%.3f" % (dmap.name, metric, scores[metric])
    dmap.openState.xform = xf

    savePath = os.path.splitext(dmap.data.path)[0] + "_rf.mrc"
    MapToMap ( dmap, dmapThr, toMap, nname="fit", saveToPath = savePath )



def MapToMap ( dmap, thr, inMap, nname="fit", saveToPath=None ) :

    fitPoints, fitWeights = fit_points_g ( dmap.data, thr )
    #print " - %s - %d fit pts" % (dmap.name, len(fitPoints))

    nv = MapToMap0 ( fitPoints, inMap, dmap, nname )

    if saveToPath != None :
        #print " --> %s" % saveToPath
        nv.write_file ( saveToPath, "mrc" )





def MapToMap0 ( points, dmap, fromMap, nname, showMesh = False, alpha=0.4 ) :

    #_contour.affine_transform_vertices ( points, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )
    #mdata = VolumeData.zone_masked_grid_data ( dmap.data, points, atomRad )

    import _contour
    points1 = numpy.copy ( points )
    _contour.affine_transform_vertices ( points1, Matrix.xform_matrix( fromMap.openState.xform ) )
    _contour.affine_transform_vertices ( points1, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )
    points0 = numpy.copy ( points1 )
    _contour.affine_transform_vertices ( points1, dmap.data.xyz_to_ijk_transform )

    bound = 5
    li,lj,lk = numpy.min ( points1, axis=0 ) - (bound, bound, bound)
    hi,hj,hk = numpy.max ( points1, axis=0 ) + (bound, bound, bound)

    n1 = hi - li + 1
    n2 = hj - lj + 1
    n3 = hk - lk + 1

    #print " - bounds - %d %d %d --> %d %d %d --> %d %d %d" % ( li,lj,lk, hi,hj,hk, n1,n2,n3 )

    #nmat = numpy.zeros ( (n1,n2,n3), numpy.float32 )
    #dmat = dmap.full_matrix()

    nstep = (dmap.data.step[0], dmap.data.step[1], dmap.data.step[2] )
    #nstep = (fmap.data.step[0]/2.0, fmap.data.step[1]/2.0, fmap.data.step[2]/2.0 )

    nn1 = int ( round (dmap.data.step[0] * float(n1) / nstep[0]) )
    nn2 = int ( round (dmap.data.step[1] * float(n2) / nstep[1]) )
    nn3 = int ( round (dmap.data.step[2] * float(n3) / nstep[2]) )

    O = dmap.data.origin
    #print " - %s origin:" % dmap.name, O
    nO = ( O[0] + float(li) * dmap.data.step[0],
           O[1] + float(lj) * dmap.data.step[1],
           O[2] + float(lk) * dmap.data.step[2] )

    #print " - new map origin:", nO

    nmat = numpy.zeros ( (nn1,nn2,nn3), numpy.float32 )
    ndata = VolumeData.Array_Grid_Data ( nmat, nO, nstep, dmap.data.cell_angles )

    #print " - fmap grid dim: ", numpy.shape ( fmap.full_matrix() )
    #print " - new map grid dim: ", numpy.shape ( nmat )

    npoints = VolumeData.grid_indices ( (nn1, nn2, nn3), numpy.single)  # i,j,k indices
    _contour.affine_transform_vertices ( npoints, ndata.ijk_to_xyz_transform )
    #_contour.affine_transform_vertices ( npoints, fromMap.openState.xform.inverse() )

    #dvals = fromMap.interpolated_values ( npoints, fromMap.openState.xform )
    dvals = fromMap.interpolated_values ( npoints, dmap.openState.xform )
    #dvals = fromMap.interpolated_values ( npoints, fromMap.openState.xform.inverse() )
    #dvals = fromMap.interpolated_values ( npoints, chimera.Xform.identity() )
    #dvals = fromMap.interpolated_values ( npoints, dmap.openState.xform.inverse() )
    #dvals = numpy.where ( dvals > threshold, dvals, numpy.zeros_like(dvals) )
    #nze = numpy.nonzero ( dvals )

    nmat = dvals.reshape( (nn3,nn2,nn1) )
    #f_mat = fmap.data.full_matrix()
    #f_mask = numpy.where ( f_mat > fmap.surface_levels[0], numpy.ones_like(f_mat), numpy.zeros_like(f_mat) )
    #df_mat = df_mat * f_mask

    ndata = VolumeData.Array_Grid_Data ( nmat, nO, nstep, dmap.data.cell_angles, name=nname )
    #try : nv = VolumeViewer.volume.add_data_set ( ndata, None )
    #except : nv = VolumeViewer.volume.volume_from_grid_data ( ndata )

    nv = VolumeViewer.volume.volume_from_grid_data ( ndata )
    nv.openState.xform = dmap.openState.xform

    nv.name = nname
    #dmap.display = False
    nv.region = ( nv.region[0], nv.region[1], [1,1,1] )
    nv.surface_levels[0] = fromMap.surface_levels[0]
    chimera.runCommand ( "vol #%d level %f" % (nv.id, fromMap.surface_levels[0]) )

    return nv



def fitProc () :

    maps = chimera.openModels.list ( modelTypes = [VolumeViewer.volume.Volume] )
    mols = chimera.openModels.list ( modelTypes = [chimera.Molecule] )

    print "%d maps, %d mols" % (len(maps), len(mols))

    if len(maps) < 0 :
        return
    if len(mols) < 0 :
        return

    dmap = maps[0]
    mol = mols[0]

    print " - map: %s" % dmap.data.path
    print " - mol: %s" % mol.openedAs[0]

    molPath, molFile = os.path.split ( mol.openedAs[0] )
    molName, molExt = os.path.splitext ( molFile )

    procI = molName.split("__")[-1]

    print ""
    print "proc %s" % procI

    procPointsPath = os.path.join ( molPath, "%s_points.txt" % procI )
    print " - reading %s" % procPointsPath
    fin = open ( procPointsPath, "r" )

    l1 = fin.readline(); thr = float(l1); print " - thr: %.2f" % thr
    l1 = fin.readline(); step = float(l1); print " - step: %.2f" % step
    l1 = fin.readline(); metric = l1.strip(); print " - metric: %s" % metric
    l1 = fin.readline(); resolution = float(l1); print " - resolution: %.2f" % resolution
    l1 = fin.readline(); nrot = int(l1); print " - num rot: %d" % nrot

    #mapPts = GetMapPts ( dmap, thr, step, numMax=None )
    mapPts = []
    for l in fin :
        ts = l.split()
        #p = lambda: None
        p = [ float(ts[0]), float(ts[1]), float(ts[2]) ]
        mapPts.append ( p )
    fin.close()

    print " - read %d map points" % len(mapPts)

    # molecule for fitting - use lower resolution for efficiency and avoiding local minima
    #molPoints, molWeights, molCenter, molMap = GetMolPts ( mol, metric, 7.0 )

    # use just atom points for fitting - there are less, and get similar results than using
    # points from model-map with weights (highest weights are at atom position)
    molPoints, molWeights, molCenter, molMap = GetMolPts ( mol, "", 7.0 )
    print " - fitting with %d points with metric %s" % (len(molPoints), metric)

    # molecule for scoring - use actual resolution specified
    molPoints_, molWeights_, molCenter_, molMap_ = GetMolPts ( mol, metric, resolution )
    print " - scoring with %d points with metric %s" % (len(molPoints_), metric)

    sys.stdout.flush()

    statf = os.path.join ( molPath, "%s_stat.txt" % procI )

    fits = fitPts ( mapPts, dmap, molPoints, molWeights, molCenter, molPoints_, molWeights_, nrot, metric, statf=statf )

    print " - writing %d top unique fits" % len(fits)

    outp = os.path.join ( molPath, "%s_fits.txt" % procI )
    print " -> %s" % outp
    #sys.stdout.flush()

    fout = open ( outp, "w" )
    for scores, xf in fits :
        xf.premultiply ( dmap.openState.xform.inverse() )
        v, ang = xf.getRotation()
        t = xf.getTranslation()
        fout.write ( "%f %f %f %f %f %f %f %f %f %f %f %f\n" % (scores["avgd"], scores["cc"], scores["ccm"], scores["olap"], scores["z-score"], ang, v[0], v[1], v[2], t[0], t[1], t[2]) )
        #print ".fit. %f %f %f %f %f %f %f %f %f %f %f %f\n" % (scores["avgd"], scores["cc"], scores["ccm"], scores["olap"], scores["z-score"], ang, v[0], v[1], v[2], t[0], t[1], t[2])
    fout.close()

    time.sleep ( 1.0 )

    print ".done."
    sys.stdout.flush()


# called by fitMPm.py - multiproc process
def fitProcM () :

    maps = chimera.openModels.list ( modelTypes = [VolumeViewer.volume.Volume] )

    print "%d maps" % ( len(maps) )

    if len(maps) < 2 :
        print " - need at least two maps open"
        return

    toMap, dmap = maps[0], maps[1]

    print " - toMap: %s" % toMap.data.path
    print " - map: %s" % dmap.data.path

    toMapPath, toMapFile = os.path.split ( toMap.data.path )
    mapPath, mapFile = os.path.split ( dmap.data.path )
    mapName = os.path.splitext ( mapFile )[0]
    procI = mapName.split("_")[-1]

    print ""
    print "proc %s" % procI

    procPointsPath = os.path.join ( mapPath, "%s_points.txt" % procI )
    print " - reading %s" % procPointsPath
    fin = open ( procPointsPath, "r" )

    l1 = fin.readline(); toMapThr = float(l1); print " - to map thr: %.2f" % toMapThr
    l1 = fin.readline(); dmapThr = float(l1); print " - map thr: %.2f" % dmapThr
    l1 = fin.readline(); step = float(l1); print " - step: %.2f" % step
    l1 = fin.readline(); metric = l1.strip(); print " - metric: %s" % metric
    l1 = fin.readline(); nrot = int(l1); print " - num rot: %d" % nrot

    #mapPts = GetMapPts ( dmap, thr, step, numMax=None )
    mapPts = []
    for l in fin :
        ts = l.split()
        #p = lambda: None
        p = [ float(ts[0]), float(ts[1]), float(ts[2]) ]
        mapPts.append ( p )
    fin.close()

    print " - read %d map points" % len(mapPts)

    fitPoints, fitWeights = fit_points_g ( dmap.data, dmapThr )
    print " - %s - %d fit pts" % (dmap.name, len(fitPoints))

    com = numpy.sum(fitPoints, axis=0) / len(fitPoints)
    comv = chimera.Vector ( *com )
    print " - center of %s: %f,%f,%f" % ( dmap.name, com[0], com[1], com[2] )

    sys.stdout.flush()

    statf = os.path.join ( mapPath, "%s_stat.txt" % procI )

    fits = fitPts ( mapPts, toMap, fitPoints, fitWeights, comv, fitPoints, fitWeights, nrot, metric, statf=statf )

    print " - writing %d unique fits" % len(fits)

    outp = os.path.join ( mapPath, "%s_fits.txt" % procI )
    print " -> %s" % outp
    #sys.stdout.flush()

    fout = open ( outp, "w" )
    for scores, xf in fits :
        xf.premultiply ( toMap.openState.xform.inverse() )
        v, ang = xf.getRotation()
        t = xf.getTranslation()
        fout.write ( "%f %f %f %f %f %f %f %f %f %f %f %f\n" % (scores["avgd"], scores["cc"], scores["ccm"], scores["olap"], scores["z-score"], ang, v[0], v[1], v[2], t[0], t[1], t[2]) )
        #print ".fit. %f %f %f %f %f %f %f %f %f %f %f %f\n" % (scores["avgd"], scores["cc"], scores["ccm"], scores["olap"], scores["z-score"], ang, v[0], v[1], v[2], t[0], t[1], t[2])
    fout.close()

    time.sleep ( 1.0 )

    print ".done."
    sys.stdout.flush()






def fitProc_pipe () :

    maps = chimera.openModels.list ( modelTypes = [VolumeViewer.volume.Volume] )
    mols = chimera.openModels.list ( modelTypes = [chimera.Molecule] )

    print "%d maps, %d mols" % (len(maps), len(mols))

    if len(maps) < 0 :
        return
    if len(mols) < 0 :
        return

    dmap = maps[0]
    mol = mols[0]

    print " - map: %s" % dmap.data.path
    print " - mol: %s" % mol.openedAs[0]

    molPath, molFile = os.path.split ( mol.openedAs[0] )
    molName, molExt = os.path.splitext ( molFile )

    procI = molName.split("__")[-1]

    print ""
    print "proc %s" % procI

    procPointsPath = os.path.join ( molPath, "%s_points.txt" % procI )
    print " - reading %s" % procPointsPath
    fin = open ( procPointsPath, "r" )

    l1 = fin.readline(); thr = float(l1); print " - thr: %.2f" % thr
    l1 = fin.readline(); step = float(l1); print " - step: %.2f" % step
    l1 = fin.readline(); metric = l1.strip(); print " - metric: %s" % metric
    l1 = fin.readline(); resolution = float(l1); print " - resolution: %.2f" % resolution
    l1 = fin.readline(); nrot = int(l1); print " - num rot: %d" % nrot

    #mapPts = GetMapPts ( dmap, thr, step, numMax=None )
    mapPts = []
    for l in fin :
        ts = l.split()
        #p = lambda: None
        p = [ float(ts[0]), float(ts[1]), float(ts[2]) ]
        mapPts.append ( p )
    fin.close()

    print " - read %d map points" % len(mapPts)

    # molecule for fitting - use lower resolution for efficiency and avoiding local minima
    #molPoints, molWeights, molCenter, molMap = GetMolPts ( mol, metric, 7.0 )

    # use just atom points for fitting - there are less, and get similar results than using
    # points from model-map with weights (highest weights are at atom position)
    molPoints, molWeights, molCenter, molMap = GetMolPts ( mol, "", 7.0 )
    print " - fitting with %d points with metric %s" % (len(molPoints), metric)

    # molecule for scoring - use actual resolution specified
    molPoints_, molWeights_, molCenter_, molMap_ = GetMolPts ( mol, metric, resolution )
    print " - scoring with %d points with metric %s" % (len(molPoints_), metric)

    sys.stdout.flush()

    statf = os.path.join ( molPath, "%s_stat.txt" % procI )

    fits = fitPts ( mapPts, dmap, molPoints, molWeights, molCenter, molPoints_, molWeights_, nrot, metric, statf=statf )

    #print " - writing %d top unique fits" % len(fits)
    #sys.stdout.flush()

    #outp = os.path.join ( molPath, "%s_fits.txt" % procI )
    #print " -> %s" % outp
    #sys.stdout.flush()

    #fout = open ( outp, "w" )
    for scores, xf in fits :
        xf.premultiply ( dmap.openState.xform.inverse() )
        v, ang = xf.getRotation()
        t = xf.getTranslation()
        #fout.write ( "%f %f %f %f %f %f %f %f %f %f %f %f\n" % (scores["avgd"], scores["cc"], scores["ccm"], scores["olap"], scores["z-score"], ang, v[0], v[1], v[2], t[0], t[1], t[2]) )
        print ".fit. %f %f %f %f %f %f %f %f %f %f %f %f\n" % (scores["avgd"], scores["cc"], scores["ccm"], scores["olap"], scores["z-score"], ang, v[0], v[1], v[2], t[0], t[1], t[2])
        #sys.stdout.flush()
    #fout.close()

    print ".done."
    sys.stdout.flush()






def GetChimeraPath () :

    # '/Users/greg/_mol/Chimera.app/Contents/Resources/share/__main__.py'
    import sys
    chimeraPath = os.path.split ( sys.argv[0] )[0]

    chimeraPath, share = os.path.split ( chimeraPath )
    chimeraPath = os.path.join ( chimeraPath, 'bin' )
    chimeraPath = os.path.join ( chimeraPath, 'chimera' )
    if os.path.isfile ( chimeraPath ) :
        #print " -- on unix/mac"
        pass
    else :
        chimeraPath += ".exe"
        if os.path.isfile ( chimeraPath ) :
            #print " -- on windows"
            pass
        else :
            print " - chimera path not found..."
            return ""

    return chimeraPath



def fitMP ( dmap, mol, numProc=None, thr=None, step=8.0, nrot=8, metric="avgd", resolution=8.0, bCutoff=None, numOut=5 ) :

    print "Fit MP %s -> %s" % (mol.name, dmap.name)

    if not hasattr ( dmap, 'surface_levels' ) :
        print " - fitting to a molecule, converting it to a map first"
        dmol = dmap

        from qscores import SetBBAts
        SetBBAts ( dmol )
        #atoms = [at for at in mol.atoms if not at.element.name == "H"]
        atoms = [at for at in dmol.atoms if at.isBB == True]

        dmapRes, dmapStep = 7.0, 2.0
        molData = MyMolMapX2 ( atoms, dmapRes, dmapStep, chimera.Xform.identity() )
        dmap = VolumeViewer.volume.volume_from_grid_data ( molData )
        dmap.openState.xform = dmol.openState.xform
        dmap.name = dmol.name + "_molMap_res%.2fA.mrc" % dmapRes
        thr = None

        npath = os.path.split(dmol.openedAs[0])[0] + "/" + dmap.name
        #print " -> ", npath
        dmap.write_file ( npath, "mrc" )
        print " -> ", dmap.data.path


    print " - in map: %s, thr %.2f" % (dmap.name, dmap.surface_levels[0])
    print " - mol: %s, %d atoms" % ( mol.name, len(mol.atoms) )
    print " - step: %.2f" % step
    print " - resolution: %.2f" % resolution
    print " - num rotations: %.2f x %.2f = %.2f, angle %.2f" % (nrot, nrot, nrot*nrot, 360/float(nrot))
    if thr != None :
        print " - threshold: %.2f" % thr
    else :
        M = dmap.data.full_matrix()
        thr = numpy.average(M)+numpy.std(M)*2.0
        print " - threshold at 2-sigma: %.2f" % thr

    #inPath = os.path.dirname(os.path.realpath(__file__))
    #print inPath
    #print __file__

    chimeraPath = GetChimeraPath ()
    print " - chimera path: %s" % chimeraPath

    dir_path = os.path.dirname(os.path.realpath(__file__))
    inDir = os.path.split(dir_path)[0]
    print " -- working dir:", inDir
    #mapQPPath = os.path.join ( inDir, 'Segger' )
    fitScriptPath = os.path.join ( dir_path, 'fitMP.py' )
    print " -- path to fitScriptPath script:", fitScriptPath

    if numProc == None :
        import multiprocessing
        numProc = multiprocessing.cpu_count() / 2
        #numProc = 2

    print " - num processors: %d" % numProc

    molName0 = os.path.splitext ( mol.name )[0]

    weight0, weight1 = 0.0, 0.0

    if bCutoff != None :
        nmol = CopyAfMol ( mol, chimera.Xform.identity(), bCutoff=bCutoff )

        chimera.openModels.add ( [nmol] )

        print " - addh for better molecular weight estimation"
        try :
            chimera.runCommand ( "addh" )
        except :
            print " - addh failed, that's ok"

        for at in mol.atoms : weight0 += at.element.number * 2.0
        for at in nmol.atoms : weight1 += at.element.number * 2.0
        weight0 /= 1000.0
        weight1 /= 1000.0

        nmol.name = os.path.splitext ( mol.name )[0] + "__%.2fkDa__af%.0f__%d__%.2fkDa__.pdb" % ( weight0, bCutoff, len(nmol.residues), weight1 )
        #mfile = os.path.join(fitsPath, nmol.name )
        #print " - saving mol %s" % mfile
        #chimera.PDBio().writePDBfile ( [nmol], mfile )
        npath = os.path.split(mol.openedAs[0])[0] + nmol.name
        nmol.openedAs = [ npath, []]
        print " - af mol: %s  (%d res - %.2f kDa), cutoff %.2f - (%d res, %.2f kDa)" % ( nmol.name, len(mol.residues), weight0, bCutoff, len(nmol.residues), weight1 )
        mol = nmol

    molPath = os.path.splitext(mol.openedAs[0])[0]
    molName = os.path.splitext(mol.name)[0]
    mapName = os.path.splitext(dmap.name)[0]
    mapPath = os.path.split ( dmap.data.path )[0]
    mapBase = os.path.splitext ( dmap.data.path )[0]

    def checkPath ( path ) :
        if not os.path.isdir ( path ) :
            try :
                os.mkdir(path)
            except :
                print " - could not make: %s" % path
                return

    fitsPath = mapBase + "__fits_res%.2f_step%.0f_nrot%d_bf%.0f" % (resolution, step, nrot, (bCutoff if bCutoff != None else 0))
    print " - fits path: %s" % fitsPath
    if checkPath (fitsPath) == False : return

    allFitsPath = fitsPath + "/all_fits"
    print " - all fits path: %s" % allFitsPath
    if checkPath (allFitsPath) == False : return

    if len(mol.residues) == 0 :
        print " - no residues in molecule?"
        fp = open ( allFitsPath + "/" + molName0 + "_all_fits.txt", "w" )
        fp.write ( " - no residues to fit\n" )
        fp.close()
        return


    # make path for temporary results
    from random import choice
    from string import ascii_uppercase, digits
    randStr = ''.join(choice(ascii_uppercase + digits) for _ in range(7))
    #tempPath = mapBase + "__ez_fit_temp__%s__" % randStr
    tempPath = mapBase + "__fitMP__res%.2f_step%.0f_nrot%d__%s" % (resolution, step, nrot, molName0)

    start0 = time.time()

    procs = []

    if os.path.isdir ( tempPath ) :
        print " - temp path exists, removing: %s" % tempPath
        from shutil import rmtree
        rmtree( tempPath )

    print " - making temp path: %s" % tempPath

    try :
        os.mkdir(tempPath)
    except :
        print " - could not make temp path"
        print "    : %s" % tempPath
        print "    - check/remove temp path manually and try again"
        print "    - or, check write permissions"
        return


    #molPoints, molWeights, molCenter, molMap = GetMolPts ( mol, metric, resolution )
    #print " - fitting %d points with metric %s" % (len(molPoints), metric)

    mapPts = GetMapPts ( dmap, thr, step, numMax=None )
    #mapPts = mapPts[0:100]

    #fits = fitPts ( mapPts, dmap, molPoints, molWeights, molCenter, nrot, metric )

    # ------------------

    n = len(mapPts)
    g = [mapPts[(n*c)/numProc:(n*(c+1))/numProc] for c in range(numProc)]
    for mi, mpts in enumerate(g) :

        #mpts = mpts[0:2]
        print " - proc %d/%d, %d map points" % (mi+1, numProc, len(mpts))

        procPointsPath = os.path.join ( tempPath, "%d_points.txt" % mi )
        fout = open ( procPointsPath, "w" )
        #fout.write ( "%s\n" % dmap.data.path )
        fout.write ( "%.2f\n" % thr )
        fout.write ( "%.2f\n" % step )
        fout.write ( "%s\n" % metric )
        fout.write ( "%.2f\n" % resolution )
        fout.write ( "%d\n" % nrot )
        for mpt in mpts :
            fout.write ( "%f %f %f\n" % (mpt[0], mpt[1], mpt[2]) )
        fout.close()

        #nmap_path = os.path.join ( tempPath, "%d_map.mrc" % mi )
        nmol_path = os.path.join ( tempPath, molName + "__%d.pdb" % mi )

        if 0 :
            nmap = MaskMapResize ( atoms1, 6, dmap, nmap_path )
        else :
            #import shutil
            #shutil.copyfile ( dmap.data.path, nmap_path )
            print " - writing mol -> %s" % nmol_path
            #shutil.copyfile ( mol.openedAs[0], nmol_path )
            chimera.PDBio().writePDBfile ( [mol], nmol_path )

        args = [chimeraPath, '--nogui', '--silent', '--nostatus', dmap.data.path, nmol_path, fitScriptPath]
        if mi == 0 :
            print "  >" + " ".join(args)

        fout = open ( os.path.join(tempPath, "%d_out.txt" % mi), "w" )
        foute = open ( os.path.join(tempPath, "%d_err.txt" % mi), "w" )

        inPath = "/Users/greg/Dropbox/_mol"
        if os.path.exists ( "/storage/toxop" ) : inPath="/storage/toxop"
        elif os.path.exists ( "/scratch/gregp/toxop" ) : inPath="/scratch/gregp/toxop"
        elif os.path.exists ( "/sdf/home/g/gregp/g24" ) : inPath="/sdf/home/g/gregp/g24"
        elif os.path.exists ( "/storage/gregp/toxop" ) : inPath="/storage/gregp/toxop"
        elif os.path.exists ( "/sdf/home/g/gregp/g24/ez" ) : inPath="/sdf/home/g/gregp/g24/ez"
        inPath = "./"
        print " - inPath: %s" % inPath

        from subprocess import Popen
        import subprocess
        p = Popen(args, stdout=fout, stderr=foute, cwd=inPath)
        procs.append ( [mi, p, fout, foute] )

        #p = Popen(args, stdout=subprocess.PIPE, stderr=foute, cwd=inPath)
        #procs.append ( [mi, p, foute] )
        #p.done = False


    print ""
    print "Waiting for fits..."
    print ""
    posMap = None
    totFits = 0
    error = False
    for mi, p, fout, foute in procs :
        p.wait()
        fout.close()
        foute.close()

        fitsf = os.path.join(tempPath, "%d_fits.txt" % mi)
        print "%d -- %s" % (mi, fitsf)

        if not os.path.isfile ( fitsf ) :
            print " -x- fits file not found -x-"
            error = True
            continue

        ff = open ( fitsf )
        fiti = 0
        for l in ff :
            ll = l.strip()
            t = ll.split()
            avgd, cc, ccm, olap, z, ang = float(t[0]), float(t[1]), float(t[2]), float(t[3]), float(t[4]), float(t[5])
            scores = {"avgd":avgd, "cc":cc, "ccm":ccm, "olap":olap, "z-score":z}
            v = chimera.Vector ( float(t[6]), float(t[7]), float(t[8]) ); v.normalize()
            t = chimera.Vector ( float(t[9]), float(t[10]), float(t[11]) )
            xf = chimera.Xform.translation ( t )
            xf.multiply ( chimera.Xform.rotation ( v, ang ) )
            xf.premultiply ( dmap.openState.xform )
            #addFit ( fits, xf, scores )
            posMap = addFit__ ( posMap, xf, scores )
            totFits += 1
            fiti += 1
        ff.close()
        print "%d -- %d fits" % (mi, fiti)

    print ""


    if 0 :
        print " - removing temp files"
        for mi, p, fout in procs :
            #print " - removing temp files for proc %d" % mi
            #os.remove ( os.path.join(tempPath, "%d_fits.txt" % mi) )
            os.remove ( os.path.join(tempPath, "%d_points.txt" % mi) )
            #os.remove ( os.path.join(tempPath, "%d_map.mrc" % mi) )
            #os.remove ( os.path.join(tempPath, "%d_out.txt" % mi) )
            os.remove ( os.path.join(tempPath, "%d_err.txt" % mi) )
            os.remove ( os.path.join(tempPath, molName + "__%d.pdb" % mi ) )

            statf = os.path.join(tempPath, "%d_stat.txt" % mi)
            if os.path.isfile ( statf ) : os.remove ( statf )

        os.rmdir ( tempPath )

    if posMap == None :
        print " -- no fits??? --"
        return

    if 1 :
        print " - removing temp path: %s" % tempPath
        from shutil import rmtree
        rmtree ( tempPath )

    if error :
        print " -- stopping with error -- "
        return

    fits = []
    for i, posMapJ in posMap.iteritems() :
        for j, posMapK in posMapJ.iteritems() :
            for k, posMapXfs in posMapK.iteritems() :
                fits.extend ( posMapXfs )

    print ""
    print " - %d unique fits from %d total" % (len(fits), totFits)

    #fitsPath = mapBase + "__fits_res%.2f_step%.0f_nrot%d_bf%.0f" % (resolution, step, nrot, bCutoff)
    print " - fits path: %s" % fitsPath
    if checkPath (fitsPath) == False : return

    allFitsPath = fitsPath + "/all_fits"
    print " - all fits path: %s" % allFitsPath
    if checkPath (allFitsPath) == False : return

    fp = open ( allFitsPath + "/" + molName0 + "_all_fits.txt", "w" )
    fp.write ( "Avg.Density\tZ(Avg.D.)\tCC\tZ(CC)\tCC-mean\tZ(CC-mean)\tAngle\tRotAx_X\tRotAx_Y\tRotAx_Z\tTrans_X\tTrans_Y\tTrans_Z\n" )

    if len(fits) < 1 :
        print " - no fits found... probably no map points above threshold?"
        fp.close()
        return

    avg, std = {}, {}
    for tp in ['avgd', 'cc', 'ccm'] :
        scs = [scores[tp] for scores, xf in fits]
        avg[tp] = numpy.average ( scs )
        std[tp] = numpy.std ( scs )

    print " - sorting fits..."
    fits = sorted(fits, key=lambda x: x[0]['ccm'], reverse=True)

    for scores, xf in fits :
        v, ang = xf.getRotation()
        t = xf.getTranslation()
        avgd, cc, ccm = scores['avgd'], scores['cc'], scores['ccm']
        # scores['olap'], scores['z-score']
        avgd_z = (avgd - avg['avgd']) / std['avgd']
        cc_z =   (cc -   avg['cc'])   / std['cc']
        ccm_z =  (ccm -  avg['ccm'])  / std['ccm']
        fp.write ( "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (avgd, avgd_z, cc, cc_z, ccm, ccm_z , ang, v[0], v[1], v[2], t[0], t[1], t[2]) )

        scores['avgd_z'] = avgd_z
        scores['cc_z'] = cc_z
        scores['ccm_z'] = ccm_z

    fp.close()


    #for mt in ["avgd", "cc", "ccm"] :
    #for mt in [ "cc", "ccm", "olap" ] :
    #for mt in [ "ccm", "cc", 'z-score' ] :
    for mt in [ "ccm" ] :

        fitsPath_ = fitsPath + "/by_%s" % mt
        if checkPath (fitsPath_) == False : return
        fitsPathN_ = fitsPath + "/by_%s_N_%d" % (mt, numOut)

        scs = [scores[mt] for scores, xf in fits]
        avg = numpy.average ( scs )
        stdev = numpy.std ( scs )
        print ""
        print "By %s" % mt
        print " - %d fits/scores, average %.2f, stdev %.4f, top:" % (len(fits), avg, stdev)

        #fits = sorted(fits, key=lambda x: x[0]['z-score'], reverse=True)
        fits = sorted(fits, key=lambda x: x[0][mt], reverse=True)

        for scores, xf in fits[0:10] :
            z = (scores[mt]-avg)/stdev
            print "%s:%.3f, z:%.3f, local-z-score %.3f" % (mt, scores[mt], z, scores['z-score'])

        scores, xf = fits[0]
        mol.openState.xform = xf

        ati = 1
        for scores, xf in fits[0:numOut] :
            z = (scores[mt] - avg)/stdev
            locz = scores['z-score']
            xfc = xf.__copy__()
            xfc.premultiply ( dmap.openState.xform.inverse() )

            nmol = CopyMolX ( mol, xfc )

            #nmol.name = "lz_%07.3f__z_%07.3f__%s_%.6f__%02d__%s.pdb" % ( locz, z, mt, scores[mt], ati, molName)
            nmol.name = "z_%07.3f__%s_%.6f__cc_%.6f__%02d__numRes_%d__numAt_%d__%s.pdb" % ( z, mt, scores[mt], scores['cc'], ati, len(mol.residues), len(mol.atoms), molName )
            #nmol.name = "%s_%03.6f__ccm_%.6f__cc_%.6f__%02d__%s.pdb" % ( mt, scores[mt], scores['ccm'], scores['cc'], ati, molName)
            ati += 1

            if numOut < 10 :
                mfile = os.path.join(fitsPath_, nmol.name )
                print " - saving mol %s" % mfile
                chimera.PDBio().writePDBfile ( [nmol], mfile )
            else :
                fitsPathN = fitsPathN_ + "/%s" % (molName)
                print " -- saving N %d to path: %s" % (numOut, fitsPathN)
                if checkPath (fitsPathN_) == False : return
                if checkPath (fitsPathN) == False : return
                mfile = os.path.join(fitsPathN, nmol.name )
                print " - saving mol %s" % mfile
                chimera.PDBio().writePDBfile ( [nmol], mfile )


    totSec = time.time() - start0
    totMin = numpy.floor ( totSec / 60.0 )
    totSec = totSec - totMin * 60.0
    print ""
    print " - done, took: %.0f min, %.0f sec" % ( totMin, totSec )


    # ------------------



def fitMPm ( dmap, toMap, numProc=None, toMapThr=None, dmapThr=None, step=8.0, useGrid=False, nrot=8, metric="ccm", numOut=10 ) :

    for m in chimera.openModels.list() :
        m.openState.xform = chimera.Xform.identity()

    print "Fit MP"
    print " - to map: %s, thr %.2f" % (toMap.name, toMap.surface_levels[0])
    print " - map: %s, thr %.2f" % ( dmap.name, dmap.surface_levels[0] )
    print " - step: %.2f" % step
    print " - num rotations: %.2f x %.2f = %.2f, angle %.2f" % (nrot, nrot, nrot*nrot, 360/float(nrot))

    if toMapThr != None :
        print " - %s threshold: %.2f" % (toMap.name, toMapThr)
    else :
        M = toMap.data.full_matrix()
        toMapThr = numpy.average(M)+numpy.std(M)*2.0
        print " - %s threshold at 2-sigma: %.2f" % (toMap.name, toMapThr)
        chimera.runCommand ( "vol #%d level %f" % (toMap.id, toMapThr) )

    if dmapThr != None :
        print " - %s threshold: %.2f" % (dmap.name, thrMap)
    else :
        M = dmap.data.full_matrix()
        dmapThr = numpy.average(M)+numpy.std(M)*2.0
        print " - %s threshold at 2-sigma: %.2f" % (dmap.name, dmapThr)
        chimera.runCommand ( "vol #%d level %f" % (dmap.id, dmapThr) )

    chimeraPath = GetChimeraPath ()
    print " - chimera path: %s" % chimeraPath

    dir_path = os.path.dirname(os.path.realpath(__file__))
    inDir = os.path.split(dir_path)[0]
    print " -- working dir:", inDir
    #mapQPPath = os.path.join ( inDir, 'Segger' )
    fitScriptPath = os.path.join ( dir_path, 'fitMPm.py' )
    print " -- path to fitScriptPath script:", fitScriptPath

    if numProc == None :
        import multiprocessing
        numProc = multiprocessing.cpu_count() / 2
        print " - auto num processors: %d" % numProc
        #numProc = 2
    else :
        print " - num processors: %d" % numProc

    #weight0, weight1 = 0.0, 0.0
    dmapPath = os.path.splitext(dmap.data.path)[0]
    dmapName = os.path.splitext(dmap.name)[0]
    mapName = os.path.splitext(toMap.name)[0]
    mapPath = os.path.split(toMap.data.path)[0]
    mapBase = os.path.splitext(toMap.data.path)[0]

    def checkPath ( path ) :
        if not os.path.isdir ( path ) :
            try :
                os.mkdir(path)
            except :
                print " - could not make: %s" % path
                return

    fitsPath = mapBase + "__fits_step%.0f_nrot%d" % (step, nrot)
    print " - fits path: %s" % fitsPath
    if checkPath (fitsPath) == False : return

    allFitsPath = fitsPath + "/all_fits"
    print " - all fits path: %s" % allFitsPath
    if checkPath (allFitsPath) == False : return

    # make path for temporary results
    from random import choice
    from string import ascii_uppercase, digits
    randStr = ''.join(choice(ascii_uppercase + digits) for _ in range(7))
    #tempPath = mapBase + "__ez_fit_temp__%s__" % randStr
    tempPath = mapBase + "__fitMPm__step%.0f_nrot%d__%s" % (step, nrot, dmapName)

    start0 = time.time()

    procs = []

    if os.path.isdir ( tempPath ) :
        print " - temp path exists, removing: %s" % tempPath
        from shutil import rmtree
        rmtree( tempPath )

    print " - making temp path: %s" % tempPath

    try :
        os.mkdir(tempPath)
    except :
        print " - could not make temp path"
        print "    : %s" % tempPath
        print "    - check/remove temp path manually and try again"
        print "    - or, check write permissions"
        return

    #molPoints, molWeights, molCenter, molMap = GetMolPts ( mol, metric, resolution )
    #print " - fitting %d points with metric %s" % (len(molPoints), metric)

    mapPts = None
    if useGrid :
        mapPts = GetMapPtsGrid ( toMap, thr=toMapThr, step=step )
        #mapPts = mapPts [ len(mapPts)/2 : len(mapPts)/2 + 2 ]
    else :
        mapPts = GetMapPts ( toMap, toMapThr, step, numMax=None )
        #mapPts = mapPts[0:18]

    #fits = fitPts ( mapPts, dmap, molPoints, molWeights, molCenter, nrot, metric )
    # ------------------

    n = len(mapPts)
    g = [mapPts[(n*c)/numProc:(n*(c+1))/numProc] for c in range(numProc)]
    for mi, mpts in enumerate(g) :

        #nmap_path = os.path.join ( tempPath, "%d_map.mrc" % mi )
        dmapMiPath = os.path.join ( tempPath, dmapName + "__%d.mrc" % mi )
        shutil.copyfile ( dmap.data.path, dmapMiPath )

        #mpts = mpts[0:2]
        print " - proc %d/%d, %d map points" % (mi+1, numProc, len(mpts))

        procPointsPath = os.path.join ( tempPath, "%d_points.txt" % mi )
        fout = open ( procPointsPath, "w" )
        #fout.write ( "%s\n" % dmap.data.path )
        fout.write ( "%.2f\n" % toMapThr )
        fout.write ( "%.2f\n" % dmapThr )
        fout.write ( "%.2f\n" % step )
        fout.write ( "%s\n" % metric )
        fout.write ( "%d\n" % nrot )
        for mpt in mpts :
            fout.write ( "%f %f %f\n" % (mpt[0], mpt[1], mpt[2]) )
        fout.close()

        args = [chimeraPath, '--nogui', '--silent', '--nostatus', toMap.data.path, dmapMiPath, fitScriptPath]
        if mi == 0 :
            print "  >" + " ".join(args)

        fout = open ( os.path.join(tempPath, "%d_out.txt" % mi), "w" )
        foute = open ( os.path.join(tempPath, "%d_err.txt" % mi), "w" )

        inPath = "/Users/greg/Dropbox/_mol"
        if os.path.exists ( "/storage/toxop" ) : inPath="/storage/toxop"
        elif os.path.exists ( "/scratch/gregp/toxop" ) : inPath="/scratch/gregp/toxop"
        elif os.path.exists ( "/sdf/home/g/gregp/g24" ) : inPath="/sdf/home/g/gregp/g24"
        elif os.path.exists ( "/storage/gregp/toxop" ) : inPath="/storage/gregp/toxop"
        elif os.path.exists ( "/sdf/home/g/gregp/g24/ez" ) : inPath="/sdf/home/g/gregp/g24/ez"
        inPath = "./"
        print " - inPath: %s" % inPath

        from subprocess import Popen
        import subprocess
        p = Popen(args, stdout=fout, stderr=foute, cwd=inPath)
        procs.append ( [mi, p, fout, foute] )

        #p = Popen(args, stdout=subprocess.PIPE, stderr=foute, cwd=inPath)
        #procs.append ( [mi, p, foute] )
        #p.done = False


    print ""
    print "Waiting for fits..."
    print ""
    posMap = None
    totFits = 0
    error = False
    for mi, p, fout, foute in procs :
        p.wait()
        fout.close()
        foute.close()

        fitsf = os.path.join(tempPath, "%d_fits.txt" % mi)
        print "%d -- %s" % (mi, fitsf)

        if not os.path.isfile ( fitsf ) :
            print " -x- fits file not found -x-"
            error = True
            continue

        ff = open ( fitsf )
        fiti = 0
        for l in ff :
            ll = l.strip()
            t = ll.split()
            avgd, cc, ccm, olap, z, ang = float(t[0]), float(t[1]), float(t[2]), float(t[3]), float(t[4]), float(t[5])
            scores = {"avgd":avgd, "cc":cc, "ccm":ccm, "olap":olap, "z-score":z}
            v = chimera.Vector ( float(t[6]), float(t[7]), float(t[8]) ); v.normalize()
            t = chimera.Vector ( float(t[9]), float(t[10]), float(t[11]) )
            xf = chimera.Xform.translation ( t )
            xf.multiply ( chimera.Xform.rotation ( v, ang ) )
            xf.premultiply ( dmap.openState.xform )
            #addFit ( fits, xf, scores )
            posMap = addFit__ ( posMap, xf, scores )
            totFits += 1
            fiti += 1
        ff.close()
        print "%d -- %d fits" % (mi, fiti)

    print ""

    if posMap == None :
        print " -- no fits??? --"
        return

    if 0 :
        print " - removing temp path: %s" % tempPath
        from shutil import rmtree
        rmtree ( tempPath )

    if error :
        print " -- stopping with error -- "
        return

    fits = []
    for i, posMapJ in posMap.iteritems() :
        for j, posMapK in posMapJ.iteritems() :
            for k, posMapXfs in posMapK.iteritems() :
                fits.extend ( posMapXfs )

    print ""
    print " - %d unique fits from %d total" % (len(fits), totFits)

    #fitsPath = mapBase + "__fits_res%.2f_step%.0f_nrot%d_bf%.0f" % (resolution, step, nrot, bCutoff)
    print " - fits path: %s" % fitsPath
    if checkPath (fitsPath) == False : return

    allFitsPath = fitsPath + "/all_fits"
    print " - all fits path: %s" % allFitsPath
    if checkPath (allFitsPath) == False : return

    fp = open ( allFitsPath + "/" + dmapName + "_all_fits.txt", "w" )
    fp.write ( "Avg.Density\tZ(Avg.D.)\tCC\tZ(CC)\tCC-mean\tZ(CC-mean)\tAngle\tRotAx_X\tRotAx_Y\tRotAx_Z\tTrans_X\tTrans_Y\tTrans_Z\n" )

    if len(fits) < 1 :
        print " - no fits found... try a lower threshold"
        fp.close()
        return

    avg, std = {}, {}
    for tp in ['avgd', 'cc', 'ccm'] :
        scs = [scores[tp] for scores, xf in fits]
        avg[tp] = numpy.average ( scs )
        std[tp] = numpy.std ( scs )

    print " - sorting fits..."
    fits = sorted(fits, key=lambda x: x[0]['ccm'], reverse=True)

    for scores, xf in fits :
        v, ang = xf.getRotation()
        t = xf.getTranslation()
        avgd, cc, ccm = scores['avgd'], scores['cc'], scores['ccm']
        # scores['olap'], scores['z-score']
        avgd_z = (avgd - avg['avgd']) / std['avgd']
        cc_z =   (cc -   avg['cc'])   / std['cc']
        ccm_z =  (ccm -  avg['ccm'])  / std['ccm']
        fp.write ( "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (avgd, avgd_z, cc, cc_z, ccm, ccm_z , ang, v[0], v[1], v[2], t[0], t[1], t[2]) )

        scores['avgd_z'] = avgd_z
        scores['cc_z'] = cc_z
        scores['ccm_z'] = ccm_z

    fp.close()


    #for mt in ["avgd", "cc", "ccm"] :
    #for mt in [ "cc", "ccm", "olap" ] :
    #for mt in [ "ccm", "cc", 'z-score' ] :
    for mt in [ "ccm", "cc" ] :

        fitsPath_ = fitsPath + "/by_%s" % mt
        if checkPath (fitsPath_) == False : return
        fitsPathN_ = fitsPath + "/by_%s_N_%d" % (mt, numOut)

        scs = [scores[mt] for scores, xf in fits]
        avg = numpy.average ( scs )
        stdev = numpy.std ( scs )
        print ""
        print "By %s" % mt
        print " - %d fits/scores, average %.2f, stdev %.4f, top:" % (len(fits), avg, stdev)

        #fits = sorted(fits, key=lambda x: x[0]['z-score'], reverse=True)
        fits = sorted(fits, key=lambda x: x[0][mt], reverse=True)

        for scores, xf in fits[0:10] :
            z = (scores[mt]-avg)/stdev
            print "%s:%.3f, z:%.3f, local-z-score %.3f" % (mt, scores[mt], z, scores['z-score'])

        scores, xf = fits[0]

        ati = 1
        for scores, xf in fits[0:numOut] :
            z = (scores[mt] - avg)/stdev
            locz = scores['z-score']
            #xfc = xf.__copy__()
            #xfc.premultiply ( toMap.openState.xform.inverse() )
            dmap.openState.xform = xf

            #nmol = CopyMolX ( mol, xfc )
            dname = "z_%07.3f__%s_%.6f__cc_%.6f__%02d__%s.mrc" % ( z, mt, scores[mt], scores['cc'], ati, dmapName )
            ati += 1

            dpath = os.path.join(fitsPath_, dname )
            print " - saving fitted map %s" % dpath
            #chimera.PDBio().writePDBfile ( [nmol], mfile )
            #nv = MapToMap0 ( fitPoints, toMap, dmap, "fit" )
            MapToMap ( dmap, dmapThr, toMap, nname=dname, saveToPath=dpath )


    totSec = time.time() - start0
    totMin = numpy.floor ( totSec / 60.0 )
    totSec = totSec - totMin * 60.0
    print ""
    print " - done, took: %.0f min, %.0f sec" % ( totMin, totSec )


    # ------------------







def fitAll ( numProc=6, metric="ccm", resolution=7.0, nrot=8, step=10, numRandFit=None, numOut=5 ) :

    maps = chimera.openModels.list (modelTypes = [VolumeViewer.volume.Volume])
    mols = chimera.openModels.list (modelTypes = [chimera.Molecule])

    print ""
    print "%d maps, %d mols" % (len(maps), len(mols))

    for dmap in maps :
        if dmap.display == True :
            print ""
            print "Map: %s" % dmap.name
            for mol in mols :
                if mol.display == True :
                    if numProc == 1 :
                        fit ( dmap, mol, metric=metric, resolution=resolution, step=step, nrot=nrot, numRandFit=numRandFit, numOut=4 )
                    else :
                        fitMP ( dmap, mol, numProc=numProc, metric=metric, resolution=resolution, step=step, nrot=nrot, bCutoff=80.0, numOut=numOut )




def fitAll2 ( numProc=6, metric="ccm", resolution=7.0, nrot=8, step=10, numRandFit=None, numOut=5, bCutoff=None ) :

    m0 = chimera.openModels.list() [0]
    m1 = chimera.openModels.list() [1]

    if numProc == 1 :
        fit ( m0, m1, metric=metric, resolution=resolution, step=step, nrot=nrot, numRandFit=numRandFit, numOut=4 )
    else :
        fitMP ( m0, m1, numProc=numProc, metric=metric, resolution=resolution, step=step, nrot=nrot, bCutoff=bCutoff, numOut=numOut )



def fitScores ( resolution ) :

    maps = chimera.openModels.list (modelTypes = [VolumeViewer.volume.Volume])
    mols = chimera.openModels.list (modelTypes = [chimera.Molecule])

    print ""
    print "%d maps, %d mols" % (len(maps), len(mols))

    for dmap in maps :
        if dmap.display == True :
            print ""
            print "Map: %s" % dmap.name
            for mol in mols :
                if mol.display == True :
                    #fit ( dmap, mol )

                    print "Mol: %s" % mol.name
                    points, weights, molCenter, molMap = GetMolPts ( mol, "cc", resolution )
                    print " - %d points at res %.2f, step %.2f" % (len(points), resolution, molMap.data.step[0])

                    if 1 and molMap != None :
                        chimera.openModels.close( [molMap] )

                    #fits = fitPts ( mapPts, dmap, molPoints, molWeights, molCenter, nrot, metric )

                    dmapXf = dmap.openState.xform # .inverse()
                    dmapXfi = dmapXf.inverse()

                    #xf = mol.openState.xform.__copy__()
                    #xf.premultiply ( dmapXfi )

                    dmapXfi.multiply ( mol.openState.xform )

                    xfm = Matrix.xform_matrix ( dmapXfi )
                    #xyz_to_ijk_tf = dmap.data.xyz_to_ijk_transform
                    xyz_to_ijk_tf = Matrix.multiply_matrices ( dmap.data.xyz_to_ijk_transform, xfm )
                    darray = dmap.data.full_matrix()

                    map_values, outside = VolumeData.interpolate_volume_data(points, xyz_to_ijk_tf, darray)
                    olap, cc, ccm = FitMap.overlap_and_correlation ( weights, map_values )
                    #print " - score ?", score, cc, ccm
                    print "overlap: %.4f" % olap
                    print "CC : %.4f" % cc
                    print "CCm: %.4f" % ccm




#
