

import chimera
import _surface
import numpy


from chimera import Vector


def Quad2Tri ( vi ) :
    t1 = (vi[0], vi[1], vi[2])
    t2 = (vi[0], vi[2], vi[3])
    return t1, t2


def SubdivideQuadRec ( pts, qt, numit ) :

    for i in range (numit) :
        nq = ()
        for q in qt :
            pts, q2 = SubdivideQuad (pts, q)
            if q2 :
                nq = nq + q2
            else :
                nq = nq + q

        #print i, len(qt), "*** Old quads" #, qt
        #print i, len(nq), "*** New quads" #, nq
        if len(nq) == len(qt) :
            break
        qt = nq
        #break

    return pts, qt


def SubdivideQuad ( pts, q ) :
    # 4 pts -
    if ( len(q) != 4 ) :
        print "SubdivideQuad needs a 4-tuple, got", q
        return pts

    #print "Points:\n", pts
    #print "Quad:\n", q

    v1 = pts[q[1]] - pts[q[0]]
    v2 = pts[q[2]] - pts[q[1]]
    l1 = numpy.sqrt(numpy.dot(v1,v1))
    l2 = numpy.sqrt(numpy.dot(v2,v2))
    # print "V1: ", v1, l1
    # print "V2: ", v2, l2

    q2 = None

    if ( l2 > 1 and l2 > l1 ) :
        m1 = pts[q[0]] + v2 * numpy.array([.5],'f')
        m2 = pts[q[1]] + v2 * numpy.array([.5],'f')
        li = len(pts)
        pts = numpy.concatenate( [pts, [m1]] )
        pts = numpy.concatenate( [pts, [m2]] )
        q2 = ( (q[0], q[1], li+1, li), (li, li+1, q[2], q[3]) )
    elif ( l1 > 1 ) :
        m1 = pts[q[0]] + v1 * numpy.array([.5],'f')
        m2 = pts[q[3]] + v1 * numpy.array([.5],'f')
        li = len(pts)
        pts = numpy.concatenate( [pts, [m1]] )
        pts = numpy.concatenate( [pts, [m2]] )
        q2 = ( (q[0], li, li+1, q[3]), (li, q[1], q[2], li+1) )

    return pts, q2



def AddWalls () :

    m = _surface.Surface_Model()
    v = numpy.array( (
        (-1,-1,-1), (-1,1,-1), (1,1,-1), (1,-1,-1),
        (-1,-1,1), (-1,1,1), (1,1,1), (1,-1,1)
        ), numpy.float32 )
    #for i in range ( len(v) ) :
    #    v[i] = v[i] * .1
    v = v * WD
    v = v.astype(numpy.float32)

    vi_floor = (
        (0,1,2), (0,2,3)
        )

    vi_walls = (
        (0,4,1), (1,4,5), (1,5,2), (2,5,6)
        )

    #vi_walls_mesh = (
    #    (2,6,3), (6,7,3), (3,4,0), (3,7,4)
    #    )

    # q = (2,3,7,6)
    # vi_walls_mesh = Quad2Tri ( q )
    # v, q2 = SubdivideQuad (v, q)

    qt = ( (2,6,7,3), )
    v_wall_1, qt = SubdivideQuadRec ( v, qt, 5 )
    vi_wall_1_mesh = ()
    for q in qt :
        vi_wall_1_mesh = vi_wall_1_mesh + Quad2Tri (q)

    qt = ( (3,7,4,0), )
    v_wall_2, qt = SubdivideQuadRec ( v, qt, 5 )
    vi_wall_2_mesh = ()
    for q in qt :
        vi_wall_2_mesh = vi_wall_2_mesh + Quad2Tri (q)

    red = (1,0,0,1)
    grey = (.7,.7,.7,1)

    g_floor = m.add_group(v, vi_floor, red)
    g_walls = m.add_group(v, vi_walls, red)
    g_wall_1_mesh = m.add_group(v_wall_1, vi_wall_1_mesh, grey)
    g_wall_2_mesh = m.add_group(v_wall_2, vi_wall_2_mesh, grey)

    g_wall_1_mesh.set_display_style(g_wall_1_mesh.Mesh)
    g_wall_2_mesh.set_display_style(g_wall_2_mesh.Mesh)

    chimera.openModels.add([m])

    return m




def AddWireWalls ( dim, ctr ) :

    dim[0] = dim[0]/2.0; dim[1] = dim[1]/2.0; dim[2] = dim[2]/2.0

    m = _surface.Surface_Model()
    v = numpy.array( (
        (ctr[0]-dim[0],ctr[1]-dim[1],ctr[2]-dim[2]),
        (ctr[0]-dim[0],ctr[1]+dim[1],ctr[2]-dim[2]),
        (ctr[0]+dim[0],ctr[1]+dim[1],ctr[2]-dim[2]),
        (ctr[0]+dim[0],ctr[1]-dim[1],ctr[2]-dim[2]),
        (ctr[0]-dim[0],ctr[1]-dim[1],ctr[2]+dim[2]),
        (ctr[0]-dim[0],ctr[1]+dim[1],ctr[2]+dim[2]),
        (ctr[0]+dim[0],ctr[1]+dim[1],ctr[2]+dim[2]),
        (ctr[0]+dim[0],ctr[1]-dim[1],ctr[2]+dim[2])
        ), numpy.float32 )

    # v = v * lengthX2
    # v = v.astype(numpy.float32)

    qt = ( (2,6,7,3), )
    qt = ( (3,7,4,0), )

    vi_floor = ((0,1,2), (0,2,3))
    vi_walls = ((0,4,1), (1,4,5), (1,5,2), (2,5,6),
                (2,6,7), (2,7,3), (3,7,4), (3,4,0) )

    red = (1,0,0,0.1)
    grey = (.2,.2,.2,0.0)

    g_floor = m.add_group(v, vi_floor, grey)
    g_walls = m.add_group(v, vi_walls, grey)

    g_floor.set_display_style(g_floor.Mesh)
    g_walls.set_display_style(g_walls.Mesh)

    chimera.openModels.add([m])

    return m




def SphereMesh (r, div, color, patchpts, pos = Vector(0,0,0)) :

    m = _surface.Surface_Model()

    v = numpy.array( [ [0+pos.x,0+pos.y,r+pos.z], ], numpy.float32 )
    vi = ()

    at = 1
    l = int ( numpy.ceil (float(div)*3.0/2.0) )
    if div < 10 : l = div*2
    print "SphereMesh:", div, 'x', l
    lat = 0

    for phi_i in range(div) :

        phi = 90.0 - ( float(phi_i+1) * 180.0/float(div+1) )
        #print "%.2f: " % phi,
        z = r * numpy.sin(phi * numpy.pi/180)
        s = r * numpy.cos(phi * numpy.pi/180)

        for psi_i in range (l) :
            psi = float(psi_i) * 360.0/float(l)

            #print "%.0f(%d)(%d)" % (psi, at, at-l),
            x = s * numpy.sin(psi * numpy.pi/180)
            y = s * numpy.cos(psi * numpy.pi/180)

            pt = numpy.array( [ [x+pos.x,y+pos.y,z+pos.z], ], numpy.float32 )
            v = numpy.concatenate ( [v, pt] )

            if phi_i == 0 :
                if psi_i > 0 :
                    vi = vi + ( (at-1, at, 0), )
                if psi_i == l-1 :
                    vi = vi + ( (at, 1, 0), )
            else :
                if psi_i > 0 :
                    tris = Quad2Tri ( [at-1, at, at-l, at-l-1] )
                    vi = vi + tris
                if psi_i == l-1 :
                    tris = Quad2Tri ( [at, at-l+1, at-l*2+1, at-l] )
                    vi = vi + tris

            if phi_i == div-1 :
                if psi_i > 0 :
                    vi = vi + ( (at, at-1, lat+l), )
                if psi_i == l-1 :
                    vi = vi + ( (at-l+1, at, lat+l), )

            at = at + 1


        lat = len ( v )

    pt = numpy.array( [ [0+pos.x,0+pos.y,-r+pos.z], ], numpy.float32 )
    v = numpy.concatenate ( [v, pt] )


    sph = m.add_group( v, vi, color )
    #sph.set_display_style(sph.Mesh)


    if patchpts :
        vcolors = ()

        for i in range ( len(v) ) :
            vp = chimera.Vector (  v[i][0], v[i][1], v[i][2] ) - pos
            inP = None
            for pt in patchpts:
                if (pt[0] - vp).length < (r/3) :
                    inP = pt[1]
                    break
            if inP :
                if inP < 0.0 :
                    vcolors = vcolors + ( (-inP*.6+.4, .4, .4, 1), )
                else :
                    vcolors = vcolors + ( (.4, .4, inP*.6+.4, 1), )
            else :
                vcolors = vcolors + ( (color), )

        print len ( vcolors ), len ( v )

        sph.set_vertex_colors( vcolors )


    chimera.openModels.add([m])

    return m




def CylinderMesh (r1, r2, Length, div, color) :

    m = _surface.Surface_Model()
    chimera.openModels.add([m])

    v = None
    vi = ()

    # print "CylinderMesh:", div

    at = 0
    for psi_i in range(div) :

        psi = float(psi_i) * 360.0/float(div)

        #print "%.0f(%d)(%d)" % (psi, at, at-l),
        x1 = r1 * numpy.sin(psi * numpy.pi/180)
        y1 = r1 * numpy.cos(psi * numpy.pi/180)

        x2 = r2 * numpy.sin(psi * numpy.pi/180)
        y2 = r2 * numpy.cos(psi * numpy.pi/180)

        if psi_i == 0 :
            v = numpy.array( [ [x1,y1,0], ], numpy.float32 )
        else :
            pt1 = numpy.array( [ [x1,y1,0], ], numpy.float32 )
            v = numpy.concatenate ( [v, pt1] )

        pt2 = numpy.array( [ [x2,y2,Length], ], numpy.float32 )
        v = numpy.concatenate ( [v, pt2] )

        at = at + 2

        if psi_i == 0 :
            pass
        else :
            tris = Quad2Tri ( [at-4, at-2, at-1, at-3] )
            vi = vi + tris

        if psi_i == div-1 :
            tris = Quad2Tri ( [at-2, 0, 1, at-1] )
            vi = vi + tris


    pt1 = numpy.array( [ [0,0,0], ], numpy.float32 )
    v = numpy.concatenate ( [v, pt1] )

    pt1 = numpy.array( [ [0,0,Length], ], numpy.float32 )
    v = numpy.concatenate ( [v, pt1] )

    if 0 and r1 > .01 :
        print "capping 1"
        vi = vi + ( (at, 0, at-2), )
        for i in range ( (at-2)/2 ) :
            vi = vi + ( (at, (i+1)*2, (i+0)*2), )

    if 0 and r2 > .01 :
        print "capping 2"
        vi = vi + ( (at+1, at-1, 1), )
        for i in range ( (at-2)/2 ) :
            vi = vi + ( (at+1, (i+0)*2+1, (i+1)*2+1), )


    sph = m.add_group( v, vi, color )
    return m


# execfile("c:\greg\chimera\Blob\gui.py")

def ReadMesh (patchpts = None, pos=Vector(0,0,0)) :

    m = _surface.Surface_Model()

    com = Vector (0,0,0)
    rad = Vector (0,0,0)
    numv = 0

    aV = []
    fp = open ( "C:\\greg\\chimera\\Blob\\psu_points.txt", 'r' )
    for line in fp :
        n = line.split(',')
        for i in range( len(n) ) :
            c = n[i].split()
            if len(c) == 3 :
                #v = numpy.array([[float(c[0]),float(c[1]),float(c[2])]], 'f' )
                aV = aV + [ [float(c[0]), float(c[1]), float(c[2])], ]
                numv = numv + 1
                v = Vector ( float(c[0]), float(c[1]), float(c[2]) )
                com = com + v
                if v.length > rad.length :
                    rad = v
    fp.close()

    v = numpy.array( aV, 'f' )
    print "COM:", com/float(numv)
    print "Rad:", rad.length
    scale = 1.0/rad.length

    xf = chimera.Xform.rotation ( Vector(0,0,1), 180 )
    xf.multiply ( chimera.Xform.rotation ( Vector(1,0,0), -90 ) )
    for i in range (numv) :
        vec = xf.apply ( Vector ( v[i][0], v[i][1], v[i][2] ) * scale )
        v[i][0] = vec.x
        v[i][1] = vec.y
        v[i][2] = vec.z


    vi = []
    vs = []
    fp = open ( "C:\\greg\\chimera\\Blob\\psu_tris.txt", 'r' )
    for line in fp :
        n = line.split(',')
        for t in n :
            try :
                ivi = int(t)
                if ivi == -1 :
                    #print 'tri', vs
                    vi = vi + [vs]
                    vs = []
                else :
                    vs = vs + [ivi]
            except:
                #print "bad token:", t
                continue
    fp.close()

    color = (0.6039, 0.8431, 0.898, 1.0)
    sph = m.add_group( v, vi, color )
    #sph.set_display_style(sph.Mesh)


    if patchpts :
        vcolors = ()

        for i in range ( len(v) ) :
            vp = chimera.Vector (  v[i][0], v[i][1], v[i][2] ) - pos
            inP = None
            for pt in patchpts:
                if (pt[0] - vp).length < (1.0/2.5) :
                    inP = pt[1]
                    break
            if inP :
                if inP < 0.0 :
                    vcolors = vcolors + ( (-inP*.6+.4, .4, .4, 1), )
                else :
                    vcolors = vcolors + ( (.4, .4, inP*.6+.4, 1), )
            else :
                vcolors = vcolors + ( (color), )

        print len ( vcolors ), len ( v )

        sph.set_vertex_colors( vcolors )


    chimera.openModels.add([m])
    return m



def ReadMesh2 (fname, m) :

    if m == None :
        m = _surface.SurfaceModel()

    com = Vector (0,0,0)
    rad = Vector (0,0,0)
    numv = 0

    aV = []
    fp = open ( fname, 'r' )

    l1 = fp.readline()
    numV, numT = l1.split ()
    numV = int ( numV )
    numT = int ( numT )

    print "%d verts, %d tris" % (numV, numT)

    verts = numpy.ones ( [numV, 3] )
    tris = numpy.ones ( [numT, 3], numpy.int )


    for vi in range ( numV ) :
        line = fp.readline()
        c = line.split(' ')
        if len(c) == 3 :
            verts[vi] = c

    for ti in range ( numT ) :
        line = fp.readline()
        n = line.split(' ')
        tris[ti] = n

    fp.close()

    color = (0.6039, 0.8431, 0.898, 1.0)
    patch = m.addPiece( verts, tris, color )
    #sph.set_display_style(sph.Mesh)
    #sph.set_vertex_colors( vcolors )

    return m



def MeshFromVertsTris (verts, tris, color=None, m=None) :

    if m == None :
        m = _surface.SurfaceModel()

    com = Vector (0,0,0)
    rad = Vector (0,0,0)
    numv = 0

    print " - mesh from %d verts, %d tris" % (len(verts), len(tris))

    if color == None :
        color = (0.6039, 0.8431, 0.898, 1.0)

    patch = m.addPiece( verts, tris, color )
    #sph.set_display_style(sph.Mesh)
    #sph.set_vertex_colors( vcolors )

    return m





def MakeColors ( n ) :

    # create n colors by interpolating between 5 basic colors:

    clrs = numpy.array ( [[.9,.3,.3], [.9,.9,.3], [.3,.9,.3],
                          [.3,.9,.9], [.3,.3,.9]] )

    if n == 1 :
        return [ clrs[0] ]

    at = 0.000001
    d = -0.00001 + float ( len(clrs) - 1  ) / float ( n - 1  )

    colors = []

    for i in range (n) :
        l = int ( numpy.floor ( at ) )
        u = int ( numpy.ceil ( at ) )
        print "color %d (%f) [%d-%d]" % (i, at, l, u),
        clr = clrs[l] * (at-l) + clrs[u] * (u-at)
        print " %f %f %f" % (clr[0], clr[1], clr[2])
        colors.append ( clr )
        at = at + d

    return colors



def AddDiffRandColor ( clist, on_black = True, tol = .5 ) :

    R,G,B = None, None, None
    random.seed()

    for tot in range(1000) :

        if on_black :
            R = random.random()*.6+.4
            G = random.random()*.6+.4
            B = random.random()*.6+.4
        else :
            R = random.random()*.7
            G = random.random()*.7
            B = random.random()*.7

        tooclose = False
        for c in clist :
            r,g,b = c[0], c[1], c[2];
            if abs(r-R) + abs(g-G) + abs(b-B) < tol :
                tooclose = True

        if not tooclose: break

    clist.append ( numpy.array ( [R,G,B] ) )
    print " - New rand color: %.3f %.3f %.3f"%(R,G,B)

    return clist

def MakeDiffRandColors ( n ) :

    clrs = []
    for i in range (n) :
        AddDiffRandColor ( clrs )

    return clrs


def rrange ( n ) :

    i = range ( n )
    for n in i :
        ri = int ( numpy.round ( random.random() * n ) )
        # print "swapping %d with %d" % (n, ri)
        t = i[n]
        i[n] = i[ri]
        i[ri] = t

    return i





def AlignXf ( pos, v ) :
    Z = v
    Z.normalize()
    dZ = Vector( random.random(), random.random(), random.random() )
    dZ.normalize()
    X = chimera.cross ( Z, dZ )
    X.normalize ()
    Y = chimera.cross ( Z, X )
    Y.normalize ()

    xf = chimera.Xform.xform (
        X.x, Y.x, Z.x, pos.x,
        X.y, Y.y, Z.y, pos.y,
        X.z, Y.z, Z.z, pos.z )

    #xf3 = chimera.Xform.xform (
    #    d, 0, 0, 0,
    #    0, d, 0, 0,
    #    0, 0, d, 0 )
    #print xf3

    return xf


def AddArrow ( pos, v, d, mxf, clr=(0,1,1,1), rad=0.2 ) :

    #d = v.length
    # mxf = a[0].molecule.openState.xform

    xf = AlignXf ( pos, v )
    if mxf != None : xf.premultiply ( mxf )
    a = CylinderMesh ( rad, rad, d-(rad*3), 10, clr )
    a.openState.xform = xf

    xf = AlignXf ( pos+(v*(d-(rad*3))), v )
    if mxf != None : xf.premultiply ( mxf )
    b = CylinderMesh ( rad*3, 0.01, rad*3, 10, clr )
    b.openState.xform = xf

    return [a, b]



class MyBBox :


    def __init__ (self) :

        self.llf = chimera.Point (1e99, 1e99, 1e99)
        self.urb = chimera.Point (-1e99, -1e99, -1e99)
        self.center = chimera.Point (0, 0, 0)
        self.dim = chimera.Point (0, 0, 0)


    def add (self, pt) :

        if pt[0] < self.llf[0] : self.llf[0] = pt[0]
        if pt[1] < self.llf[1] : self.llf[1] = pt[1]
        if pt[2] < self.llf[2] : self.llf[2] = pt[2]

        if pt[0] > self.urb[0] : self.urb[0] = pt[0]
        if pt[1] > self.urb[1] : self.urb[1] = pt[1]
        if pt[2] > self.urb[2] : self.urb[2] = pt[2]

        self.center[0] = (self.urb[0] + self.llf[0]) / 2.0
        self.center[1] = (self.urb[1] + self.llf[1]) / 2.0
        self.center[2] = (self.urb[2] + self.llf[2]) / 2.0

        self.dim = self.urb - self.llf


    def AddModel (self) :

        self.model = AddWireWalls ( self.dim, self.center )

        #T = chimera.Xform.translation ( self.center.toVector() )
        #S = chimera.Xform.translation ( 1,2,3 )
        #self.model.openState.xform = xf.translation(1,2,3)




def Measure ( mol ) :

    B = MyBBox ()

    for at in mol.atoms :
        B.add ( at.coord() )

    print "B-box: (%.2f %.2f %.2f) - (%.2f %.2f %.2f), center (%.2f %.2f %.2f), dim (%.2f %.2f %.2f)" % (
        B.llf.x, B.llf[1], B.llf[2], B.urb[0], B.urb[1], B.urb[2],
        B.center[0], B.center[1], B.center[2], B.dim[0], B.dim[1], B.dim[2] )
