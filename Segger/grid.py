

import chimera
import numpy
import time
import _contour
import _multiscale



class Grid (object) :

    boxes = None # array of boxes
    O = None # origin
    D = None # dimension of boxes
    tree = None # for comparison

    def __init__ (self) :
        self.boxes = None # array of boxes
        self.O = None # origin
        self.D = None # dimension of boxes
        self.tree = None # for comparison


    def FromMols (self, mols, maxD) :

        ats = []
        for m in mols :
            ats.extend ( m.atoms )

        self.FromAtoms ( ats, maxD )


    def FromAtoms ( self, atoms, maxD ) :

        minP, maxP = chimera.Point(1e9,1e9,1e9), chimera.Point(-1e9,-1e9,-1e9)

        if 1 :

            startt = time.time()

            for at in atoms :
                C = at.xformCoord()
                if C[0] < minP[0] : minP[0] = C[0]
                if C[1] < minP[1] : minP[1] = C[1]
                if C[2] < minP[2] : minP[2] = C[2]

                if C[0] > maxP[0] : maxP[0] = C[0]
                if C[1] > maxP[1] : maxP[1] = C[1]
                if C[2] > maxP[2] : maxP[2] = C[2]

            #print minP, " -> ", maxP

            minP = minP - chimera.Vector(maxD,maxD,maxD)
            maxP = maxP + chimera.Vector(maxD,maxD,maxD)

        else :

            points1 = numpy.copy ( points )
            _contour.affine_transform_vertices ( points1, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )
            points0 = numpy.copy ( points1 )
            _contour.affine_transform_vertices ( points1, dmap.data.xyz_to_ijk_transform )

            bound = maxD
            li,lj,lk = numpy.min ( points1, axis=0 ) - (bound, bound, bound)
            hi,hj,hk = numpy.max ( points1, axis=0 ) + (bound, bound, bound)

            minP, maxP = chimera.Point(li,lj,lk), chimera.Point(hi,hj,hk)


        n0 = int ( numpy.ceil ( (maxP[0] - minP[0]) / maxD ) )
        n1 = int ( numpy.ceil ( (maxP[1] - minP[1]) / maxD ) )
        n2 = int ( numpy.ceil ( (maxP[2] - minP[2]) / maxD ) )

        n = [ n0, n1, n2 ]
        N = n0 * n1 * n2

        if N > 1e8 :
            print " - grid max mem?"
            return

        dur = time.time() - startt

        print n, " -> ", N
        #print "0 -> %.3f sec" % dur


        if 0 :
            startt = time.time()
            boxes = [ [] ] * N
            dur = time.time() - startt
            print "alloc all -> %.3f sec" % dur

        startt = time.time()
        boxes = [ [] ] * n2
        for k in range(n2) :
            jboxes = [ [] ] * n1
            boxes[k] = jboxes
            for j in range(n1) :
                iboxes = [ [] ] * n0
                jboxes[j] = iboxes
                for i in range(n0) :
                    aboxes = []
                    iboxes[i] = aboxes

        dur = time.time() - startt
        #print "alloc [] -> %.3f sec" % dur


        startt = time.time()
        for at in atoms :
            C = at.xformCoord() - minP
            i = int ( numpy.floor ( C[0]/maxD ) )
            j = int ( numpy.floor ( C[1]/maxD ) )
            k = int ( numpy.floor ( C[2]/maxD ) )
            boxes[k][j][i].append  ( at )

        dur = time.time() - startt
        #print "put -> %.3f sec" % dur

        if 0 :
            ats = []
            for m in mols :
                ats.extend ( m.atoms )

            import _multiscale
            from CGLutil.AdaptiveTree import AdaptiveTree

            startt = time.time()
            points = _multiscale.get_atom_coordinates ( ats, transformed = True )
            self.tree = AdaptiveTree ( points.tolist(), ats, maxD)
            dur = time.time() - startt
            #print "tree -> %.3f sec" % dur

        self.boxes = boxes
        self.O = minP
        self.D = maxD
        self.n = n
        self.ni = n0
        self.nj = n1
        self.nk = n2




    def FromAtomsLocal (self, atoms, maxD) :

        points = _multiscale.get_atom_coordinates ( atoms, transformed = False )
        #_contour.affine_transform_vertices ( points1, Matrix.xform_matrix( dmap.openState.xform.inverse() ) )
        #points0 = numpy.copy ( points1 )

        li,lj,lk = numpy.min ( points, axis=0 ) - (maxD, maxD, maxD)
        hi,hj,hk = numpy.max ( points, axis=0 ) + (maxD, maxD, maxD)

        minP, maxP = chimera.Point(li,lj,lk), chimera.Point(hi,hj,hk)
        #print minP, " -> ", maxP

        n0 = int ( numpy.ceil ( (maxP[0] - minP[0]) / maxD ) )
        n1 = int ( numpy.ceil ( (maxP[1] - minP[1]) / maxD ) )
        n2 = int ( numpy.ceil ( (maxP[2] - minP[2]) / maxD ) )

        n = [ n0, n1, n2 ]
        N = n0 * n1 * n2

        if N > 1e8 :
            print " - grid max mem?"
            return

        #dur = time.time() - startt
        #print n, " -> ", N
        #print "0 -> %.3f sec" % dur


        if 0 :
            startt = time.time()
            boxes = [ [] ] * N
            dur = time.time() - startt
            print "alloc all -> %.3f sec" % dur

        startt = time.time()
        boxes = [ [] ] * n2
        for k in range(n2) :
            jboxes = [ [] ] * n1
            boxes[k] = jboxes
            for j in range(n1) :
                iboxes = [ [] ] * n0
                jboxes[j] = iboxes
                for i in range(n0) :
                    aboxes = []
                    iboxes[i] = aboxes

        dur = time.time() - startt
        #print "alloc [] -> %.3f sec" % dur


        startt = time.time()
        for at in atoms :
            #C = at.xformCoord() - minP
            C = at.coord() - minP
            i = int ( numpy.floor ( C[0]/maxD ) )
            j = int ( numpy.floor ( C[1]/maxD ) )
            k = int ( numpy.floor ( C[2]/maxD ) )

            boxes[k][j][i].append  ( at )

        dur = time.time() - startt
        #print "put -> %.3f sec" % dur


        if 0 :
            from CGLutil.AdaptiveTree import AdaptiveTree
            startt = time.time()
            points = _multiscale.get_atom_coordinates ( atoms, transformed = False )
            self.tree = AdaptiveTree ( points.tolist(), atoms, maxD)
            dur = time.time() - startt
            #print "tree -> %.3f sec" % dur

        self.boxes = boxes
        self.O = minP
        self.D = maxD
        self.n = n
        self.ni = n0
        self.nj = n1
        self.nk = n2



    def AtsNearPt ( self, P ) :

        #startt = time.time()

        ats = []
        #atsByDist = []

        C = P - self.O
        i = int ( numpy.floor ( C[0]/self.D ) )
        j = int ( numpy.floor ( C[1]/self.D ) )
        k = int ( numpy.floor ( C[2]/self.D ) )

        for kk in (k-1, k, k+1) :
            if kk < 0 or kk >= self.nk :
                continue

            jboxes = self.boxes[kk]

            for jj in (j-1, j, j+1) :
                if jj < 0 or jj >= self.nj :
                     continue

                iboxes = jboxes[jj]

                for ii in (i-1, i, i+1) :
                    if ii < 0 or ii >= self.ni :
                        continue

                    atsInBox = iboxes[ii]
                    #ats.extend ( atsInBox )

                    for at in atsInBox :
                        v = at.xformCoord() - P
                        if v.length < self.D :
                            ats.append ( at )
                            #atsByDist.append ( [v.length, at] )


        #dur = time.time() - startt
        #print "boxes -> %.5f sec - %d" % (dur, len(ats))

        #atsByDist.sort()
        #for d, at in atsByDist : print "  %.3f -- %s.%d.%s.%s" % (d, at.name,at.residue.id.position,at.residue.type,at.residue.id.chainId)

        if self.tree :
            startt = time.time()

            #pt = at.coord()
            #vPt = numpy.array ( P.data() )
            atsNear = self.tree.searchTree ( [P[0], P[1], P[2]], self.D )

            cats = []
            catsByDist = []
            for at in atsNear :
                v = at.xformCoord() - P
                if v.length < self.D :
                    cats.append ( at )
                    catsByDist.append ( [v.length, at] )

            dur = time.time() - startt
            print "tree -> %.5f sec - %d" % (dur, len(cats))

        #catsByDist.sort()
        #for d, at in catsByDist : print "  %.3f -- %s.%d.%s.%s" % (d, at.name,at.residue.id.position,at.residue.type,at.residue.id.chainId)

        return ats





    def AtsNearPtLocal ( self, P ) :

        #startt = time.time()

        ats = []
        #atsByDist = []

        C = P - self.O
        i = int ( numpy.floor ( C[0]/self.D ) )
        j = int ( numpy.floor ( C[1]/self.D ) )
        k = int ( numpy.floor ( C[2]/self.D ) )

        for kk in (k-1, k, k+1) :
            if kk < 0 or kk >= self.nk :
                continue

            jboxes = self.boxes[kk]

            for jj in (j-1, j, j+1) :
                if jj < 0 or jj >= self.nj :
                     continue

                iboxes = jboxes[jj]

                for ii in (i-1, i, i+1) :
                    if ii < 0 or ii >= self.ni :
                        continue

                    atsInBox = iboxes[ii]
                    #ats.extend ( atsInBox )

                    for at in atsInBox :
                        v = at.coord() - P
                        if v.length < self.D :
                            ats.append ( [at, v] )
                            #atsByDist.append ( [v.length, at] )

        return ats








# end
