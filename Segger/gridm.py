

import chimera
import numpy
import time
import _contour
import _multiscale



class Grid (object) :

    boxes = {} # map of boxes
    D = None # dimension of boxes
    tree = None # for comparison

    def __init__ (self) :
        self.boxes = {} # array of boxes
        self.D = None # dimension of boxes
        self.tree = None # for comparison


    def FromMols (self, mols, maxD) :

        ats = []
        for m in mols :
            ats.extend ( m.atoms )

        self.FromAtoms ( ats, maxD )


    def FromAtoms ( self, atoms, maxD ) :

        boxes = {}
        if 0 :
            startt = time.time()
            boxes = [ [] ] * N
            dur = time.time() - startt
            print "alloc all -> %.3f sec" % dur

        startt = time.time()
        for at in atoms :
            C = at.xformCoord()
            i = int ( numpy.floor ( C[0]/maxD ) )
            j = int ( numpy.floor ( C[1]/maxD ) )
            k = int ( numpy.floor ( C[2]/maxD ) )

            boxesk = None
            if k in boxes :
                boxesk = boxes[k]
            else :
                boxesk = {}
                boxes[k] = boxesk

            boxesj = None
            if j in boxesk :
                boxesj = boxesk[j]
            else :
                boxesj = {}
                boxesk[j] = boxesj

            boxesi = None
            if i in boxesj :
                boxesi = boxesj[i]
            else :
                boxesi = []
                boxesj[i] = boxesi

            boxesi.append  ( at )

        #dur = time.time() - startt
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
        self.D = maxD




    def FromAtomsLocal (self, atoms, maxD) :

        boxes = {}
        #startt = time.time()

        for at in atoms :
            #C = at.xformCoord() - minP
            C = at.coord()
            i = int ( numpy.floor ( C[0]/maxD ) )
            j = int ( numpy.floor ( C[1]/maxD ) )
            k = int ( numpy.floor ( C[2]/maxD ) )

            boxesk = None
            if k in boxes :
                boxesk = boxes[k]
            else :
                boxesk = {}
                boxes[k] = boxesk

            boxesj = None
            if j in boxesk :
                boxesj = boxesk[j]
            else :
                boxesj = {}
                boxesk[j] = boxesj

            boxesi = None
            if i in boxesj :
                boxesi = boxesj[i]
            else :
                boxesi = []
                boxesj[i] = boxesi

            boxesi.append  ( at )

        #dur = time.time() - startt
        #print "put -> %.3f sec" % dur


        if 0 :
            from CGLutil.AdaptiveTree import AdaptiveTree
            startt = time.time()
            points = _multiscale.get_atom_coordinates ( atoms, transformed = False )
            self.tree = AdaptiveTree ( points.tolist(), atoms, maxD)
            dur = time.time() - startt
            #print "tree -> %.3f sec" % dur

        self.boxes = boxes
        self.D = maxD



    def AddAtomsLocal (self, atoms) :

        for at in atoms :
            #C = at.xformCoord() - minP
            C = at.coord()
            i = int ( numpy.floor ( C[0]/self.D ) )
            j = int ( numpy.floor ( C[1]/self.D ) )
            k = int ( numpy.floor ( C[2]/self.D ) )

            boxesk = None
            if k in self.boxes :
                boxesk = self.boxes[k]
            else :
                boxesk = {}
                boxes[k] = boxesk

            boxesj = None
            if j in boxesk :
                boxesj = boxesk[j]
            else :
                boxesj = {}
                boxesk[j] = boxesj

            boxesi = None
            if i in boxesj :
                boxesi = boxesj[i]
            else :
                boxesi = []
                boxesj[i] = boxesi

            boxesi.append  ( at )



    def AtsNearPt ( self, C ) :

        #startt = time.time()

        ats = []
        #atsByDist = []

        i = int ( numpy.floor ( C[0]/self.D ) )
        j = int ( numpy.floor ( C[1]/self.D ) )
        k = int ( numpy.floor ( C[2]/self.D ) )

        for kk in (k-1, k, k+1) :

            if not kk in self.boxes : continue
            jboxes = self.boxes[kk]

            for jj in (j-1, j, j+1) :

                if not jj in jboxes : continue
                iboxes = jboxes[jj]

                for ii in (i-1, i, i+1) :

                    if not ii in iboxes : continue
                    atsInBox = iboxes[ii]
                    #ats.extend ( atsInBox )

                    for at in atsInBox :
                        v = at.xformCoord() - C
                        if v.length < self.D :
                            ats.append ( [at, v] )
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





    def AtsNearPtLocal ( self, C ) :

        #startt = time.time()
        ats = []
        #atsByDist = []

        i = int ( numpy.floor ( C[0]/self.D ) )
        j = int ( numpy.floor ( C[1]/self.D ) )
        k = int ( numpy.floor ( C[2]/self.D ) )

        for kk in (k-1, k, k+1) :

            if not kk in self.boxes : continue
            jboxes = self.boxes[kk]

            for jj in (j-1, j, j+1) :

                if not jj in jboxes : continue
                iboxes = jboxes[jj]

                for ii in (i-1, i, i+1) :

                    if not ii in iboxes : continue
                    atsInBox = iboxes[ii]
                    #ats.extend ( atsInBox )

                    for at in atsInBox :
                        v = at.coord() - C
                        if v.length < self.D :
                            ats.append ( [at, v] )
                            #atsByDist.append ( [v.length, at] )

        return ats








# end
