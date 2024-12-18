
import chimera
import numpy


class Quaternion :

    def __init__ ( self, s=1.0, v=chimera.Vector(0,0,0) ) :
        self.s = s
        self.v = v

    def length (self) :
        return numpy.sqrt ( (self.s*self.s) + self.v.sqlength() )


    def rotation (self, angDegrees, axis) :
        angRad = 0.5 * angDegrees * numpy.pi / 180.0
        self.s = numpy.cos ( angRad )
        self.v = axis * numpy.sin ( angRad )


    def inverse ( self ) :
        return Quaternion ( self.s, self.v * -1.0 )


    def fromXform ( self, xf ) :

        axis, angle = xf.getRotation ()
        if angle >= -180.0 and angle <= 180.0 :
            self.rotation ( angle, axis )
        elif angle < -180.0 :
            blah
            self.rotation ( angle, axis*-1.0 )
        else :
            blah
            self.rotation ( angle, axis*-1.0 )

        m = numpy.reshape ( xf.getOpenGLMatrix(), (4,4) )
        m = numpy.transpose ( m )
        self.fromMatrix ( m )


    def dot ( self, q ) :
        return self.s * q.s + self.v * q.v

    def angleTo ( self, q2 ) :
        self.normalize()
        q2.normalize()
        return 2.0 * numpy.arccos ( self * q2 )


    def normalize (self) :
        l = self.length()
        if (l > 1e-4) :
            self.s = self.s / l
            self.v = self.v / l
        else :
            raise ("quaternion normalization error")

    def __mul__(self, x) :
        if type(x) == type(1.0) or type(x) == numpy.float64 :
            return Quaternion ( self.s*x, self.v*x )
        else :
            return self.dot ( x )

    def __add__(self, x) :
        return Quaternion ( self.s + x.s, self.v + x.v )

    def __sub__(self, x) :
        return Quaternion ( self.s - x.s, self.v - x.v )

    def __copy__ (self) :
        return Quaternion ( self.s, self.v.__copy__() )

    def Xform (self) :
        #self.normalize()
        s = self.s
        v = self.v
        return chimera.Xform.xform (
            1-2*v.y*v.y-2*v.z*v.z, 2*v.x*v.y-2*s*v.z, 2*v.x*v.z+2*s*v.y, 0,
            2*v.x*v.y+2*s*v.z, 1-2*v.x*v.x-2*v.z*v.z, 2*v.y*v.z-2*s*v.x, 0,
            2*v.x*v.z-2*s*v.y, 2*v.y*v.z+2*s*v.x, 1-2*v.x*v.x-2*v.y*v.y, 0
        )

    def matrix (self) :
        #self.normalize()
        s = self.s
        v = self.v
        return [
            [1-2*v.y*v.y-2*v.z*v.z, 2*v.x*v.y-2*s*v.z, 2*v.x*v.z+2*s*v.y],
            [2*v.x*v.y+2*s*v.z, 1-2*v.x*v.x-2*v.z*v.z, 2*v.y*v.z-2*s*v.x],
            [2*v.x*v.z-2*s*v.y, 2*v.y*v.z+2*s*v.x, 1-2*v.x*v.x-2*v.y*v.y],
        ]


    def fromMatrix ( self, rkRot ) :
        # Algorithm in Ken Shoemake's article in 1987 SIGGRAPH course notes
        # article "Quaternion Calculus and Fast Animation".

        fTrace = rkRot[0,0] + rkRot[1,1] + rkRot[2,2]
        fRoot = 0.0
        if fTrace > 0.0 :
            # |w| > 1/2, may as well choose w > 1/2
            fRoot = numpy.sqrt (fTrace + 1.0)  # 2w
            self.s = 0.5 * fRoot;
            fRoot = 0.5 / fRoot;  # 1/(4w)
            self.v[0] = (rkRot[2,1]-rkRot[1,2])*fRoot;
            self.v[1] = (rkRot[0,2]-rkRot[2,0])*fRoot;
            self.v[2] = (rkRot[1,0]-rkRot[0,1])*fRoot;

        else :
            # |w| <= 1/2
            i = 0
            if rkRot[1,1] > rkRot[0,0] :
                i = 1
            if rkRot[2,2] > rkRot[i,i] :
                i = 2

            j = (i + 1) % 3  # ms_iNext[i];
            k = (j + 1) % 3  # ms_iNext[j];

            fRoot = numpy.sqrt(rkRot[i,i]-rkRot[j,j]-rkRot[k,k]+1.0);

            # Real* apfQuat[3] = { &m_afTuple[1], &m_afTuple[2], &m_afTuple[3] };
            self.v[i] = 0.5 * fRoot # *apfQuat[i] = ((Real)0.5)*fRoot;

            fRoot = 0.5 / fRoot
            self.s = (rkRot[k,j]-rkRot[j,k])*fRoot
            self.v[j] = (rkRot[j,i]+rkRot[i,j])*fRoot  # *apfQuat[j]
            self.v[k] = (rkRot[k,i]+rkRot[i,k])*fRoot  # *apfQuat[k]


def mult (a, b) :
    return Quaternion (a.s*b.s - a.v*b.v, b.v*a.s + a.v*b.s + chimera.cross(a.v,b.v))


def slerp0 (p, q, t) :

    cs = p.dot(q)
    angle = numpy.arccos ( cs )

    if abs (angle) > 0.0 :
        sn = numpy.sin ( angle )
        invSn = 1.0 / sn;
        tAngle = t*angle;
        c0 = numpy.sin(angle - tAngle)*invSn;
        c1 = numpy.sin(tAngle)*invSn;

        #mTuple[0] = coeff0*p.mTuple[0] + coeff1*q.mTuple[0];
        #mTuple[1] = coeff0*p.mTuple[1] + coeff1*q.mTuple[1];
        #mTuple[2] = coeff0*p.mTuple[2] + coeff1*q.mTuple[2];
        #mTuple[3] = coeff0*p.mTuple[3] + coeff1*q.mTuple[3];
        return Quaternion (p.s*c0+q.s*c1, p.v*c0 + q.v*c1)

    else :
        return Quaternion (p.s, chimera.Vector(p.v[0], p.v[1], p.v[2]))


def slerp (v0, v1, t) :

    # http://number-none.com/product/Understanding%20Slerp,%20Then%20Not%20Using%20It/

    #; Inputs are: unit vectors v0 and v1, scalar t
    #; v0 and v1 are linearly independent

    # Quaternion slerp(Quaternion const &v0, Quaternion const &v1, double t) {
    #     // v0 and v1 should be unit length or else
    #     // something broken will happen.
    #
    #     // Compute the cosine of the angle between the two vectors.
    #     double dot = dot_product(v0, v1);
    #
    #     const double DOT_THRESHOLD = 0.9995;
    #     if (dot > DOT_THRESHOLD) {
    #         // If the inputs are too close for comfort, linearly interpolate
    #         // and normalize the result.
    #
    #         Quaternion result = v0 + t*(v1 - v0)
    #         result.normalize();
    #         return result;
    #     }
    #
    #     Clamp(dot, -1, 1);           // Robustness: Stay within domain of acos()
    #     double theta_0 = acos(dot);  // theta_0 = angle between input vectors
    #     double theta = theta_0*t;    // theta = angle between v0 and result
    #
    #     Quaternion v2 = v1 - v0*dot
    #     v2.normalize();              // { v0, v2 } is now an orthonormal basis
    #
    #     return v0*cos(theta) + v2*sin(theta);



    dot = v0.dot(v1)
    #print dot

    if 1 or dot > 0.9995 :
        r = v0 + (v1-v0) * t
        r.normalize()
        return r

    if dot < -1.0 : dot = -1.0
    if dot > 1.0 : dot = 1.0

    theta_0 = numpy.arccos ( dot )
    theta = theta_0*t

    v2 = v1 - v0 * dot
    v2.normalize()

    r = v0 * numpy.cos(theta) + v2 * numpy.sin(theta)

    if 0 :
        # from http://graphics.cs.cmu.edu/nsp/course/15-464/Fall05/assignments/p245-shoemake.pdf
        a0 = numpy.sin( (1-t) * theta_0 ) / numpy.sin(theta_0)
        a1 = numpy.sin ( t * theta_0 ) / numpy.sin ( theta_0 )
        r = v0 * a0 + v1 * a1

    return r
