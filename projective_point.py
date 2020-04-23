from finite_field import F

from affine_point import AffinePoint


def ProjectivePoint(a, b, mod):
    class ProjectivePoint:
        __slots__ = ['x', 'y', 'z', 'field']
        def __init__(self, x, y, z):
            self.field = F(mod)
            if isinstance(x, int):
                self.x = self.field(x)
            else:
                self.x = x
            if isinstance(y, int):
                self.y = self.field(y)
            else:
                self.y = y
            if isinstance(z, int):
                self.z = self.field(z)
            else:
                self.z = z

        def is_zero(self):
            return self.x is None

        def is_on_curve(self):
            return not self.is_zero() and \
                   self.y * self.y * self.z == \
                   self.x * self.x * self.x + self.field(a) * self.x * self.z * self.z + self.field(b) * self.z * self.z * self.z

        def __add__(self, other):
            if self.is_zero():
                return other
            elif other.is_zero():
                return self
            sx, ox = self.x.int, other.x.int
            sy, oy = self.y.int, other.y.int
            sz, oz = self.z.int, other.z.int
            t0 = sy * oz
            t1 = oy * sz
            u0 = sx * oz
            u1 = ox * sz
            if u0 % mod == u1 % mod:
                if t0 % mod == t1 % mod:
                    return self.double()
                else:
                    return ProjectivePoint(None, None, None)
            else:
                t = t0 - t1
                u = u0 - u1
                u2 = u * u
                v = sz * oz
                w = t * t * v - u2 * (u0 + u1)
                u3 = u * u2
                rx = u * w
                ry = t * (u0 * u2 - w) - t0 * u3
                rz = u3 * v
                return ProjectivePoint(self.field(rx), self.field(ry), self.field(rz))

        def double(self):
            if self.is_zero() or self.y == self.field(0):
                return ProjectivePoint(None, None, None)
            else:
                two = 2
                sx = self.x.int
                sy = self.y.int
                sz = self.z.int
                t = sx * sx * 3 + a * sz * sz
                u = sy * sz * two
                v = u * sx * sy * two
                w = t * t - v * two
                rx = u * w
                ry = t * (v - w) - u * u * sy * sy * two
                rz = u * u * u
                return ProjectivePoint(self.field(rx), self.field(ry), self.field(rz))

        def __neg__(self):
            if self.is_zero():
                return self
            else:
                return ProjectivePoint(self.x, -self.y, self.z)

        def __sub__(self, other):
            return self + -other

        def __mul__(self, n):
            if n < 0:
                return -self * -n
            result = ProjectivePoint(None, None, None)
            temp = self
            while n != 0:
                if n & 1 != 0:
                    result += temp
                temp = temp.double()
                n >>= 1
            return result

        __rmul__ = __mul__

        def __eq__(self, other):
            if self.is_zero() or other.is_zero():
                return self.is_zero() and other.is_zero()
            else:
                return (self.x * other.z, self.y * other.z) == (other.x * self.z, other.y * self.z)

        def __ne__(self, other):
            return not (self == other)

        def __str__(self):
            return f'{self.x} {self.y} {self.z}'

        def to_affine_point(self):
            if self.is_zero():
                return AffinePoint(a, b, mod)(None, None)
            else:
                return AffinePoint(a, b, mod)(self.x / self.z, self.y / self.z)

    return ProjectivePoint
