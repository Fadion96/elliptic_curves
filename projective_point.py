from finite_field import F

from affine_point import AffinePoint


def ProjectivePoint(a, b, mod):
    class ProjectivePoint:
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
            t0 = self.y * other.z
            t1 = other.y * self.z
            u0 = self.x * other.z
            u1 = other.x * self.z
            if u0 == u1:
                if t0 == t1:
                    return self.double()
                else:
                    return ProjectivePoint(None, None, None)
            else:
                t = t0 - t1
                u = u0 - u1
                u2 = u * u
                v = self.z * other.z
                w = t * t * v - u2 * (u0 + u1)
                u3 = u * u2
                rx = u * w
                ry = t * (u0 * u2 - w) - t0 * u3
                rz = u3 * v
                return ProjectivePoint(rx, ry, rz)

        def double(self):
            if self.is_zero() or self.y == self.field(0):
                return ProjectivePoint(None, None, None)
            else:
                two = self.field(2)
                t = self.x * self.x * self.field(3) + self.field(a) * self.z * self.z
                u = self.y * self.z * two
                v = u * self.x * self.y * two
                w = t * t - v * two
                rx = u * w
                ry = t * (v - w) - u * u * self.y * self.y * two
                rz = u * u * u
                return ProjectivePoint(rx, ry, rz)

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

        def to_affine_point(self):
            if self.is_zero():
                return AffinePoint(a, b, mod)(None, None)
            else:
                return AffinePoint(a, b, mod)(self.x / self.z, self.y / self.z)

    return ProjectivePoint
