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
            s_1 = self.y * other.z  # s_1 = y_1z_2
            s_2 = other.y * self.z  # s_2 = y_2z_1
            u_1 = self.x * other.z  # u_1 = x_1z_2
            u_2 = other.x * self.z  # u_2 = x_2z_1
            if u_1 == u_2:  # u
                if s_1 == s_2:  # s
                    return self.double()
                else:
                    return ProjectivePoint(None, None, None)
            else:
                w = self.z * other.z  # w = z_1z_2
                r = s_2 - s_1  # r = s_2 - s_1
                p = u_2 - u_1   
                p_sqr = p * p
                rrw = r * r * w
                p_sqr_u1_u2 = p_sqr * (u_1 + u_2)
                p_cube = p * p_sqr
                rx = p * (rrw - p_sqr_u1_u2)
                ry = r * (self.field(-2)*rrw + self.field(3) * p_sqr_u1_u2) - p_cube * (s_1 + s_2)
                rz = p_cube * w
                return ProjectivePoint(self.field(2)*rx, ry, self.field(2)*rz)

        def double(self):
            if self.is_zero() or self.y == self.field(0):
                return ProjectivePoint(None, None, None)
            else:
                two = self.field(2)
                w = self.x * self.x * self.field(3) + self.field(a) * self.z * self.z
                s = self.y * self.z * two
                b = s * self.x * self.y * two
                h = w * w - b * two
                rx = h * s  
                ry = w * (b - h) - s * s * self.y * self.y * two
                rz = s * s * s
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
