from finite_field import F



def AffinePoint(a, b, mod):
    class AffinePoint:
        __slots__ = ['field', 'x', 'y']
        def __init__(self, x, y):
            self.field = F(mod)
            if isinstance(x, int):
                self.x = self.field(x)
            else:
                self.x = x
            if isinstance(y, int):
                self.y = self.field(y)
            else:
                self.y = y            

        def is_zero(self):
            return self.x is None

        def __add__(self, other):
            if self.is_zero():
                return other
            elif other.is_zero():
                return self
            elif self.x == other.x:
                if self.y == other.y:
                    return self.double()
                else:
                    return AffinePoint(None, None)
            else:
                temp = (other.y - self.y) / (other.x - self.x)
                rx = temp * temp - self.x - other.x
                ry = - self.y + temp * (self.x - rx) 
                return AffinePoint(rx, ry)
            
        def double(self):
            if self.is_zero() or self.y == self.field(0):
                return AffinePoint(None, None)
            else:
                temp = (self.x * self.x * self.field(3) + self.field(a)) / (self.y * self.field(2))
                rx = temp * temp - self.x * self.field(2)
                ry = temp * (self.x - rx) - self.y
                return AffinePoint(rx, ry)

        def __neg__(self):
            if self.is_zero():
                return self
            else:
                return AffinePoint(self.x, -self.y)

        def __mul__(self, n):
            result = AffinePoint(None, None)
            temp = self
            while n != 0:
                if n & 1 != 0:
                    result += temp
                temp = temp.double()
                n >>= 1
            return result

        def __str__(self):
            if self.is_zero():
                return "(Zero)"
            else:
                return f"({self.x}, {self.y})"

        __rmul__ = __mul__
        __repr__ = __str__

        def __eq__(self, other):
            if self.is_zero() or other.is_zero():
                return self.is_zero() and other.is_zero()
            else:
                return (self.x, self.y) == (other.x, other.y)

    return AffinePoint
