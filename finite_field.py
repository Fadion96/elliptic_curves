# def egcd(a, b):
#     lastremainder, remainder = abs(a), abs(b)
#     x, lastx, y, lasty = 0, 1, 1, 0
#     while remainder:
#         lastremainder, (quotient, remainder) = remainder, divmod(lastremainder, remainder)
#         x, lastx = lastx - quotient*x, x
#         y, lasty = lasty - quotient*y, y
#     return lastremainder, lastx * (-1 if a < 0 else 1), lasty * (-1 if b < 0 else 1)

# def modinv(b, n):
# 	g, x, _ = egcd(b, n)
# 	if g == 1:
# 		return x % n



def F(p):
    class F:
        def __init__(self, x):
            self.int = x % p

        def __str__(self):
            return str(self.int)

        __repr__ = __str__

        def __eq__(self, b):
            return self.int == b.int

        def __ne__(self, b):
            return self.int != b.int

        def __add__(self, b):
            return F(self.int + b.int)

        def __sub__(self, b):
            return F(self.int - b.int)

        def __mul__(self, b):
            return F(self.int * b.int)

        def __neg__(self):
            return F(p - self.int)

        def __truediv__(self, b):
            return self*F(pow(b.int, p-2, p))

    return F