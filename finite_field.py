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
