from finite_field import F
from affine_point import AffinePoint
import random


def step(A, alpha, beta, g, y, F_ord_p):
    if A.is_zero() or A.x.int % 3 == 1:
        return (A + y), alpha, (beta + F_ord_p(1))
    elif A.x.int % 3 == 0:
        return (A + A), (F_ord_p(2) * alpha), (F_ord_p(2) * beta)
    elif A.x.int % 3 == 2:
        return (g + A), (alpha + F_ord_p(1)), beta

def pollard_rho(g, y, F_ord_p):
    A = g
    B = g
    alpha_A, beta_A = F_ord_p(1), F_ord_p(0)
    alpha_B, beta_B = F_ord_p(1), F_ord_p(0)
    while True:
        A, alpha_A, beta_A = step(A, alpha_A, beta_A, g, y, F_ord_p)
        B, alpha_B, beta_B = step(B, alpha_B, beta_B, g, y, F_ord_p)
        B, alpha_B, beta_B = step(B, alpha_B, beta_B, g, y, F_ord_p)
        if A == B:
            x = (alpha_B - alpha_A) / (beta_A - beta_B)
            return x 

a = 40798
b = 14047
mod = 62071
x = 2928
y = 42354
ord_ec = 62039


F_ord_ec = F(ord_ec)

EC = AffinePoint(a, b, mod)

P = EC(x, y)
s = random.randint(2, ord_ec - 1)
Q = s * P
print(s)
print(Q)
x = pollard_rho(P, Q, F_ord_ec)
print(x)