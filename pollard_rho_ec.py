import argparse
from datetime import datetime
import json
import os

from finite_field import F
from affine_point import AffinePoint
import random

from projective_point import ProjectivePoint


def step1(A, alpha, beta, g, y, F_ord_p):
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
        A, alpha_A, beta_A = step1(A, alpha_A, beta_A, g, y, F_ord_p)
        B, alpha_B, beta_B = step1(B, alpha_B, beta_B, g, y, F_ord_p)
        B, alpha_B, beta_B = step1(B, alpha_B, beta_B, g, y, F_ord_p)
        if A == B:
            x = (alpha_B - alpha_A) / (beta_A - beta_B)
            return x 


parser = argparse.ArgumentParser()
parser.add_argument('bit_length')
parser.add_argument('-p', dest='projective', const=True, default=False, nargs='?', help='Projective point')
parser.add_argument('-f', dest='file', const=True, default=False, nargs='?', help='file')
args = parser.parse_args()

# a = 40798
# b = 14047
# mod = 62071
# x = 2928
# y = 42354
# ord_ec = 62039

if not args.file:
    os.system(f'sage ec-prime-order.sage {args.bit_length}')

with open(f'curve_params_{args.bit_length}.json', "r") as f:
    params = json.load(f)
print(params)

F_ord_ec = F(params["curveOrder"])
for s in [6]:
    print(s)
    if args.projective:
        print('Projective')
        PP = ProjectivePoint(*params['invariants'], params['fieldOrder'])
        P1 = PP(*params["basePoint"])
        Q1 = s * P1
        print(P1.x, P1.y, P1.z, Q1.x,Q1.y, Q1.z)
        start = datetime.now()
        x = pollard_rho(P1, Q1, F_ord_ec)
    else:
        print('Affine')
        EC = AffinePoint(*params['invariants'], params['fieldOrder'])
        P = EC(*params["basePoint"][:-1])
        Q = s * P
        start = datetime.now()
        x = pollard_rho(P, Q, F_ord_ec)

    print(datetime.now() - start)
    assert s == x.int
    print(x)
