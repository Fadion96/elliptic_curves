import argparse
from datetime import datetime
import json
import os
import numpy as np 
from finite_field import F
from affine_point import AffinePoint
import random
import operator
from functools import reduce
from projective_point import ProjectivePoint

parser = argparse.ArgumentParser()
parser.add_argument('bit_length')
parser.add_argument('-p', dest='projective', const=True, default=False, nargs='?', help='Projective point')
parser.add_argument('-f', dest='file', const=True, default=False, nargs='?', help='file')
args = parser.parse_args()



# if not args.file:
#     os.system(f'sage ec-prime-order.sage {args.bit_length}')

with open(f'curve_params_{args.bit_length}.json', "r") as f:
    params = json.load(f)
print(params)

F_ord_ec = F(params["curveOrder"])
EC = AffinePoint(*params['invariants'], params['fieldOrder'])
P = EC(*params["basePoint"][:-1])

# s = random.randint(2, params["curveOrder"] - 1)
s = random.randint(2, 2**520 - 1)
print(s)
bit_length = 520
h = 10
v = 13
a = bit_length // h
b = a // v
print(a, b)

Gs = []
for j in range(0, v):
    if j == 0:
        g_j = [EC(None, None)]
        for u in range(1, 2 ** h):
            u_bits = [0] * (h - u.bit_length()) + [int(x) for x in bin(u)[2:]]
            u_bits.reverse()
            XD = []
            XD = [(P*(2**(i * a)))*u_bit for i, u_bit in enumerate(u_bits)]
            DD = reduce(operator.add, XD)
            g_j.append(DD)
        Gs.append(g_j)
    else:
        g_j = [EC(None, None)]
        for u in range(1, 2 ** h):
            DD = Gs[0][u] * (2**(j*b))
            g_j.append(DD)
        Gs.append(g_j)
# print(Gs)
print("Done precomp")
bits = [0] * (bit_length - s.bit_length()) + [int(x) for x in bin(s)[2:]]
bits.reverse()
# print(bits)
# print(len(bits))
a_bits = np.array_split(bits, h)
# print(a_bits)
b_bits = [np.array_split(a, v) for a in a_bits]
b_bits = [[list(zzz)for zzz in ccc] for ccc in b_bits]
# print(b_bits)
start = datetime.now()
R = EC(None, None)
for k in range(b - 1, -1, -1):
    R = 2 * R
    for j in range(v - 1, -1, -1):
        I_jk = sum([b_bits[i][j][k] * (2**i) for i in range(0, h)])
        # print(I_jk)
        R = R + Gs[j][I_jk]
print(datetime.now() - start)
print(R)
start = datetime.now()
Q = P * s
print(datetime.now() - start)
assert R == Q
        
# if args.projective:
# print('Projective')
# PP = ProjectivePoint(*params['invariants'], params['fieldOrder'])
# P1 = PP(*params["basePoint"])
# Q1 = s * P1

# print('Affine')
# EC = AffinePoint(*params['invariants'], params['fieldOrder'])
# P = EC(*params["basePoint"][:-1])
# Q = s * P
