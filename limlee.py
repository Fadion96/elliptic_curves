import argparse
import math
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
import pprint
from multiprocessing.pool import ThreadPool

def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def find_opt(l,storage):
    opt_a, opt_b = None, None
    min_comp = np.inf
    for a in range(2, l + 1):
        for b in range(1, a):
            h = math.ceil(l / a)
            v = math.ceil(a / b)
            a_last = l - a *(h-1)
            v_last = math.ceil(a_last / b)
            comp = (b - 1)  + ((2**(h - 1) - 1)/(2**(h - 1))) * (a - a_last) + ((2**h  - 1)/(2**h )) * a_last - 1
            number_of_precomp = (2 ** h - 1) * v_last + (2 ** (h-1) - 1) * (v - v_last)
            if comp < min_comp and number_of_precomp <= storage:
                opt_a, opt_b = a, b
                min_comp = comp
    return opt_a, opt_b

def _precomp0(P, v, h, a, start, end):
    g_j = []
    for u in range(start, end):
        u_bits = [0] * (h - u.bit_length()) + [int(x) for x in bin(u)[2:]]
        u_bits.reverse()
        R = reduce(operator.add, [(P*(2**(i * a)))*u_bit for i, u_bit in enumerate(u_bits)])
        g_j.append(R)
    return g_j

def _precomp(P, v, h, a, start, end, G_0, j):
    g_j = []
    for u in range(start, end):
        R = G_0[u] * (2**(j*b))
        g_j.append(R)
    return g_j

def precomp(P, v, h, a):
    G = []
    ranges = np.array_split(list(range(1,2**h + 1)), 4)
    ranges = [(ranges[i][0], ranges[i+1][0]) for i in range(len(ranges) - 1)] + [(ranges[-1][0], 2**h)]
    pool = ThreadPool(processes=4)
    params = [(P, v, h, a, *rang) for rang in ranges]
    g_j = pool.starmap(_precomp0, params)
    g_j = reduce(operator.iconcat, g_j, [EC(None, None)])
    G.append(g_j)
    for j in range(1, v_last):
        G_0 = G[0]
        params = [(P, v, h, a, *rang, G_0, j) for rang in ranges]
        g_j = pool.starmap(_precomp, params)
        g_j = reduce(operator.iconcat, g_j, [EC(None, None)])
        G.append(g_j)
    ranges = np.array_split(list(range(1,2**(h-1) + 1)), 4)
    ranges = [(ranges[i][0], ranges[i+1][0]) for i in range(len(ranges) - 1)] + [(ranges[-1][0], 2**(h-1))]
    for j in range(0, v - v_last):
        G_0 = G[0]
        params = [(P, v, h-1, a, *rang, G_0, v_last + j) for rang in ranges]
        g_j = pool.starmap(_precomp, params)
        g_j = reduce(operator.iconcat, g_j, [EC(None, None)])
        G.append(g_j)
    return G

def multiplication(P, e, a, b, v, v_last, b_last, G):
    bits = [0] * (bit_length - s.bit_length()) + [int(x) for x in bin(s)[2:]]
    bits.reverse()
    a_bits = list(chunks(bits, a))
    b_bits = [list(chunks(a, b)) for a in a_bits]
    if b == b_last:
        b_bits[-1] += [[None]] * (v - v_last)
    else:
        b_bits[-1][-1] += [None] * (b - b_last)
        b_bits[-1] += [[None]*b] 
    start = datetime.now()
    R = EC(None, None)
    for k in range(b - 1, -1, -1):
        R = 2 * R
        for j in range(v - 1, -1, -1):
            I_jk = 0
            for i in range(0,h):
                if b_bits[i][j][k] is not None:
                    I_jk += b_bits[i][j][k] * (2**i)
                else:
                    continue
            R = R + G[j][I_jk]
    print("Lim/lee:", datetime.now() - start)
    return R

parser = argparse.ArgumentParser()
parser.add_argument('bit_length')
parser.add_argument('-p', dest='projective', const=True, default=False, nargs='?', help='Projective point')
parser.add_argument('-f', dest='file', const=True, default=False, nargs='?', help='file')
args = parser.parse_args()


with open(f'curve_params_{args.bit_length}.json', "r") as f:
    params = json.load(f)
print(params)

F_ord_ec = F(params["curveOrder"])
EC = AffinePoint(*params['invariants'], params['fieldOrder'])
P = EC(*params["basePoint"][:-1])

s = random.randint(2, params["curveOrder"] - 1)
print("Scalar:", s)

bit_length = int(args.bit_length)
a, b = find_opt(bit_length, 1000)
h = math.ceil(bit_length / a)
v = math.ceil(a / b)
a_last = bit_length - a *(h-1)
v_last = math.ceil(a_last / b)
b_last = a_last - b * (v_last - 1)
print(f'a: {a}, b: {b}, h: {h}, v: {v}, a_last: {a_last}, b_last: {b_last}, v_last:{v_last}')
G = precomp(P, v, h, a)
R = multiplication(P, s, a, b, v, v_last, b_last, G)
start = datetime.now()
Q = P * s
print("Normal:", datetime.now() - start)
assert R == Q
