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



def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def find_opt(l,storage):
    opt_a, opt_b = None, None
    min_comp = math.inf
    opt_storage = None
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
                opt_storage = number_of_precomp
    print("Minimum number of operations:", min_comp)
    print("Predicted storage:", opt_storage)
    return opt_a, opt_b

def precomp(P, v, h, a):
    G = []
    G_0 = [None]
    for u in range(1, 2**h):
        u_bits = [0] * (h - u.bit_length()) + [int(x) for x in bin(u)[2:]]
        u_bits.reverse()
        R = reduce(operator.add, [(P*(2**(i * a)))*u_bit for i, u_bit in enumerate(u_bits)])
        G_0.append(R)
    G.append(G_0)
    for j in range(1, v_last):
        g_j = [None]
        for u in range(1, 2**h):
            R = G_0[u] * (2**(j*b))
            g_j.append(R)
        G.append(g_j)
    for j in range(0, v - v_last):
        g_j = [None]
        for u in range(1, 2**(h-1)):
            R = G_0[u] * (2**((v_last+j)*b))
            g_j.append(R)
        G.append(g_j)
    return G


def multiplication(P, e, a, b, v, v_last, b_last, G):
    bits = [0] * (bit_length - s.bit_length()) + [int(x) for x in bin(s)[2:]]
    bits.reverse()
    a_bits = list(chunks(bits, a))
    b_bits = [list(chunks(a, b)) for a in a_bits]
    start = datetime.now()
    R = EC.get_inf()
    adding = 0
    mult = 0
    for k in range(b - 1, -1, -1):
        R = 2 * R
        mult += 1
        for j in range(v - 1, -1, -1):
            I_jk = 0
            for i in range(0, h):
                try:
                    I_jk += b_bits[i][j][k] * (2**i)                   
                except:
                    continue
            try:
                R = R + G[j][I_jk]
                adding += 1                   
            except:
                continue
           
    print("Lim/lee:", datetime.now() - start)
    print("Additions:", adding)
    print("Multiplications:", mult)
    return R

def normal_mult(P, n):
    adding = 0
    mult = 0
    start = datetime.now()
    result = EC.get_inf()
    temp = P
    while n != 0:
        if n & 1 != 0:
            result += temp
            adding += 1
        temp = temp.double()
        mult += 1
        n >>= 1
    print("Normal:", datetime.now() - start)
    print("Additions:", adding)
    print("Multiplications:", mult)
    return result

parser = argparse.ArgumentParser()
parser.add_argument('bit_length', type=int)
parser.add_argument('-p', dest='projective', const=True, default=False, nargs='?', help='Projective point')
parser.add_argument('-s', dest='storage', type=int, default= 100, help='storage')
args = parser.parse_args()


if __name__ == "__main__":    
    with open(f'curve_params_{args.bit_length}.json', "r") as f:
        params = json.load(f)
    print(params)

    F_ord_ec = F(params["curveOrder"])
    if args.projective:
        EC = ProjectivePoint(*params['invariants'], params['fieldOrder'])
        P = EC(*params["basePoint"])
    else:
        EC = AffinePoint(*params['invariants'], params['fieldOrder'])
        P = EC(*params["basePoint"][:-1])   

    bit_length = args.bit_length
    a, b = find_opt(bit_length, args.storage)
    h = math.ceil(bit_length / a)
    v = math.ceil(a / b)
    a_last = bit_length - a * (h-1)
    v_last = math.ceil(a_last / b)
    b_last = a_last - b * (v_last - 1)
    print(f'a: {a}, b: {b}, h: {h}, v: {v}, a_last: {a_last}, b_last: {b_last}, v_last:{v_last}')

    G = precomp(P, v, h, a)
    print("Used storage:", sum([len(g) for g in G]) - v)
    for _ in range(10):
        print("================")
        s = random.randint(2**(bit_length-1), 2**(bit_length) -1)
        # s = random.randint(2**(bit_length-1), params["curveOrder"] - 1)
        print("Scalar:", s)
        R = multiplication(P, s, a, b, v, v_last, b_last, G)
        print("================")
        Q = normal_mult(P, s)
        assert R == Q
    

