from finite_field import F
import gensafeprime
import random

def step(A, alpha, beta, g, y, F_ord_p):
    if A.int % 3 == 0:
        return (A * A), (F_ord_p(2) * alpha), (F_ord_p(2) * beta)
    elif A.int % 3 == 1:
        return (A * y), alpha, (beta + F_ord_p(1))
    elif A.int % 3 == 2:
        return (g * A), (alpha + F_ord_p(1)), beta

def pollard_rho(g, y, F_p, F_ord_p):
    A = F_p(1)
    B = F_p(1)
    alpha_A, beta_A = F_ord_p(0), F_ord_p(0)
    alpha_B, beta_B = F_ord_p(0), F_ord_p(0)
    while True:
        A, alpha_A, beta_A = step(A, alpha_A, beta_A, g, y, F_ord_p)
        B, alpha_B, beta_B = step(B, alpha_B, beta_B, g, y, F_ord_p)
        B, alpha_B, beta_B = step(B, alpha_B, beta_B, g, y, F_ord_p)
        if A == B:
            x = (alpha_B - alpha_A) / (beta_A - beta_B)
            return x 
        

for i in range(50):
    safeprime = gensafeprime.generate(30)
    F_p = F(safeprime)
    F_ord_p = F((safeprime-1)//2)
    y = random.randint(2, safeprime - 1) ** 2
    x = pollard_rho(F_p(4), F_p(y), F_p, F_ord_p)
    print(i , F_p(pow(4, x.int, safeprime)) == F_p(y))

