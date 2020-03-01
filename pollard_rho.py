def egcd(a, b):
    lastremainder, remainder = abs(a), abs(b)
    x, lastx, y, lasty = 0, 1, 1, 0
    while remainder:
        lastremainder, (quotient, remainder) = remainder, divmod(lastremainder, remainder)
        x, lastx = lastx - quotient*x, x
        y, lasty = lasty - quotient*y, y
    return lastremainder, lastx * (-1 if a < 0 else 1), lasty * (-1 if b < 0 else 1)

def modinv(b, n):
	g, x, _ = egcd(b, n)
	if g == 1:
		return x % n

def step(A, alpha, beta, g, y, p, ord_p):
    if A % 3 == 0:
        return (A * A) % p, (2 * alpha) % ord_p, (2 * beta) % ord_p
    elif A % 3 == 2:
        return (A * y) % p, alpha, (beta + 1) % ord_p
    elif A % 3 == 1:
        return (g * A) % p, (alpha + 1) % ord_p, beta

def pollard_rho(g, y, p, ord_p):
    A = 1
    B = 1
    alpha_A, beta_A = 0, 0
    alpha_B, beta_B = 0, 0
    while True:
        A, alpha_A, beta_A = step(A, alpha_A, beta_A, g, y, p, ord_p)
        B, alpha_B, beta_B = step(B, alpha_B, beta_B,g, y, p, ord_p)
        B, alpha_B, beta_B = step(B, alpha_B, beta_B,g, y, p, ord_p)
        if A == B:
            x = (alpha_B - alpha_A) * modinv(beta_A - beta_B, ord_p) % ord_p
            return x 
        
    
x = pollard_rho(2, 5, 1019, 509)
assert pow(2,x, 1019) == 5


    