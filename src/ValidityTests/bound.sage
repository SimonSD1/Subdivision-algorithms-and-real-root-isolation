import sys
from sage.all import *

def Lagrange_bound(poly) :
    sum = 0
    length = poly.degree()
    for i in range(length) :
        sum += abs(poly[i])
    
    potentialBound = ceil(sum/abs(poly[length]))

    return max(1,potentialBound)



def Cauchy_bound(poly) :

    coeffList = poly.list()
    max_coeff_abs = max(map(abs, coeffList))

    return 1 + max_coeff_abs/abs(poly[poly.degree()])



def local_max_bound_implementation(poly):
    deg = poly.degree()
    if deg < 1:
        return 0  
    j = deg
    t = 1

    if poly[deg] < 0:
        poly = -poly  

    tempub = 0
    ub = 0

    for i in range(deg, -1, -1):
        coef_i = poly[i]
        coef_j = poly[j]

        if coef_i < 0:
            tempub = (2^t * abs(coef_i)) / coef_j
            k = j - i

            tempub = floor(tempub^(1/k)) + 1
            
            if tempub > ub:
                ub = tempub
                t += 1
        else:
            if coef_i > coef_j:
                j = i
                t = 1
    return ub







# Read polynomial from command line
if len(sys.argv) < 2:
    print("Usage: sage script.sage 'polynomial'")
    sys.exit(1)

R = PolynomialRing(ZZ, 'x')
x = R.gen()

# Convert input string into polynomial
input_poly_str = sys.argv[1]
poly = R(input_poly_str)

print("Input Polynomial:", poly)
print("Lagrange Bound:", Lagrange_bound(poly))
print("Cauchy Bound:", Cauchy_bound(poly))
print("Local Max Bound:", local_max_bound_implementation(poly))