import sys
from sage.all import *
from cypari2 import pari

pari.allocatemem(2**31)

R = PolynomialRing(ZZ, 'x')
x = R.gen()

# Lire le polynôme depuis le fichier
with open("src/bin/tmp_bound.txt", "r") as f:
    poly_str = f.read().strip()

P = R(poly_str)

try:
    roots = P.roots(ring=RR, multiplicities=False)
    positive_roots = [r for r in roots if r > 0]
    if positive_roots:
        max_positive_root = max(positive_roots)
    else:
        max_positive_root = 0
except Exception as e:
    print("Erreur dans le calcul des racines :", e)
    max_positive_root = -1

# Sauvegarde dans un autre fichier pour éviter l'écrasement du polynôme
with open("src/bin/sage_bound.txt", "w") as o:
    o.write(str(ceil(max_positive_root)))
