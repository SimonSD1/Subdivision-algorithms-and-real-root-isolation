restart:

with(RootFinding):

filename := "src/bin/tmp_bound.txt":
f := readline(filename):
P := parse(f):

sols := Isolate(P, x):

# Ouvrir le fichier de sortie
out := fopen("src/bin/maple_roots_bound.txt", WRITE):

# Exporter les valeurs numériques (droites) des racines
for sol in sols do
    fprintf(out, "%.16f\n", evalf(rhs(sol))):
end do:

fclose(out):
