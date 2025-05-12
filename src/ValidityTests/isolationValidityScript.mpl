restart:
randomize():

WriteUFile := proc(f, str)
    local fd, deg, i;
    deg := degree(f);
    fd := fopen(str, WRITE);
    fprintf(fd, "%d  ", deg + 1);
    for i from 0 to deg do
        fprintf(fd, "%a ", coeff(f, x, i));
    end do;
    fclose(fd);
end proc:

ExportRoots := proc(roots, filepath)
    local fd, r;
    fd := fopen(filepath, WRITE);
    for r in roots do
        fprintf(fd, "%.20f\n", evalf(rhs(r)));
    end do;
    fclose(fd);
end proc:

with(RootFinding):
f := randpoly([x], degree = 2000, dense, coeffs = rand(-2^1500 .. 2^1500)):
WriteUFile(f, "D:/Users/Joshua/Documents/cours-FAC/M1-S2/Subdivision-algorithms-and-real-root-isolation/src/bin/test.out"):
sols := Isolate(f, x):
print(sols);
ExportRoots(sols, "D:/Users/Joshua/Documents/cours-FAC/M1-S2/Subdivision-algorithms-and-real-root-isolation/src/bin/maple_roots.txt"):

