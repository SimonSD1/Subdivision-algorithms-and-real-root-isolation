import matplotlib.pyplot as plt

evolution = [pow(2,i) for i in range(0,17)]

# Temps multiplication fixed size
tps_mult_classic_fixedSize = [1.00, 1.50, 2.20, 3.20, 5.10, 9.50, 23.40, 71.70, 261.60, 1018.30, 3651.40, 13329.50, 51955.60, 215687.10, 942212.75, 3656463.85, 14097578.40]
tps_mult_karatsuba_fixedSize = [5.15, 6.80, 8.90, 14.75, 29.35, 64.55, 167.30, 505.10, 1502.00, 3931.90, 11107.05, 32435.55, 96988.25, 284210.45, 850355.00, 2614987.65, 7821176.30]
tps_mult_choice_fixedSize = [2.80, 5.55, 8.25, 10.30, 12.35, 14.95, 18.55, 53.90, 70.95, 100.55, 159.50, 270.25, 479.85, 909.75, 1738.15, 3511.40, 7138.10]

# Memoire multiplication fixed size
mem_mult_classic_fixedSize = [1.60, 4.80, 10.40, 20.80, 40.80, 80.00, 157.60, 234.40, 388.00, 695.20, 1309.60, 2538.40, 4996.00, 922337203685477632.00, 922337203685477632.00, 922337203685477632.00, 922337203685477632.00]
mem_mult_karatsuba_fixedSize = [1.60, 17.60, 20.00, 25.60, 34.40, 53.60, 92.00, 168.80, 322.40, 629.60, 1244.00, 2472.80, 4930.40, 9845.60, 19676.00, 39336.80, 78658.40]
mem_mult_choice_fixedSize = [1.60, 3.20, 5.60, 10.40, 20.00, 39.20, 77.60, 12458.40, 13022.40, 14148.00, 16400.80, 21726.40, 32376.80, 53676.80, 96276.00, 181473.60, 351868.00]

# Temps multiplication fixed degree
tps_mult_classic_fixedDeg = [1.45, 2.90, 4.50, 5.95, 7.35, 15.20, 18.90, 31.45, 42.25, 58.65, 96.10, 178.45, 425.55, 1139.55, 3276.40, 8650.95, 22013.70]
tps_mult_karatsuba_fixedDeg = [10.50, 22.65, 33.00, 43.50, 54.70, 66.65, 87.95, 123.60, 163.75, 213.25, 282.80, 407.40, 726.00, 1479.70, 3481.10, 8847.00, 22486.10]
tps_mult_choice_fixedDeg = [2.40, 4.35, 6.35, 8.25, 10.25, 12.60, 17.80, 24.85, 35.90, 63.65, 107.25, 200.00, 390.45, 821.15, 1652.95, 3362.60, 6760.20]

# Memoire multiplication fixed degree
mem_mult_classic_fixedDeg = [13.60, 13.60, 13.60, 13.60, 13.60, 11741.60, 11741.60, 11781.60, 11866.40, 12078.40, 12599.20, 14353.60, 18752.00, 21496.80, 24769.60, 31227.20, 44199.20]
mem_mult_karatsuba_fixedDeg = [13.60, 13.60, 13.60, 13.60, 13.60, 13.60, 13.60, 128.80, 355.20, 804.80, 1636.80, 2767011611056436224.00, 12912720851596689408.00, 21213755684765990912.00, 30437127721620766720.00, 38738162554790068224.00, 47039197387959386112.00]
mem_mult_choice_fixedDeg = [13.60, 13.60, 13.60, 13.60, 13.60, 13.60, 42.40, 75.20, 134.40, 247.20, 477.60, 899.20, 3689348814741912064.00, 11990383647911211008.00, 20291418481080516608.00, 28592453314249830400.00, 35971150943733678080.00]


plt.figure(figsize=(10, 5))
plt.plot(evolution, tps_mult_classic_fixedSize, label="Classic multiplication", marker='o')
plt.plot(evolution, tps_mult_karatsuba_fixedSize, label="Karatsuba algorithm", marker='o')
plt.plot(evolution, tps_mult_choice_fixedSize, label="Adapted algorithm", marker='o')
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Degree size")
plt.ylabel("Execution time (s)")
plt.title("Polynomial multiplication execution time with fixed coeff size=16")
plt.legend()
plt.grid(True)
plt.show()

plt.figure(figsize=(10, 5))
plt.plot(evolution, tps_mult_classic_fixedDeg, label="Classic multiplication", marker='o')
plt.plot(evolution, tps_mult_karatsuba_fixedDeg, label="Karatsuba algorithm", marker='o')
plt.plot(evolution, tps_mult_choice_fixedDeg, label="Adapted algorithm", marker='o')
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Coefficient size (bit)")
plt.ylabel("Execution time (s)")
plt.title("Polynomial multiplication execution time with fixed degree=16")
plt.legend()
plt.grid(True)
plt.show()


plt.figure(figsize=(10, 5))
plt.plot(evolution, mem_mult_classic_fixedSize, label="Classic multiplication", marker='o')
plt.plot(evolution, mem_mult_karatsuba_fixedSize, label="Karatsuba algorithm", marker='o')
plt.plot(evolution, mem_mult_choice_fixedSize, label="Adapted algorithm", marker='o')
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Degree size")
plt.ylabel("Memory used (byte)")
plt.title("Polynomial multiplication memory used with fixed coeff size=16")
plt.legend()
plt.grid(True)
plt.show()

plt.figure(figsize=(10, 5))
plt.plot(evolution, mem_mult_classic_fixedDeg, label="Classic multiplication", marker='o')
plt.plot(evolution, mem_mult_karatsuba_fixedDeg, label="Karatsuba algorithm", marker='o')
plt.plot(evolution, mem_mult_choice_fixedDeg, label="Adapted algorithm", marker='o')
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Degree size")
plt.ylabel("Memory used (byte)")
plt.title("Polynomial multiplication memory used with fixed degree=16")
plt.legend()
plt.grid(True)
plt.show()