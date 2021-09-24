import matplotlib.pyplot as plt
from numpy.typing import _256Bit
from models import model_prokaryotic_all
from standard_parameters import standard_parameters


fig, ax = plt.subplots()

parameters, constants, initial_conditions = standard_parameters(
    dna_conc=5*10**-3, vit_conc=2)
results2 = model_prokaryotic_all(parameters, constants, initial_conditions)

parameters, constants, initial_conditions = standard_parameters(
    dna_conc=5*10**-3, vit_conc=0.5)
results = model_prokaryotic_all(parameters, constants, initial_conditions)

time2 = results2[0]
e2 = results2[10]
s2 = results2[11]
p2 = results2[12]

umrna2 = results2[2]
umrna_vit2 = results2[6]
cmrna2 = results2[3]
total = umrna2 + umrna_vit2 + cmrna2

frac_umrna2 = umrna2 / total
frac_umrna_vit2 = umrna_vit2 / total
frac_cmrna2 = cmrna2 / total


time = results[0]
e = results[10]
s = results[11]
p = results[12]
umrna = results[2]
umrna_vit = results[6]
cmrna = results[3]
total = umrna + umrna_vit + cmrna

frac_umrna = umrna / total
frac_umrna_vit = umrna_vit / total
frac_cmrna = cmrna / total

ax.plot(time, frac_umrna, label="umrna")
ax.plot(time, frac_umrna_vit, label="umrna_vit")
ax.plot(time, frac_cmrna, label="cmrna")


ax.plot(time, frac_umrna2, label="umrna2")
ax.plot(time, frac_umrna_vit2, label="umrna_vit2")
ax.plot(time, frac_cmrna2, label="cmrna2")
# ax.plot(time, s, label="enzyme")
# ax.plot(time2, s2, label="enzyme2")
ax.legend()

fig2, ax2 = plt.subplots()
ax2.plot(time, p, label="p")
ax2.plot(time, p2, label="p2")
plt.show()
