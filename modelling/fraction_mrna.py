import matplotlib.pyplot as plt
from numpy.typing import _256Bit
from models import model_prokaryotic_all
from standard_parameters import standard_parameters


fig, ax = plt.subplots()

parameters, constants, initial_conditions = standard_parameters(
    dna_conc=5*10**-3, vit_conc=1)
results = model_prokaryotic_all(parameters, constants, initial_conditions)

parameters, constants, initial_conditions = standard_parameters(
    dna_conc=5*10**-4, vit_conc=0.5)
results2 = model_prokaryotic_all(parameters, constants, initial_conditions)

time2 = results2[0]
e2 = results2[10]
s2 = results2[11]

time = results[0]
e = results[10]
s = results[11]
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
# ax.plot(time, s, label="enzyme")
# ax.plot(time2, s2, label="enzyme2")
ax.legend()
plt.show()
