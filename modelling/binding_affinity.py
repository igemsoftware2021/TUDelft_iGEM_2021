import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

vit = np.linspace(0, 0.5, 100)
mrna_i = 0.005
k_d = 1


def aptamer_binding(vit, mrna_i, k_d):
    mrna = np.zeros(vit.shape[0])
    mrna_vit = np.zeros(vit.shape[0])
    for ii in range(vit.shape[0]):
        a = 1
        b = - (vit[ii] + mrna_i + k_d)
        c = vit[ii] * mrna_i
        all_roots = np.roots([a, b, c])
        mrna_vit[ii] = np.amin(all_roots)
        mrna[ii] = mrna_i - mrna_vit[ii]
    mrna_ratio = np.divide(mrna, mrna_vit + mrna)  # fraction free
    return mrna_ratio


mrna_frac_free = aptamer_binding(vit, mrna_i, k_d)

# Create the figure and the line
fig, ax = plt.subplots()
line, = plt.plot(vit*1000, mrna_frac_free, lw=2)
ax.set_xlabel('Vitamin concentration (nM)')
ax.set_ylabel("Fraction free")


axcolor = 'lightgoldenrodyellow'
ax.margins(x=0)
plt.subplots_adjust(left=0.25, bottom=0.25)

mrna_slider = Slider(
    ax=plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor),
    label="mrna_i",
    valmin=0.1 * mrna_i,
    valmax=10 * mrna_i,
    valinit=mrna_i,
)


def update(val):
    mrna_i = mrna_slider.val
    mrna_frac_free = aptamer_binding(vit, mrna_i, k_d)
    line.set_ydata(mrna_frac_free)
    fig.canvas.draw_idle()


mrna_slider.on_changed(update)
plt.show()
