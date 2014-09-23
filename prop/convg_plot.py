
import numpy as np
import matplotlib.pyplot as plt

### Prob def

# * 400x400 total points
# * $N_t$ = 200
# * 4x4 regions
# * 16 waveform iterations

raw_times = [210.36, 110.96, 56.34, 28.83]
core_mult = [1, 2, 4, 8]

errors = [.130607, .1230728, .0208692, .030100126, .00595161,
          .00736622524755, .00112752678098, .00044018865,
          .00037872, .000589796, .00011843557,
          4.896327e-5, 4.596327e-5, 6.890354e-5,
          9.085948e-6, 9.4602816e-6]

# for t, m in zip(raw_times, core_mult):
#     print m, round(raw_times[0]/t, 2), round(raw_times[0]/t/m, 2)

fig = plt.figure()#figsize=(6,5))
ax = fig.add_subplot(111)
ax.plot(range(1, len(errors)+1), np.log(errors), color='k')
ax.set_xlim([1, len(errors)+1])
ax.set_ylim([-12, -1])
ax.set_xlabel("Waveform iteration", fontsize=18)
ax.set_ylabel(r"$\ln(||u^{[k]}(T)-u(T)||_\infty)$", fontsize=22)

fig.savefig("convg_plot.png")

