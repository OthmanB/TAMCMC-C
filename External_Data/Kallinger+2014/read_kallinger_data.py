import csv
import numpy as np
import matplotlib.pyplot as plt

def read_csv_file(filename):
    data = []
    with open(filename, 'r') as f:
        csv_reader = csv.DictReader(f)
        for row in csv_reader:
            data.append({key: float(value) if value.replace('.', '', 1).isdigit() else value for key, value in row.items()})
    return data

def a_fct(numax, M):
    t=-0.26
    s=-0.613
    k=3710
    return k * numax**s * M**t

filename='data/APOKASC.summary.csv'
data = read_csv_file(filename)

numax = np.asarray([row['numax'] for row in data], dtype=float)
a1 = np.asarray([row['a1'] for row in data], dtype=float)
a2 = np.asarray([row['a2'] for row in data], dtype=float)

M=[0.8, 1, 1.2, 1.4, 1.6]
a_th=[a_fct(numax, M[0])]
a_th.append(a_fct(numax, M[1]))
a_th.append(a_fct(numax, M[2]))
a_th.append(a_fct(numax, M[3]))
a_th.append(a_fct(numax, M[4]))

fig, ax= plt.subplots()
ax.loglog(numax, a1, marker='o', label=r'$a_1$', linestyle='none')
ax.loglog(numax, a2, marker='o', label=r'$a_2$', linestyle='none')
ax.loglog(numax, a_th[0], label=r'$a_{th}$ (0.8 $M_\odot$)')
ax.loglog(numax, a_th[1], label=r'$a_{th}$ ($1.0 M_\odot$)')
ax.loglog(numax, a_th[2], label=r'$a_{th}$ (1.2 $M_\odot$)')
ax.loglog(numax, a_th[3], label=r'$a_{th}$ (1.4 $M_\odot$)')
ax.loglog(numax, a_th[4], label=r'$a_{th}$ (1.6 $M_\odot$)')
ax.set_xlabel(r'$\nu_{max}$ ($\mu$Hz)')
ax.set_ylabel(r'$a_j$ (ppm)')
ax.legend()
fig.savefig('plots/a_vs_numax.jpg', dpi=300)


fig, ax= plt.subplots()
#delta=(a1-a2)/a1
i=2
delta=(a1-a2)/a_th[i]
ax.plot(numax, delta, marker='o', label=r'$(a_1 - a_2)/a(M=$'+ str(M[i]) + r')', linestyle='none')
ax.plot(numax, np.mean(delta)*np.ones(len(numax)), label=r'$\bar{\delta}$')
ax.plot(numax, (np.mean(delta) + np.std(delta))*np.ones(len(numax)), label=r'$\bar{\delta}$')
ax.plot(numax, (np.mean(delta) - np.std(delta))*np.ones(len(numax)), label=r'$\bar{\delta}$')
ax.set_xlabel(r'$\nu_{max}$ ($\mu$Hz)')
ax.set_ylabel(r'$(a_1 - a_2)/a_1$')
ax.legend()
fig.savefig('plots/delta_a_vs_numax.jpg', dpi=300)
