import numpy as np
import matplotlib.pyplot as plt

import codecs

from pyMie import PyMieScatterForAll

# creating an miePara object
mie = PyMieScatterForAll()

data_out = np.zeros(3)

filecp = codecs.open("./Data_MiePlot/488nm_n_1_7.txt",encoding='cp1252')
data_mieplot = np.loadtxt(filecp,skiprows=14)

# using MiePara's method

lmb = 0.488
radius0 = 0.001
fv = 1.0

n_i_real = 0.050
n_i_imag = 3.019917
n_medium = 1.7

num_cal = 500

Qext = np.zeros(num_cal)
Qabs = np.zeros(num_cal)
Qsca = np.zeros(num_cal)


radius = np.linspace(radius0, radius0*200, num_cal)

for i in range(1,num_cal):
    
    print i

    mie.mieCal(lmb, radius[i], fv, n_i_real, n_i_imag, n_medium, data_out)

    Qext[i] = data_out[0]
    Qabs[i] = data_out[1]
    Qsca[i] = data_out[2]



plt.figure(0)
plt.subplot(2,2,1)
plt.plot(radius*1e3, Qabs, 'y-s',label='Qabs')
plt.plot(data_mieplot[:,0]*1e3, data_mieplot[:,3], 'k--',label='Qabs MiePlot')
plt.xlabel(r'Radius ($\mathrm{\mu m}$)')
plt.legend()

plt.subplot(2,2,2)
plt.plot(radius*1e3, Qsca, 'g-o',label='Qsca')
plt.plot(data_mieplot[:,0]*1e3, data_mieplot[:,2], 'k--',label='Qsca MiePlot')
plt.xlabel(r'Radius ($\mathrm{\mu m}$)')
plt.legend()

# data from Mie plot
plt.subplot(2,2,3)
plt.plot(radius*1e3, Qext, 'r-^',label='Qext')
plt.plot(data_mieplot[:,0]*1e3, data_mieplot[:,1], 'k--',label='Qext MiePlot')
plt.xlabel(r'Radius ($\mathrm{\mu m}$)')
plt.legend()

plt.tight_layout(pad=0.5)
plt.savefig("Fig_compare")
#plt.show()


