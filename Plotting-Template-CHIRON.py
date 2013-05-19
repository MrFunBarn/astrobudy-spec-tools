#! /usr/bin/python
import numpy as np
import pyfits
import matplotlib.pyplot as plt
import sys 
import os

                                                              
order = sys.argv[1]
filename = sys.argv[2]

spec, header = pyfits.getdata(filename, 0, header=True)                     
objectn = header['object']
objectn = objectn.split('-')[0]
objectd = header['UTSHUT']
objectd = objectd.split('T')[0] #Keeps only date, not the time of exp.

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
xmin = spec[order,0,0]                                                
xmax = spec[order,800,0]
ymax = 1.1*(np.amax(spec[order,:,1]))
ymin = (np.amin(spec[order,:,1])-(0.1*(np.amin(spec[order,:,1]))))
ax.set_ylim(0,ymax)
ax.set_xlim(xmin, xmax)


plt.title(objectn)
plt.xlabel(r'Wavelength ($\AA$)')
plt.ylabel('Flux')
ax.text(.84,0.85,'CHIRON 2013 A Data',
        horizontalalignment='center',transform=ax.transAxes)

ax.plot(spec[order,:,0],spec[order,:,1], label=objectd)

# Example Code for a line ID at H alpha (order 40). w = wavelength of line,
# lf = line flux (hight of bottom part of line.), o = length of line.
#w=6563
#lf=np.amax(spec[order,:,1])
#o=.06*np.amax(spec[order,:,1])
#plt.annotate(r'H$\alpha$', 
#             xy=(w,lf), 
#             xytext=(w,(lf + o)), 
#             horizontalalignment='center', 
#             arrowprops=dict(facecolor='black', arrowstyle='-'))

ax.legend()

#plotname = objectn+'-'+objectd+'-'+order+'.png'
#plt.savefig(plotname, dpi=300)
plt.show()
