#!/usr/bin/python
import sys
import os
import datetime
import math


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

import pyfits # From the space telescope science institute.

import InteractiveNorm
# Normalize-CHIRON version 2.0.3 -Beta
# Written by Brandon Bell April 2013
#
# Normalization-CHIRON calls an IDL routine to interactivly normalize each
# order of a file and then dumps all the fit data, normalized data, and plots
# of the fit into .tgz archives. It then makes a .fits of the data with
# apprpriate headers

# When re-normalizing the median file comes first then the data be
# re-normalized you don't have to overwrite it if you don't want.
renorm = False
filename = sys.argv[1]
if len(sys.argv) == 3:
    renormfile = sys.argv[2]
    renorm = True
print('Normalize-CHIRON Version 2.0.3 -Beta\n')
print("""    Normalize-CHIRON will normalize the orders you want, and then make tarballs
containing two columb text files of the normalized data and fit function at
each order, plots of the fit function at each order, and finally a whole
spectrum plot and fits file containing the new data. Final output should
produce only 5 files.\n""")
if renorm == True:
    print('  -> Enter the orders to re-normalize' \
          ' separated by commas and hit enter.')
else:
    print('  -> Enter the orders to normalize' \
          ' separated by commas and hit enter.')
print('  -> To normalize all orders just hit enter.')
print("  -> Order #'s begin at 0 not 1.")
print("  -> Order #'s don't need to be in order, and spaces are ok.")
print('  -> No trailing commas please.\n')

# Get the list of orders to normalize or set up to do all the orders.
try:
    inp = input("Orders: ")
    try:
        inpl=list(inp)
    except:
        inps=str(inp)
    if len(str(inp)) <= 2:
        inps = inps.split(',')[0]
        orders = int(inps)
        orders = [orders]
    else:
        o = [int(x) for x in inpl]
        orders = o
except:
    orders = range(62)

# Read in the fits file(s) and retirve header info for latter use.
spec, header = pyfits.getdata(filename, 0, header=True)
if renorm == True:
    respec, reheader = pyfits.getdata(renormfile, 0, header=True)
objectn = header['object']
objectn = objectn.split('-')[0] #in case of combined name and date.
objectd = header['UTSHUT']
objectd = objectd.split('T')[0] #Keeps only date, not the time of exp.
root = filename.split('/')[0]

# Variables for the whole spec plot tick labeling.
majorLocator   = MultipleLocator(100)
minorLocator   = MultipleLocator(25)

# Define the new array to be filled with the normalized values.
normspec = np.copy(spec)

def norm(idl):
    "norm calls idl on the array of flux values to interactivly normalize it."
    # Run the IDL routine to normalize the function.
    cmd = 'idl "specnorm-norm.pro" -quiet -args ' + idl
    os.system(cmd)
    return

def yesnocheck(yn):
    if yn =='y':
        return yn
    elif yn=='n':
        return yn
    else:
        print('Please enter yes (y) or no (n)')
        yn = raw_input('(y/n): ')
        yesnocheck(yn)
    return yn

def smooth3(x):
   global sm
   nelm = np.shape(x)[1]
   sm = np.zeros( (62,int(math.floor(nelm/3)),2) )
   for o in range(62):
       for i in range(int(math.floor(nelm/3))):
          n=3*i
          sm[o,i,0] = ( x[o,n,0] + x[o,n+1,0] + x[o,n+2,0] ) / 3
          sm[o,i,1] = ( x[o,n,1] + x[o,n+1,1] + x[o,n+2,1] ) / 3
   return sm;

for order in orders:
    ro = str(int(order))
    # Make array to pass to IDL for normalization.
    idlspec = np.zeros( (801,2) )
    idlspec[:,0] = normspec[order,:,0]
    idlspec[:,1] = normspec[order,:,1]

    # Create nessesary file names for data transfer file and the text file
    # that will hold normalized dat for saving and write the IDL transfer file.
    normspecname = objectn+'-'+objectd+'-'+'Order'+ro+'-Normalized-Data.dat'
    fitfilename = objectn+'-'+objectd+'-'+'Order'+ro+'-Fit-Function.dat'
    idlpassfile = objectn+'-'+objectd+'-'+ro+'.dat'
    normfilename = "norm-"+idlpassfile
    np.savetxt(idlpassfile, idlspec)

    norm(idlpassfile)

    # Read in the fit function .dat file from IDL. The filename is defined
    # in the IDL procedure specnorm-norm.pro written by Brandon.
    try:     # If IDL fails to wrote the file try again.
        fit = np.genfromtxt(normfilename)
    except:
        norm(idlpassfile)

    # Normalize the flux array at this order.
    normspec[order,:,1] = normspec[order,:,1] / fit[:,1]

    # Make a line array for plotting on the normaliztion fit plots.
    redline = np.copy(normspec[order,:,:])
    redline[:,1] = 1

    # Make two colum array of normalized flux to write into a text file, then
    # write the normalized flux and the fit function to .dat files.
    x = np.copy(idlspec)
    x[:,1] = normspec[order,:,1]
    np.savetxt(normspecname, x)
               #header="# The normalzed data in the form (wavelength, flux)")
    np.savetxt(fitfilename, fit)
               #header='# The normalization fit function (wavelength, flux)')

    # Remove the .dat files used to pass the data to and from IDL.
    cmd = 'rm ' + idlpassfile
    os.system(cmd)
    cmd = 'rm ' + normfilename
    os.system(cmd)

    # Make a plot of the fit and un-normalized data to check the fit at latter
    # times if something seems fishy in the data.
    xmin = normspec[order,0,0]
    xmax = normspec[order,800,0]
    ymax = np.amax(spec[order,:,1])
    fig=plt.figure(1)
    fig.suptitle(objectn+' '+objectd+' Order '+ro)
    plt.subplot(211)
    plt.ylabel('Flux')
    plt.title('Data and Fit')
    plt.plot(spec[order,:,0],spec[order,:,1])
    plt.plot(fit[:,0],fit[:,1])
    plt.ylim(0,(1.25*ymax))
    plt.xlim(xmin, xmax)
    plt.subplot(212)
    plt.title('Normalized Spectra')
    plt.ylabel('Normalized Flux')
    plt.xlabel(r'Wavelength $(\AA)$')
    plt.plot(normspec[order,:,0],normspec[order,:,1])
    plt.plot(redline[:,0], redline[:,1])
    plt.ylim(0.5, 1.5)
    plt.xlim(xmin, xmax)
    pltfilename = objectn+'-'+objectd+'Reference-Data-Fit-'+'Order'+ro+'.png'
    plt.savefig(pltfilename, dpi=300)
    plt.clf()

     # Replace the normalized data if your renormalizing.
    if renorm == True:
        respec[order,:,1] = np.copy(normspec[order,:,1])

    # Clear terminal and print current order for sanities sake.
    os.system('clear')
    print('\n')
    print('  Finished Order '+ro+'  :)\n')
    try:
        zc = int(orders.index(order) +1)
        print('  Working on Order '+str(orders[zc])+' now.\n')
    except IndexError:
        print('All specified orders completed\n')

# Check for normalized file, if found ask if you want to overwrite then,
# write out the final nomalized data into a new fits file and add the
# modifications made and time of modifications. Also write the header to a
# text file for inclusion with the data text files.
newfilename = root+'/'+objectn+'-'+objectd+'-Normalized.fits'
skip = False
try:
    if renorm == True:
        normspec = np.copy(respec)
    with open(newfilename): pass
    print(newfilename+' Allready exists')
    yn = raw_input('Overwrite Files? (y/n):')
    yesnocheck(yn)
    if yn=='y':
        print('Writing .fits File')
        hdulist = pyfits.open(filename)
        data = hdulist[0].data
        prihdr = hdulist[0].header
        now = datetime.datetime.now()

        # Add renorm data to header if needed.
        if renorm == True:
            hdu = pyfits.open(renormfile)
            hdu[0].header['history'] = 'Orders '+str(orders).strip('[]')+ \
                                       ' renormalized on ' \
                                       +now.strftime("%Y-%m-%d %H:%M")
            hdu[0].header.totextfile(objectn+'-'+objectd+'-Header.dat',
                                     clobber=True)
            # Make HDUd and write the files.
            hdu[0].data = normspec
            hdu.writeto(newfilename, output_verify='ignore', clobber=True)
            hdu.close()
        else:
            prihdr['history'] = 'All orders normalized on ' \
                                +now.strftime("%Y-%m-%d %H:%M")
            hdulist[0].header.totextfile(objectn+'-'+objectd+'-Header.dat',
                                         clobber=True)
            # Make HDUd and write the files.
            hdulist[0].data = normspec
            hdulist.writeto(newfilename, output_verify='ignore', clobber=True)
            hdulist.close()

        # Define file and directory names for file operations.
        fitplot = root+'/'+objectn+'-'+objectd+'Normaliztion-Fit-Plots'
        datatxt = root+'/'+objectn+'-'+objectd+'Normalized-Data-Text-Files'
        fitdata = root+'/'+objectn+'-'+objectd+'Fit-Function-Text-Files'

        print('Moving Files')
        os.system('mv *Order*.png '+fitplot)
        os.system('mv *Normalized-Data.dat *Header.dat '+datatxt)
        os.system('mv *Fit-Function.dat '+fitdata)

        print('Smoothing')
        smooth3(normspec)

        # Draw a plot of the entire spectrum for referance.
        print('Generating Whole Spectrum Plot')
        xmin = normspec[0,0,0]
        xmax = normspec[61,800,0]
        ymax = 1.5
        fig = plt.figure(figsize=(16, 2))
        ax = fig.add_subplot(1,1,1,autoscale_on=False,
                             ylim=(0,ymax),xlim=(xmin, xmax))
        plt.title(objectn+' '+objectd, fontsize=9)
        plt.xlabel(r'Wavelength ($\AA$)',fontsize=8)
        plt.ylabel('Flux', fontsize=8)
        ax.tick_params(axis='both', which='major', labelsize=6)
        ax.xaxis.set_minor_locator(minorLocator)
        ax.xaxis.set_major_locator(majorLocator)
        for order in range(62):
            ro = order
            if order%2==0:    # Check wiether order is odd or even.
                ax.plot(sm[order,:,0],sm[order,:,1],
                       color='green', linewidth=.07)
                ax.text(normspec[order,400,0], .9, ro,
                        fontsize=4, color='green',
                        horizontalalignment='center',
                        va='center')
            else:
                ax.plot(sm[order,:,0],sm[order,:,1],
                        color='blue', linewidth=.07)
                ax.text(normspec[order,400,0], 1.1, ro,
                        fontsize=4, color='blue',
                        horizontalalignment='center',
                       va='center')
        pltfilename = root+'/'+objectn+'-'+objectd+'-Whole-Spectra.png'
        plt.savefig(pltfilename, bbox_inches='tight', dpi=600)
        print('Finished Normalization, Have Fun!')

        skip = True
    else:
        # Preform all the above operations with a new filename.
        print('Type new filename addition without .fits')
        print('Example: '+objectn+'-'+objectd+'-Normalized-blablabla')
        newfilext = raw_input(objectn+'-'+objectd+'-Normalized-')
        newfilename = root+'/'+objectn+'-'+objectd+'-Normalized-' \
                    +newfilext+'.fits'

        fitplot = root+'/'+objectn+'-'+objectd+ \
                  'Normaliztion-Fit-Plots-'+newfilext
        datatxt = root+'/'+objectn+'-'+objectd+ \
                  'Normalized-Data-Text-Files-'+newfilext
        fitdata = root+'/'+objectn+'-'+objectd+ \
                  'Fit-Function-Text-Files-'+newfilext

        print('Writing .fits File')
        hdulist = pyfits.open(filename)
        data = hdulist[0].data
        prihdr = hdulist[0].header
        now = datetime.datetime.now()
        if renorm == True:
            hdu = pyfits.open(renormfile)
            hdu[0].header['history'] = 'Orders '+str(orders).strip('[]')+ \
                                       ' renormalized on ' \
                                       +now.strftime("%Y-%m-%d %H:%M")
            hdu[0].header.totextfile(objectn+'-'+objectd+'-Header-'
                                     +newfilext+'.dat', clobber=True)
            # Make HDUd and write the files.
            hdu[0].data = normspec
            hdu.writeto(newfilename, output_verify='ignore', clobber=True)
            hdu.close()
        else:
            prihdr['history'] = 'All orders normalized on ' \
                                +now.strftime("%Y-%m-%d %H:%M")
            hdulist[0].header.totextfile(objectn+'-'+objectd+
                                         '-Header-'+newfilext+'.dat')
            hdulist[0].data = normspec
            hdulist.writeto(newfilename, output_verify='ignore')
            hdulist.close()

        print('Building Directory Tree and Moving Files!')

        os.system('mkdir '+fitplot+' && mv *Order*.png '+fitplot)

        os.system('mkdir '+datatxt+
                  ' && mv *Normalized-Data.dat *Header-'
                  +newfilext+'.dat '+datatxt)

        os.system('mkdir '+fitdata+' && mv *Fit-Function.dat '+fitdata)

        print('Smoothing')
        smooth3(normspec)

        # Draw a plot of the entire spectrum for referance.
        print('Generating Whole Spectrum Plot')
        xmin = normspec[0,0,0]
        xmax = normspec[61,800,0]
        ymax = 1.5
        fig = plt.figure(figsize=(16, 2))
        ax = fig.add_subplot(1,1,1,autoscale_on=False,
                             ylim=(0,ymax),xlim=(xmin, xmax))
        plt.title(objectn+' '+objectd, fontsize=9)
        plt.xlabel(r'Wavelength ($\AA$)',fontsize=8)
        plt.ylabel('Flux', fontsize=8)
        ax.tick_params(axis='both', which='major', labelsize=6)
        ax.xaxis.set_minor_locator(minorLocator)
        ax.xaxis.set_major_locator(majorLocator)
        for order in range(62):
            ro = order
            if order%2==0:
                ax.plot(sm[order,:,0],sm[order,:,1],
                       color='green', linewidth=.07)
                ax.text(normspec[order,400,0], .9, ro,
                        fontsize=4, color='green',
                        horizontalalignment='center',
                        va='center')
            else:
                ax.plot(sm[order,:,0],sm[order,:,1],
                        color='blue', linewidth=.07)
                ax.text(normspec[order,400,0], 1.1, ro,
                        fontsize=4, color='blue',
                        horizontalalignment='center',
                       va='center')
        pltfilename = root+'/'+objectn+'-'+objectd+ \
                      '-Whole-Spectra-'+newfilext+'.png'
        plt.savefig(pltfilename, bbox_inches='tight', dpi=600)
        print('Finished Normalization, Have Fun!')
        skip = True
except IOError:
    print('Writing .fits File')
    hdulist = pyfits.open(filename)
    data = hdulist[0].data
    prihdr = hdulist[0].header
    now = datetime.datetime.now()
    prihdr['history'] = 'All orders normalized on ' \
                        +now.strftime("%Y-%m-%d %H:%M") \
                        +' By Brandon'
    hdulist[0].header.totextfile(objectn+'-'+objectd+'-Header.dat')
    hdulist[0].data = normspec
    hdulist.writeto(newfilename, output_verify='ignore')
    hdulist.close()

# Create .tgz archive of all the fit and data .dat files and the normaliztion
# referance images and remove all the .dat and .png files outside of archive.
if skip == False:
    print('Building Directory Tree and Moving Files!')
    fitplot = root+'/'+objectn+'-'+objectd+'Normaliztion-Fit-Plots'
    os.system('mkdir '+fitplot+' && mv *Order*.png '+fitplot)
    datatxt = root+'/'+objectn+'-'+objectd+'Normalized-Data-Text-Files'
    os.system('mkdir '+datatxt+
              ' && mv *Normalized-Data.dat *Header.dat '+datatxt)
    fitdata = root+'/'+objectn+'-'+objectd+'Fit-Function-Text-Files'
    os.system('mkdir '+fitdata+' && mv *Fit-Function.dat '+fitdata)

ct = raw_input('Create Tarballs (y/n)? ')
yesnocheck(ct)
if ct == 'y':
    print('Creating Tarballs')
    os.system('tar -C '+root+" -czf"+fitplot+
              '.tgz '+fitplot.split('/')[1])
    os.system('tar -C '+root+" -czf "+datatxt+
              '.tgz '+datatxt.split('/')[1])
    os.system('tar -C '+root+" -czf "+fitdata+
              '.tgz '+fitdata.split('/')[1])

# Draw a plot of the entire spectrum for referance.
if skip == False:
    print('Smoothing')
    smooth3(normspec)

    print('Generating Whole Spectrum Plot')
    xmin = normspec[0,0,0]
    xmax = normspec[61,800,0]
    ymax = 1.5
    fig = plt.figure(figsize=(16, 2))
    ax = fig.add_subplot(1,1,1,autoscale_on=False,
                         ylim=(0,ymax),xlim=(xmin, xmax))
    plt.title(objectn+' '+objectd, fontsize=9)
    plt.xlabel(r'Wavelength ($\AA$)',fontsize=8)
    plt.ylabel('Flux', fontsize=8)
    ax.tick_params(axis='both', which='major', labelsize=6)
    ax.xaxis.set_minor_locator(minorLocator)
    ax.xaxis.set_major_locator(majorLocator)
    for order in range(62):
        ro = order
        if order%2==0:
            ax.plot(sm[order,:,0],sm[order,:,1],
                   color='green', linewidth=.07)
            ax.text(normspec[order,400,0], .9, order,
                    fontsize=4, color='green',
                    horizontalalignment='center',
                    va='center')
        else:
            ax.plot(sm[order,:,0],sm[order,:,1],
                    color='blue', linewidth=.07)
            ax.text(normspec[order,400,0], 1.1, ro,
                    fontsize=4, color='blue',
                    horizontalalignment='center',
                    va='center')
    pltfilename = root+'/'+objectn+'-'+objectd+'-Whole-Spectra.png'
    plt.savefig(pltfilename, bbox_inches='tight', dpi=800)

    print('Finished Normalization, Have Fun!')
