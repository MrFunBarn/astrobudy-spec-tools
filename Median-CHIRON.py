#! /usr/bin/python
import numpy as np
import pyfits
import matplotlib.pyplot as plt
import sys 
import os 
import datetime                                                           

#     Median.py simply median combines the flux array of .fits data from the
# CHIRON IDL pipeline and saves the resulting data to a new file.
if len(sys.argv) >= 3:
    exps = []
    for exp in sys.argv: exps.append(exp)
else:
    print('Not Enough Exposure, check that each exposure name is separated by \
          a space.')

# Define the input files and read in each spectra and it's header with pyfits.
if len(sys.argv) == 3:
    specfile1 = sys.argv[1]                                             
    specfile2 = sys.argv[2]
    spec1, header1 = pyfits.getdata(specfile1, 0, header=True)                     
    spec2, header2 = pyfits.getdata(specfile2, 0, header=True)
    medianflux = np.array([spec1[:,:,1],spec2[:,:,1]])                                
    medianblank = np.zeros( (62, 801, 2) )
    medianflux = np.median(medianflux, axis=0)
    medianblank[:,:,1] = np.copy(medianflux[:,:])
    medianblank[:,:,0] = np.copy(spec1[:,:,0])
    medianspec = medianblank

elif len(sys.argv) == 4:
    specfile1 = sys.argv[1]                                             
    specfile2 = sys.argv[2]
    specfile3 = sys.argv[3]
    spec1, header1 = pyfits.getdata(specfile1, 0, header=True)                     
    spec2, header2 = pyfits.getdata(specfile2, 0, header=True)
    spec3, header3 = pyfits.getdata(specfile3, 0, header=True)
    medianflux = np.array([spec1[:,:,1],spec2[:,:,1],spec3[:,:,1]])                                
    medianblank = np.zeros( (62, 801, 2) )
    medianflux = np.median(medianflux, axis=0)
    medianblank[:,:,1] = np.copy(medianflux[:,:])
    medianblank[:,:,0] = np.copy(spec1[:,:,0])
    medianspec = medianblank

elif len(sys.argv) == 5:
    specfile1 = sys.argv[1]                                             
    specfile2 = sys.argv[2]
    specfile3 = sys.argv[3]
    specfile4 = sys.argv[4]
    spec1, header1 = pyfits.getdata(specfile1, 0, header=True)                     
    spec2, header2 = pyfits.getdata(specfile2, 0, header=True)
    spec3, header3 = pyfits.getdata(specfile3, 0, header=True)
    spec4, header4 = pyfits.getdata(specfile4, 0, header=True)
    medianflux = np.array([spec1[:,:,1],spec2[:,:,1],spec3[:,:,1] \
            ,spec4[:,:,1]])                                
    medianblank = np.zeros( (62, 801, 2) )
    medianflux = np.median(medianflux, axis=0)
    medianblank[:,:,1] = np.copy(medianflux[:,:])
    medianblank[:,:,0] = np.copy(spec1[:,:,0])
    medianspec = medianblank

elif len(sys.argv) == 6:
    specfile1 = sys.argv[1]                                             
    specfile2 = sys.argv[2]
    specfile3 = sys.argv[3]
    specfile4 = sys.argv[4]
    specfile5 = sys.argv[5]
    spec1, header1 = pyfits.getdata(specfile1, 0, header=True)                     
    spec2, header2 = pyfits.getdata(specfile2, 0, header=True)
    spec3, header3 = pyfits.getdata(specfile3, 0, header=True)
    spec4, header4 = pyfits.getdata(specfile4, 0, header=True)
    spec5, header5 = pyfits.getdata(specfile5, 0, header=True)
    medianflux = np.array([spec1[:,:,1],spec2[:,:,1],spec3[:,:,1] \
            ,spec4[:,:,1],spec5[:,:,1]])                                
    medianblank = np.zeros( (62, 801, 2) )
    medianflux = np.median(medianflux, axis=0)
    medianblank[:,:,1] = np.copy(medianflux[:,:])
    medianblank[:,:,0] = np.copy(spec1[:,:,0])
    medianspec = medianblank

# Make a file name for the median combined data automatically using the object
# name and UT date of exposure. USES UT DATE of exposure start.
objectn = header1['object']
objectn = objectn.split('-')[0] #in case of combined name and date.
objectd = header1['UTSHUT']
objectd = objectd.split('T')[0] #Keeps only date, not the time of exp.
newfilename = objectn+'-'+objectd+'-Median.fits'

# Print the exposres being median combined and how much of each constiuant is
# in the final exposure.
print('Median Combining Exposures')
for ex in exps[1:]:
    print(ex)

# Build the new fits file from the data and header of 1st exposure and then
# write the file to the newfilename from above. The output_verify='ignore' is
# important because there are non-standard Keys in the CHIRON fits header.Also
# stamps a record of the modifications and mod date/time into the header.
hdulist = pyfits.open(specfile1)
data = hdulist[0].data
prihdr = hdulist[0].header
now = datetime.datetime.now()
prihdr['history'] = 'Median Combined on ' \
                    +now.strftime("%Y-%m-%d %H:%M") \
                    +' From '+str(exps[1:]).strip("[]").replace(',',' ')+ ' Header is pulled from First Exposure.'
hdulist[0].data = medianspec
hdulist.writeto(newfilename, output_verify='ignore')
hdulist.close()


