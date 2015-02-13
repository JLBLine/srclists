#!/usr/bin/python
'''Script to plot srclists over specified RA,Dec ranges over the top'''

from astropy.wcs import WCS
from wcsaxes import WCSAxes
from astropy import units as units
from numpy import *
import matplotlib.pyplot as plt
import optparse

parser = optparse.OptionParser()

parser.add_option('-s',"--srclist", default=-1, help="srclist used to peel")

parser.add_option('-c','--coords',default='0,-27.0', help ="Enter centre of field RA,Dec in deg,deg. Default=0,-27")

options, args = parser.parse_args()

RA,DEC = map(float, options.coords.split(','))

##A fits image header with which to create a wcs with
header = { 'NAXIS'  : 2,             ##Number of data axis
'NAXIS1' : 10,                  ##Length of X axis
'CTYPE1' : 'RA---SIN',           ##Projection type of X axis
'CRVAL1' : RA,        ##Central X world coord value
'CRPIX1' : 5,                    ##Central X Pixel value
'CUNIT1' : 'deg',                ##Unit of X axes
'CDELT1' : -1*cos(RA*(pi/180.0)),              ##Size of pixel in world co-ord
'NAXIS2' : 10,                  ##Length of X axis
'CTYPE2' : 'DEC--SIN',           ##Projection along Y axis
'CRVAL2' : DEC,                   ##Central Y world coord value
'CRPIX2' : 5,                    ##Central Y Pixel value
'CUNIT2' : 'deg',                ##Unit of Y world coord
'CDELT2' : +1      		     ##Size of pixel in deg
} 

fig = plt.figure(figsize=(15,15))
#ax = fig.add_subplot(111)

wcs = WCS(header=header)
ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=wcs)
fig.add_axes(ax)
tr_fk5 = ax.get_transform("fk5")

class rts_calibrator:
	 def __init__(self):
		 self.name = -1
		 self.rtsflux = -1
		 self.srcflux = -1
		 self.srcra = -1
		 self.srcdec = -1


rts_sources=[]
rts_components=[]
if options.srclist!=-1:
	srclist = open(options.srclist).read().split('ENDSOURCE')
	del srclist[-1]
	for entry in srclist:
		
		lines=entry.split('\n')
		if lines[0]=='':
			del lines[0]
		line=lines[0].split()
		source=rts_calibrator()
		source.name=line[1]
		source.srcra=float(line[2])*15.0
		source.srcdec=float(line[3])
		#if (r1<=source.srcra<=r2) or (r3<=source.srcra<=r4):
			#if dcut[0]<=source.srcdec<=dcut[1]:
					
		rts_sources.append(source)
		for line in lines:
			if 'COMPONENT' in line and 'ENDCOMPONENT' not in line:
				comp=rts_calibrator()
				line=line.split()
				comp.srcra=float(line[1])*15.0
				comp.srcdec=float(line[2])
				rts_components.append(comp)
						

source_ras = [source.srcra for source in rts_sources]
source_decs = [source.srcdec for source in rts_sources]

comp_ras = [source.srcra for source in rts_components]
comp_decs = [source.srcdec for source in rts_components]

ax.plot(comp_ras,comp_decs,'o',color='#FFFF00',transform=tr_fk5,label="COMPONENT",markersize=10,markeredgecolor='k',markeredgewidth=2.0)
ax.plot(source_ras,source_decs,'*',color='c',transform=tr_fk5,label="BASE SOURCE",markersize=30,markeredgecolor='k',markeredgewidth=1.0)

ax.set_title(options.srclist)

#ax.set_xlim(wcs.wcs_world2pix(30,-27.0,0)[0],wcs.wcs_world2pix(-30,-27.0,0)[0])
#ax.set_ylim(wcs.wcs_world2pix(0,-45.0,0)[1],wcs.wcs_world2pix(0,-9,0)[1])
		
ra_ax = ax.coords[0]
dec_ax = ax.coords[1]
ra_ax.set_axislabel('RAJ2000')
dec_ax.set_axislabel('DECJ2000')
ra_ax.set_major_formatter('hh:mm:ss')#,font='Computer Modern Typewriter')
dec_ax.set_major_formatter('dd:mm:ss')

ra_ax.set_ticks(spacing=600 * units.arcmin)
dec_ax.set_ticks(spacing=600 * units.arcmin)

g = ax.coords.grid(linewidth=1.0,alpha=1.0)

plt.show()