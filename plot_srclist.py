#!/usr/bin/python
'''Script to plot srclists with a specified field centre'''

from astropy.wcs import WCS
from wcsaxes import WCSAxes
from astropy import units as units
from numpy import *
import matplotlib.pyplot as plt
import optparse

parser = optparse.OptionParser()

parser.add_option('-s',"--srclist", default=-1, help="srclist used to calibrate or peel")

parser.add_option('-p',"--parent_srclist", default=0, help="srclist used to create srclist for plotting")

parser.add_option('-m',"--metafits", default=0, help="metafits to plot the beam with (optional)")

parser.add_option('-c','--coords',default='0,-27.0', help ="Enter centre of field RA,Dec in deg,deg. Default=0,-27")

parser.add_option('-w','--write',default=False, action='store_true',help ="Switch on to save figure instead of pop up plot")

options, args = parser.parse_args()

RA,DEC = map(float, options.coords.split(','))

def arcdist(RA1,RA2,Dec1,Dec2):
	'''calculates distance between two celestial points in degrees'''
	dr = pi/180.0
	in1 = (90.0 - Dec1)*dr
	in2 = (90.0 - Dec2)*dr
	RA_d = (RA1 - RA2)*dr
	cosalpha = cos(in1)*cos(in2) + sin(in1)*sin(in2)*cos(RA_d)
	alpha = arccos(cosalpha)
	return alpha/dr

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
		 self.extended = 'no'

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
				
parent_sources=[]
if options.parent_srclist!=0:
	srclist = open(options.parent_srclist).read().split('ENDSOURCE')
	del srclist[-1]
	for entry in srclist:
		lines=entry.split('\n')
		if lines[0]=='':
			del lines[0]
		line=lines[0].split()
		source=rts_calibrator()
		source.name=line[1]
		if 'EXT' in line[1]: source.extended='yes'
		source.srcra=float(line[2])*15.0
		source.srcdec=float(line[3])
		#if (r1<=source.srcra<=r2) or (r3<=source.srcra<=r4):
			#if dcut[0]<=source.srcdec<=dcut[1]:
		for line in lines:
			if 'COMPONENT' in line and 'ENDCOMPONENT' not in line:
				source.extended='yes'
		parent_sources.append(source)

comp_ras = [source.srcra for source in rts_components]
comp_decs = [source.srcdec for source in rts_components]
src_ras = [source.srcra for source in rts_sources]
src_decs = [source.srcdec for source in rts_sources]

if options.parent_srclist==0:
	ax.plot(comp_ras,comp_decs,'o',color='#FFFF00',transform=tr_fk5,label="COMPONENT",markersize=10,markeredgecolor='k',markeredgewidth=2.0)
	ax.plot(src_ras,src_decs,'*',color='c',transform=tr_fk5,label="BASE SOURCE",markersize=20,markeredgecolor='k',markeredgewidth=1.0)
	#ax.plot(comp_ras,comp_decs,'o',color='#FFFF00',label="COMPONENT",markersize=10,markeredgecolor='k',markeredgewidth=2.0)
	#ax.plot(src_ras,src_decs,'*',color='c',label="BASE SOURCE",markersize=20,markeredgecolor='k',markeredgewidth=1.0)

else:
	def split_sources(sources):
		for source in sources:
			offsets = [arcdist(source.srcra,par.srcra,source.srcdec,par.srcdec) for par in parent_sources]
			parent_source = parent_sources[offsets.index(min(offsets))]
			source.extended = parent_source.extended
		single_sources = [source for source in sources if source.extended=='no']
		extended_sources = [source for source in sources if source.extended=='yes']
		single_ras = [source.srcra for source in single_sources]
		single_decs = [source.srcdec for source in single_sources]
		extended_ras = [source.srcra for source in extended_sources]
		extended_decs = [source.srcdec for source in extended_sources]
	
		return single_ras,single_decs,extended_ras,extended_decs

	sing_comp_ras, sing_comp_decs, ext_comp_ras, ext_comp_decs = split_sources(rts_components)
	sing_src_ras, sing_src_decs, ext_src_ras, ext_src_decs = split_sources(rts_sources)

	if len(sing_comp_ras)!=0: ax.plot(sing_comp_ras,sing_comp_decs,'o',color='#FFFF00',transform=tr_fk5,label="COMPONENT",markersize=10,markeredgecolor='k',markeredgewidth=2.0)
	if len(ext_comp_ras)!=0: ax.plot(ext_comp_ras,ext_comp_decs,'*',color='r',transform=tr_fk5,label="EXT COMPONENT",markersize=20,markeredgecolor='k',markeredgewidth=1.0)
	if len(sing_src_ras)!=0: ax.plot(sing_src_ras,sing_src_decs,'*',color='c',transform=tr_fk5,label="BASE SOURCE",markersize=30,markeredgecolor='k',markeredgewidth=1.0)
	if len(ext_src_ras)!=0: ax.plot(ext_src_ras,ext_src_decs,'*',color='g',transform=tr_fk5,label="EXT SOURCE",markersize=30,markeredgecolor='k',markeredgewidth=1.0)

ax.set_title(options.srclist)

#ax.set_xlim(wcs.wcs_world2pix(30,-27.0,0)[0],wcs.wcs_world2pix(-30,-27.0,0)[0])
#ax.set_ylim(wcs.wcs_world2pix(0,-45.0,0)[1],wcs.wcs_world2pix(0,-9,0)[1])


if options.metafits!=0:
	print 'check here'
	
	try:
		import astropy.io.fits as pyfits
	except ImportError:
		import pyfits
	from mwapy import ephem_utils
	from mwapy.pb import primary_beam
	from numpy import *
	from mwapy.pb import mwa_tile
	import sys
	try:
		f=pyfits.open(options.metafits)
	except Exception,e:
		print 'Unable to open metafits file %s: %s' % (options.metafits,e)
		sys.exit(1)
	if not 'DELAYS' in f[0].header.keys():
		print 'Cannot find DELAYS in %s' % options.metafits
		sys.exit(1)
	if not 'LST' in f[0].header.keys():
		print 'Cannot find LST in %s' % options.metafits
		sys.exit(1)
	if not 'GPSTIME' in f[0].header.keys():
		print 'Cannot find GPSTIME in %s' % options.metafits
		sys.exit(1)
	if not 'FREQCENT' in f[0].header.keys():
		print 'Cannot find FREQCENT in %s' % options.metafits
		sys.exit(1)
	if not 'RA' in f[0].header.keys():
		print 'Cannot find RA in %s' % options.metafits
		sys.exit(1)
	if not 'DEC' in f[0].header.keys():
		print 'Cannot find DEC in %s' % options.metafits
		sys.exit(1)
	
	delays=array(map(int,f[0].header['DELAYS'].split(',')))
	normal_delays = map(int,f[0].header['DELAYS'].split(','))
	LST = float(f[0].header['LST'])
	obsID = f[0].header['GPSTIME']
	freqcent = f[0].header['FREQCENT']*1e+6
	ra_point = f[0].header['RA']
	dec_point = f[0].header['DEC']
	
	##Delays have to be in a (2,16) shaped array for the beam functions to work
	if delays.shape == (16,):
		try:
			delays=repeat(reshape(delays,(1,16)),2,axis=0)
		except:
			print 'Unable to convert delays (shape=%s) to (2,16)' % (delays.shape)
	assert delays.shape == (2,16), "Delays %s have unexpected shape %s" % (delays,delays.shape)

	##This method for plotting the contour map is gleamed from MWA_Tools/mwapy/pb/primarybeammap
	contourlevels=[0.01, 0.1, 0.25, 0.5, 0.75]
	
	def scale_ras(in_list,out_list):
		for ra in in_list:
			if ra>180.0:
				ra -= 360.0
			out_list.append(ra)
	
	src_decs = [source.srcdec for source in rts_sources]
	all_decs = [source.srcdec for source in rts_sources]
	
	min_dec = min([min(comp_decs),min(src_decs)])
	max_dec = max([max(comp_decs),max(src_decs)])
	
	comp_scaled = []
	src_scaled = []
	
	scale_ras(comp_ras,comp_scaled)
	scale_ras(src_ras,src_scaled)
	
	min_ra = min([min(comp_scaled),min(src_scaled)])
	max_ra = max([max(comp_scaled),max(src_scaled)])
	
	if min_dec - 10 < -90:
		min_dec_mesh = -90.0
	else:
		min_dec_mesh = min_dec - 10
		
	if max_dec - 10 < -90:
		max_dec_mesh = -90.0
	else:
		max_dec_mesh = max_dec + 10
		
	#if LST-max_ra-10 < -180.0:
		#min_HA_mesh = -180
	#else:
		#min_HA_mesh = LST-max_ra-10
		
	#print LST - min_ra+10
		
	#if LST - min_ra+10 > 360.0:
		#max_HA_mesh = LST - min_ra+10 - 360.0
	#else:
		#max_HA_mesh = LST - min_ra + 10
		
	max_HA_mesh = LST - min_ra+10
	min_HA_mesh = LST - max_ra-10
	
	#print min_HA_mesh,max_HA_mesh,min_dec_mesh,max_dec_mesh
	
	#Set up a HA,Dec axes that covers more than the image
	HA,Dec=meshgrid(arange(min_HA_mesh,max_HA_mesh+0.5,0.5),arange(min_dec_mesh,max_dec_mesh+0.5,0.5))
	
	#Covert to Az,Alt
	mwa=ephem_utils.Obs[ephem_utils.obscode['MWA']]
	mwa_lat = mwa.lat
	Az,Alt=ephem_utils.eq2horz(HA,Dec,mwa_lat)

	za=(90-Alt)*pi/180
	az=Az*pi/180

	d = mwa_tile.Dipole(type='lookup')
	tile = mwa_tile.ApertureArray(dipoles=[d]*16)

	j = tile.getResponse(az,za,freqcent,delays=delays)
	##Convert that in to XX,YY responses
	vis = mwa_tile.makeUnpolInstrumentalResponse(j,j)
	##This is the power in XX,YY - this is taken from primary_beam.MWA_Tile_advanced - prints out
	##lots of debugging messages so have pulled it out of the function
	XX,YY = vis[:,:,0,0].real,vis[:,:,1,1].real
	r = sqrt(XX**2+YY**2)
	
	# normalize
	r/=nanmax(r)

	if (r is None):
		print 'Oh noes!'
	Z2=where(r>=min(contourlevels),r,0)
	
	CS = ax.contour(LST - HA,Dec,Z2,contourlevels,colors='k',transform=tr_fk5,alpha=0.7)
	ax.clabel(CS, inline=1, fontsize=16)
	
	
	RA_up = LST - max_HA_mesh
	RA_down = LST - min_HA_mesh
	
	mid_dec = (max_dec - min_dec) / 2.0
	
	if max_ra - min_ra < 200:
		ax.set_xlim(wcs.wcs_world2pix(max_ra + 5,DEC,0)[0],wcs.wcs_world2pix(min_ra - 5,DEC,0)[0])
	else:
		ax.set_xlim(wcs.wcs_world2pix(max_ra + 10,DEC,0)[0],wcs.wcs_world2pix(min_ra - 10,DEC,0)[0])
		
	if min_dec_mesh<-88:
		ax.set_ylim(wcs.wcs_world2pix(RA,min_dec_mesh,0)[1],wcs.wcs_world2pix(RA,max_dec_mesh-7.0,0)[1])
	else:
		ax.set_ylim(wcs.wcs_world2pix(RA,min_dec_mesh,0)[1],wcs.wcs_world2pix(RA,max_dec_mesh,0)[1])
	
ra_ax = ax.coords[0]
dec_ax = ax.coords[1]
ra_ax.set_axislabel('RAJ2000')
dec_ax.set_axislabel('DECJ2000')
ra_ax.set_major_formatter('hh')#,font='Computer Modern Typewriter')
dec_ax.set_major_formatter('dd:mm:ss')

ra_ax.set_ticks(spacing=900 * units.arcmin)
dec_ax.set_ticks(spacing=900 * units.arcmin)

g = ax.coords.grid(linewidth=1.0,alpha=1.0)

ax.legend(loc='best')

if options.write: fig.savefig("%s.png" %options.srclist.split('.txt')[0],bbox_edges='tight')
else: plt.show()