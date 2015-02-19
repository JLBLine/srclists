#!/usr/bin/python
'''Script to create a single mega-patch calibrator based on the primary beam
convolved with the source fluxes'''
try:
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits
from mwapy import ephem_utils
from mwapy.pb import primary_beam
from numpy import *
from mwapy.pb import mwa_tile
import sys
from optparse import OptionParser,OptionGroup
import matplotlib.pyplot as plt
import numpy as n


##TO DO: 
## - Add precess options?

parser = OptionParser()
parser.add_option('-m','--metafits', 
	help="metafits file to get obs info from")
parser.add_option('-s','--srclist',
	help="Base srclist to get source info from")
parser.add_option('-n','--num_sources',
	help="Number of sources to put in the mage patch")
parser.add_option('-p', '--plot',action='store_true', 
	help='Plot the sources included - NEEDS THE MODULE wcsaxes')
parser.add_option('-c', '--cutoff', default=20.0,
	help='Distance from the pointing centre within which to accept source (cutoff in deg). Default is 20deg')
parser.add_option('-o', '--order', default='flux',
	help='Criteria with which to order the output sources - "flux" for brightest first, "distance" for closest to pointing centre first, "experimental" for a combination. Default = "flux". To use default experimental, "experimental". Other, enter "experimental=flux,distance" with flux cutoff in Jy and distance cutoff in deg. ')

(options, args) = parser.parse_args()

cutoff = float(options.cutoff)

def hour_to_deg(time):      #converts hh:mm:ss.ss in to degrees, must input as a string
	negtest=time[0]
	time=time.split(':')
	degr=float(time[0])*15.0
	mins=float(time[1])*(15.0/60.0)
	secs=float(time[2])*(15.0/(3600.0))
	if negtest=='-':
		deg=degr-mins-secs
	if negtest!='-':
		deg=degr+mins+secs
	return deg

##Try opening the metafits file and complain if certain data is missing
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
	
##Gather the useful info
delays=array(map(int,f[0].header['DELAYS'].split(',')))
LST = float(f[0].header['LST'])
obsID = f[0].header['GPSTIME']
freqcent = f[0].header['FREQCENT']*1e+6
ra_point = f[0].header['RA']
dec_point = f[0].header['DEC']

##Class to store source information with - set to lists
##to store component info in the same place
class rts_source():
	def __init__(self):
		self.name = ''
		self.ras = []
		self.has = []
		self.decs = []
		self.freqs = []
		self.fluxs = []
		self.extrap_fluxs = []
		self.weighted_flux = -1
		self.offset = -1
		
def extrap_flux(freqs,fluxs,extrap_freq):
	'''f1/f2 = (nu1/n2)**alpha
	   alpha = ln(f1/f2) / ln(nu1/nu2)
	   f1 = f2*(nu1/nu2)**alpha'''
	alpha = log(fluxs[0]/fluxs[1]) / log(freqs[0]/freqs[1])
	extrap_flux = fluxs[0]*(extrap_freq/freqs[0])**alpha
	return extrap_flux
	
def arcdist(RA1,RA2,Dec1,Dec2):
	'''calculates distance between two celestial points in degrees'''
	dr = n.pi/180.0
	in1 = (90.0 - Dec1)*dr
	in2 = (90.0 - Dec2)*dr
	RA_d = (RA1 - RA2)*dr
	cosalpha = n.cos(in1)*n.cos(in2) + n.sin(in1)*n.sin(in2)*n.cos(RA_d)
	alpha = n.arccos(cosalpha)
	return alpha/dr

##Find the mwa latitde for converting ha,dec to Az,Alt
mwa=ephem_utils.Obs[ephem_utils.obscode['MWA']]
mwa_lat = mwa.lat

##Lookup an mwa_title (I don't really know precisely what it's doing)
d = mwa_tile.Dipole(type='lookup')
tile = mwa_tile.ApertureArray(dipoles=[d]*16)

#delays = array([0]*16)
#delays_init = array([0]*16)
#delays = array([6,4,2,0,6,4,2,0,6,4,2,0,6,4,2,0])
#delays_init = array([6,4,2,0,6,4,2,0,6,4,2,0,6,4,2,0])

##Delays have to be in a (2,16) shaped array for the beam functions to work
if delays.shape == (16,):
	try:
		delays=repeat(reshape(delays,(1,16)),2,axis=0)
	except:
		print 'Unable to convert delays (shape=%s) to (2,16)' % (delays.shape)
assert delays.shape == (2,16), "Delays %s have unexpected shape %s" % (delays,delays.shape)

sources = []

##Read in the srclist information
rts_srcs = open(options.srclist,'r').read().split('ENDSOURCE')
del rts_srcs[-1]

##Go through all sources in the source list, gather their information, extrapolate
##the flux to the central frequency and weight by the beam at that position
for split_source in rts_srcs:
	
	source = rts_source()
	
	##Find the primary source info - even if no comps, this will isolate
	##primary source infomation
	
	primary_info = split_source.split('COMPONENT')[0].split('\n')
	#print primary_info
	primary_info = [info for info in primary_info if info!='']
	#print primary_info[0]
	meh,prim_name,prim_ra,prim_dec = primary_info[0].split()
	##IF LOOP HERE TO DO DISTANCE CUTOFF========================================
	##==========================================================================
	
	offset = arcdist(float(ra_point),float(prim_ra)*15.0,float(dec_point),float(prim_dec))
	if offset <= cutoff:
		
		##Put in to the source class
		source.name = prim_name
		source.ras.append(float(prim_ra))
		source.decs.append(float(prim_dec))
		source.offset = offset
		##Find the fluxes and append to the source class
		prim_freqs = []
		prim_fluxs = []
		for line in primary_info:
			if 'FREQ' in line:
				prim_freqs.append(float(line.split()[1]))
				prim_fluxs.append(float(line.split()[2]))
		source.freqs.append(prim_freqs)
		source.fluxs.append(prim_fluxs)
		
		##Split all info into lines and get rid of blank entries
		lines = split_source.split('\n')
		lines = [line for line in lines if line!='']
		##If the are components to the source, see where the components start and end
		comp_starts = [lines.index(line) for line in lines if 'COMPONENT' in line and 'END' not in line]
		comp_ends = [i for i in xrange(len(lines)) if lines[i]=='ENDCOMPONENT']
		
		##For each component, go through and find ra,dec,freqs and fluxs
		for start,end in zip(comp_starts,comp_ends):
			freqs = []
			fluxs = []
			for line in lines[start:end]:
				if 'COMPONENT' in line:
					source.ras.append(float(line.split()[1]))
					source.decs.append(float(line.split()[2]))
				elif 'FREQ' in line:
					freqs.append(float(line.split()[1]))
					fluxs.append(float(line.split()[2]))
			source.fluxs.append(fluxs)
			source.freqs.append(freqs)
		
		##For each set of source infomation, calculate and extrapolated flux at the centra flux value
		for freqs,fluxs in zip(source.freqs,source.fluxs):
			##If extrapolating below known freqs, choose two lowest frequencies
			if min(freqs)>freqcent:
				ext_flux = extrap_flux([freqs[0],freqs[1]],[fluxs[0],fluxs[1]],freqcent)
			##If extrapolating above known freqs, choose two highest frequencies
			elif max(freqs)<freqcent:
				ext_flux = extrap_flux([freqs[-2],freqs[-1]],[fluxs[-2],fluxs[-1]],freqcent)
			##Otherwise, choose the two frequencies above and below, and extrap between them
			else:
				for i in xrange(len(freqs)-1):
					if freqs[i]<freqcent and freqs[i+1]>freqcent:
						ext_flux = extrap_flux([freqs[i],freqs[i+1]],[fluxs[i],fluxs[i+1]],freqcent)
			source.extrap_fluxs.append(ext_flux)
			
		beam_weights = []
		
		for ra,dec in zip(source.ras,source.decs):
			##HA=LST-RA in def of ephem_utils.py
			ha = LST - ra*15.0  ##RTS stores things in hours
			##Convert to zenith angle, azmuth in rad
			
			Az,Alt=ephem_utils.eq2horz(ha,dec,mwa_lat)
			za=(90-Alt)*pi/180
			az=Az*pi/180
			
			##Get the tile response for the given sky position,frequency and delay
			##Needs za,az in 2D arrays ([[]] means (1,1)) for za,az
			j = tile.getResponse(array([[az]]),array([[za]]),freqcent,delays=delays)
			##Convert that in to XX,YY responses
			vis = mwa_tile.makeUnpolInstrumentalResponse(j,j)
			##This is the power in XX,YY - this is taken from primary_beam.MWA_Tile_advanced - prints out
			##lots of debugging messages so have pulled it out of the function
			XX,YY = vis[:,:,0,0].real,vis[:,:,1,1].real
			overall_power = n.sqrt(XX[0]**2+YY[0]**2)   ###CHECK THIS - go for rms as this is how Stokes I is made
			beam_weights.append(overall_power[0])
		
		source.beam_weights = beam_weights
		##Dot the weights and the extra fluxes together to come up with a weighted sum 
		##of all components
		source.weighted_flux = dot(array(beam_weights),array(source.extrap_fluxs))
		sources.append(source)
	else:
		pass
	
##Make a list of all of the weighted_fluxes and then order the sources according to those
all_weighted_fluxs = [source.weighted_flux for source in sources]
weighted_sources = [source for flux,source in sorted(zip(all_weighted_fluxs,sources),key=lambda pair: pair[0],reverse=True)][:int(options.num_sources)]

if options.order=='flux':
	ordered_sources = weighted_sources
elif options.order=='distance':
	ordered_offsets = [source.offset for source in weighted_sources]
	ordered_sources = [source for offset,source in sorted(zip(ordered_offsets,weighted_sources),key=lambda pair: pair[0])]
	
elif 'experimental' in options.order:
	if len(options.order.split('='))==1:
		flux_cut,dist_cut = 10.0,1.0
	else:
		flux_cut,dist_cut = options.order.split('=')[1].split(',')
	
	close_fluxs = []
	close_dists = []
	dist_cut = float(dist_cut)
	
	##Try to find all sources within distance cutoff above flux threshold - if none exist, extend search
	##radii by 0.5 deg
	while len(close_fluxs)==0:
		for flux,offset in zip([src.weighted_flux for src in weighted_sources],[src.offset for src in weighted_sources]):
			if flux>float(flux_cut) and offset<dist_cut:
				close_fluxs.append(flux)
				close_dists.append(offset)
		dist_cut+=0.5

	##This is the brightest source within the base source cutoff distance
	brightest_close_flux = sorted(close_fluxs)[0]
	
	##This is where the bright close source appears in the beam weighted source list
	weighted_fluxes = [source.weighted_flux for source in weighted_sources]
	brightest_ind = weighted_fluxes.index(brightest_close_flux)
	
	##Use the positional offset as an identifier as well in case so other source has the
	##same flux
	weighted_offsets = [source.offset for source in weighted_sources]
	brightest_close_offset = weighted_offsets[brightest_ind]
	
	print "Base source convolved flux is %.3fJy at a distance of %.3fdeg from point centre" %(brightest_close_flux,brightest_close_offset)
	
	##Put this source at the top of the ordered list, and then append all other sources after
	##NOTE - this means that apart from the top source, all other sources are flux ordered.
	ordered_sources = [weighted_sources[brightest_ind]]
	for source in weighted_sources:
		if source.weighted_flux!=brightest_close_flux and source.offset!=brightest_close_offset:
			ordered_sources.append(source)
	
##Make a new single patch source based on the user specified number of components
output_name = "%s_%s_cal%s.txt" %(options.srclist.split('/')[-1].split('.')[0],obsID,options.num_sources)
out_file = open(output_name,'w+')

##Print out the strongest source as the primary calibator
for source in ordered_sources[:1]:
	out_file.write('SOURCE %s %.10f %.10f' %(obsID,source.ras[0],source.decs[0]))
	for flux,freq in zip(source.fluxs[0],source.freqs[0]):
		out_file.write("\nFREQ %.4e %.5f 0 0 0" %(freq,flux))
	##Cycle through any components in that primary calibator
	for i in range(1,len(source.ras)):
		out_file.write('\nCOMPONENT %.10f %.10f' %(source.ras[i],source.decs[i]))
		for flux,freq in zip(source.fluxs[i],source.freqs[i]):
			out_file.write("\nFREQ %.4e %.5f 0 0 0" %(freq,flux))
		out_file.write('\nENDCOMPONENT')
		
##For all other sources, add all information as COMPONENTS
for source in ordered_sources[1:int(options.num_sources)]:
	for i in xrange(len(source.ras)):
		out_file.write('\nCOMPONENT %.10f %.10f' %(source.ras[i],source.decs[i]))
		for flux,freq in zip(source.fluxs[i],source.freqs[i]):
			out_file.write("\nFREQ %.4e %.5f 0 0 0" %(freq,flux))
		out_file.write('\nENDCOMPONENT')
		
out_file.write('\nENDSOURCE')
out_file.close()

print "Created %s" %output_name
##Finito!!

if options.plot:
	from astropy.wcs import WCS
	from wcsaxes import WCSAxes
	from astropy import units as units


	##A fits image header with which to create a wcs with
	header = { 'NAXIS'  : 2,             ##Number of data axis
	'NAXIS1' : 10,                  ##Length of X axis
	'CTYPE1' : 'RA---SIN',           ##Projection type of X axis
	'CRVAL1' : 0.0,        ##Central X world coord value
	'CRPIX1' : 5,                    ##Central X Pixel value
	'CUNIT1' : 'deg',                ##Unit of X axes
	'CDELT1' : -1*cos(0.0*(pi/180.0)),              ##Size of pixel in world co-ord
	'NAXIS2' : 10,                  ##Length of X axis
	'CTYPE2' : 'DEC--SIN',           ##Projection along Y axis
	'CRVAL2' : -27.0,                   ##Central Y world coord value
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

	for source in ordered_sources:
		#print source.ras[0],source.decs[0],source.weighted_flux
		ra = source.ras[0]*15.0
		#if ra>180.0: ra-=360.0
		ax.scatter(ra,source.decs[0],marker='o',transform=tr_fk5)

	delays_init = array(map(int,f[0].header['DELAYS'].split(',')))
	title = "delays = %s\nnum_sources = %s" %(str(delays_init),options.num_sources)
	ax.set_title(title)
		
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
