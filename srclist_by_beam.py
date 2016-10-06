#!/usr/bin/python
'''Script to create a single mega-patch calibrator based on the primary beam
convolved with the source fluxes'''
try:
    import astropy.io.fits as pyfits
except ImportError:
    import pyfits
from mwapy import ephem_utils
from mwapy.pb import primary_beam,beam_full_EE
from numpy import *
from mwapy.pb import mwa_tile
from mwapy import pb
import sys
from optparse import OptionParser,OptionGroup
import matplotlib.pyplot as plt
import numpy as n
import subprocess
##TO DO: 
## - Add precess options?

parser = OptionParser()
parser.add_option('-x','--no_patch', action='store_true',
	help="Switch on to not put the sources in to a patch")
parser.add_option('-m','--metafits', 
	help="metafits file to get obs info from")
parser.add_option('-s','--srclist',
	help="Base srclist to get source info from")
parser.add_option('-n','--num_sources',
	help="Number of sources to put in the mega patch")
parser.add_option('-p', '--plot',action='store_true', 
	help='Plot the sources included - NEEDS THE MODULE wcsaxes')
parser.add_option('-c', '--cutoff', default=20.0,
	help='Distance from the pointing centre within which to accept source (cutoff in deg). Default is 20deg')
parser.add_option('-o', '--order', default='flux',
	help='Criteria with which to order the output sources - "flux" for brightest first, "distance" for closest to pointing centre first, "experimental" for a combination, "name=*sourcename*" to force a calibrator. Default = "flux". {To use default experimental, "experimental". Other, enter "experimental=flux,distance" with flux cutoff in Jy and distance cutoff in deg.} ')
parser.add_option('-a', '--outside',action='store_true', default=False,
	help='Switch on to only consider sources OUTSIDE of the cutoff, rather than inside')

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
		self.shapelets = []
		self.shapelet_indexes = []
		self.shapelet_coeffs = []
		self.gaussians = []
		self.gaussian_indexes = []
		self.beam_inds = []
		
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
#d = mwa_tile.Dipole(type='lookup')
#tile = mwa_tile.ApertureArray(dipoles=[d]*16)

delays=repeat(reshape(delays,(1,16)),2,axis=0)
tile=beam_full_EE.ApertureArray(pb.h5file,freqcent)
mybeam=beam_full_EE.Beam(tile, delays)

##Read in the srclist information
rts_srcs = open(options.srclist,'r').read().split('ENDSOURCE')
del rts_srcs[-1]

def create_source(prim_name=None, prim_ra=None, prim_dec=None, offset=None, primary_info=None,beam_ind=None):
	'''Takes the information for a source and puts it into an rts_source class for later use'''
	
	source = rts_source()
	source.name = prim_name
	source.ras.append(float(prim_ra))
	source.decs.append(float(prim_dec))
	all_ras.append(float(prim_ra))
	all_decs.append(float(prim_dec))
	source.beam_inds.append(beam_ind)
	beam_ind += 1
	
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
	##If there are components to the source, see where the components start and end
	comp_starts = [lines.index(line) for line in lines if 'COMPONENT' in line and 'END' not in line]
	comp_ends = [i for i in xrange(len(lines)) if lines[i]=='ENDCOMPONENT']
	
	##Check to see if the primary source is a gaussian or shapelet
	for line in primary_info:
		###Check here to see if primary source is a gaussian:
		if 'GAUSSIAN' in line:
			source.gaussians.append(line)
			source.gaussian_indexes.append(0)
		##As shapelet line comes after the 
		elif 'SHAPELET' in line:
			coeffs = []
			source.shapelets.append(line)
			source.shapelet_indexes.append(0)
			##If the primary source is a shapelet, search for shapelet coeffs in primary data,
			##gather in a list and append to the source class
			for line in primary_info:
				if 'COEFF' in line: coeffs.append(line)
			source.shapelet_coeffs.append(coeffs)
		
	##For each component, go through and find ra,dec,freqs and fluxs
	##Also check here if the component is a gaussian or shapelet
	for start,end in zip(comp_starts,comp_ends):
		freqs = []
		fluxs = []
		coeffs = []
		for line in lines[start:end]:
			if 'COMPONENT' in line:
				source.ras.append(float(line.split()[1]))
				source.decs.append(float(line.split()[2]))
				all_ras.append(float(line.split()[1]))
				all_decs.append(float(line.split()[2]))
				source.beam_inds.append(beam_ind)
				beam_ind += 1
				
			elif 'FREQ' in line:
				freqs.append(float(line.split()[1]))
				fluxs.append(float(line.split()[2]))
			elif 'GAUSSIAN' in line:
				source.gaussians.append(line)
				gaus_ind = comp_starts.index(start) + 1
				source.gaussian_indexes.append(gaus_ind)
			elif 'SHAPELET' in line:
				source.shapelets.append(line)
				shap_ind = comp_starts.index(start) + 1
				source.shapelet_indexes.append(shap_ind)
			elif 'COEFF' in line:
				coeffs.append(line)
		source.fluxs.append(fluxs)
		source.freqs.append(freqs)
		if len(coeffs) > 0:
			source.shapelet_coeffs.append(coeffs)
			
	##For each set of source infomation, calculate and extrapolated flux at the centra flux value
	for freqs,fluxs in zip(source.freqs,source.fluxs):

		##If only one freq, extrap with an SI of -0.8:
		if len(freqs)==1:
			##f1 = c*v1**-0.8
			##c = f1 / (v1**-0.8)
			c = fluxs[0] / (freqs[0]**-0.8)
			ext_flux = c*freqcent**-0.8
		##If extrapolating below known freqs, choose two lowest frequencies
		elif min(freqs)>freqcent:
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
		
	sources.append(source)
		
def get_beam_weights(ras=None,decs=None):
	'''Takes ra and dec coords, and works out the overall beam power
	at that location using the 2016 spherical harmonic beam code from mwapy'''
	
	##For each component, work out it's position, convolve with the beam and sum for the source
	#for ra,dec in zip(source.ras,source.decs):
	##HA=LST-RA in def of ephem_utils.py
	has = LST - array(ras)*15.0  ##RTS stores things in hours
	##Convert to zenith angle, azmuth in rad
	
	Az,Alt=ephem_utils.eq2horz(has,array(decs),mwa_lat)
	za=(90-Alt)*pi/180
	az=Az*pi/180
	
	##Get the tile response for the given sky position,frequency and delay
	##Add this to interpolate over angular response of the beam
	#j = mybeam.get_interp_response(array([az]),array([za]), 5)
	j = mybeam.get_response(array([az]),array([za]))

	j = tile.apply_zenith_norm_Jones(j)
	
	if len(j.shape)==4:
		j=n.swapaxes(n.swapaxes(j,0,2),1,3)
	elif len(j.shape)==3: #1-D
		j=n.swapaxes(n.swapaxes(j,1,2),0,1)
	else: #single value
		pass

#	print 'JSHAPE', j.shape
	
	##Convert that in to XX,YY responses
	vis = mwa_tile.makeUnpolInstrumentalResponse(j,j)
	##This is the power in XX,YY - this is taken from primary_beam.MWA_Tile_advanced - prints out
	##lots of debugging messages so have pulled it out of the function
	XX,YY = vis[:,:,0,0].real,vis[:,:,1,1].real
	
	beam_weights = (XX[0]+YY[0]) / 2.0
	##OLd way of combining XX and YY - end up with beam values greater than 1, not good!
	#beam_weights = n.sqrt(XX[0]**2+YY[0]**2) 
	
	return beam_weights

sources = []
all_ras = []
all_decs = []

beam_ind = 0
##Go through all sources in the source list, gather their information, extrapolate
##the flux to the central frequency and weight by the beam at that position
for split_source in rts_srcs:
	
	source = rts_source()
	
	##Find the primary source info - even if no comps, this will isolate
	##primary source infomation
	primary_info = split_source.split('COMPONENT')[0].split('\n')
	primary_info = [info for info in primary_info if info!='']
	meh,prim_name,prim_ra,prim_dec = primary_info[0].split()
	
	##Check if the primary RA,Dec is below the horizon - it will crash the RTS otherwise
	##Skip if so
	ha_prim = LST - float(prim_ra)*15.0
	Az_prim,Alt_prim = ephem_utils.eq2horz(ha_prim,float(prim_dec),mwa_lat)
	
	if Alt_prim < 0.0:
		pass
	else:
		offset = arcdist(float(ra_point),float(prim_ra)*15.0,float(dec_point),float(prim_dec))
		
		if options.outside:
			if offset > cutoff:
				create_source(prim_name=prim_name, prim_ra=prim_ra, prim_dec=prim_dec, offset=offset, primary_info=primary_info,beam_ind=beam_ind)
				beam_ind += 1
			else:
				pass
		else:
			if offset <= cutoff:
				create_source(prim_name=prim_name, prim_ra=prim_ra, prim_dec=prim_dec, offset=offset, primary_info=primary_info,beam_ind=beam_ind)
				beam_ind += 1
			else:
				pass
			
##Need to work out all beam weightings in one single calculation,
##as in each instance of the beam ~40s to run
beam_weights = get_beam_weights(ras=all_ras,decs=all_decs)

##Go through all the sources, and apply the beam weights to all
##components in the source. Dot the weights and fluxes to get 
##a total weighted flux
for source in sources:
	source_weights = beam_weights[source.beam_inds]
	source.weighted_flux = dot(array(source_weights),array(source.extrap_fluxs))
	
##Make a list of all of the weighted_fluxes and then order the sources according to those
all_weighted_fluxs = [source.weighted_flux for source in sources]
weighted_sources = [source for flux,source in sorted(zip(all_weighted_fluxs,sources),key=lambda pair: pair[0],reverse=True)][:int(options.num_sources)]

if options.no_patch:
	print "++++++++++++++++++++++++++++++++++++++\nCreating weighted srclist - not mega-patching the sources"
	output_name = "%s_%s_peel%s.txt" %(options.srclist.split('/')[-1].split('.')[0],obsID,options.num_sources)
	out_file = open(output_name,'w+')
	for source in weighted_sources[:1]:
		out_file.write('SOURCE %s %.10f %.10f' %(source.name,source.ras[0],source.decs[0]))
		##If base source was a gaussian put it in:
		if len(source.gaussians) > 0 and 0 in source.gaussian_indexes:
			out_file.write('\n'+source.gaussians[0])
		##write out the fluxes
		for flux,freq in zip(source.fluxs[0],source.freqs[0]):
			out_file.write("\nFREQ %.4e %.5f 0 0 0" %(freq,flux))
		##Shapelet coeffs come after the frequencies
		##If base source was a shapelet, put in and it's coefficients in
		#print source.shapelets,source.shapelet_indexes,source.shapelet_coeffs
		if len(source.shapelets) > 0 and 0 in source.shapelet_indexes:
			#print 'here',source.shapelets[0]
			out_file.write('\n'+source.shapelets[0])
			for coeff in source.shapelet_coeffs[0]:
				out_file.write('\n'+coeff)
		##Cycle through any components in that primary calibator
		for i in range(1,len(source.ras)):
			out_file.write('\nCOMPONENT %.10f %.10f' %(source.ras[i],source.decs[i]))
			##Cycle through gaussian IDs and add the gaussian line if applicable
			for gaus_ind,gaus_line in zip(source.gaussian_indexes,source.gaussians):
				if gaus_ind == i: out_file.write('\n'+gaus_line)
			##add the fluxes
			for flux,freq in zip(source.fluxs[i],source.freqs[i]):
				out_file.write("\nFREQ %.4e %.5f 0 0 0" %(freq,flux))
			##Have a look and see if there is a shapelet in this source
			##and write out if applicable
			try:
				shap_ind = source.shapelet_indexes.index(i)
				out_file.write('\n'+source.shapelets[shap_ind])
				for line in source.shapelet_coeffs[shap_ind]:
					out_file.write('\n'+line)
			except:
				pass
			out_file.write('\nENDCOMPONENT')
	out_file.write('\nENDSOURCE')
	
	for source in weighted_sources[1:int(options.num_sources)]:
		out_file.write('\nSOURCE %s %.10f %.10f' %(source.name,source.ras[0],source.decs[0]))
		##If base source was a gaussian put it in:
		if len(source.gaussians) > 0 and 0 in source.gaussian_indexes:
			out_file.write('\n'+source.gaussians[0])
		##write out the fluxes
		for flux,freq in zip(source.fluxs[0],source.freqs[0]):
			out_file.write("\nFREQ %.4e %.5f 0 0 0" %(freq,flux))
		##Shapelet coeffs come after the frequencies
		##If base source was a shapelet, put in and it's coefficients in
		if len(source.shapelets) > 0 and 0 in source.shapelet_indexes:
			out_file.write('\n'+source.shapelets[0])	
			for coeff in source.shapelet_coeffs[0]:
				out_file.write('\n'+coeff)
		##Cycle through any components in that primary calibator
		for i in range(1,len(source.ras)):
			out_file.write('\nCOMPONENT %.10f %.10f' %(source.ras[i],source.decs[i]))
			##Cycle through gaussian IDs and add the gaussian line if applicable
			for gaus_ind,gaus_line in zip(source.gaussian_indexes,source.gaussians):
				if gaus_ind == i: out_file.write('\n'+gaus_line)
			##add the fluxes
			for flux,freq in zip(source.fluxs[i],source.freqs[i]):
				out_file.write("\nFREQ %.4e %.5f 0 0 0" %(freq,flux))
			##Have a look and see if there is a shapelet in this source
			##and write out if applicable
			try:
				shap_ind = source.shapelet_indexes.index(i)
				out_file.write('\n'+source.shapelets[shap_ind])
				for line in source.shapelet_coeffs[shap_ind]:
					out_file.write('\n'+line)
			except:
				pass
			out_file.write('\nENDCOMPONENT')
		out_file.write('\nENDSOURCE')
	out_file.close()

else:
	if options.order=='flux':
		ordered_sources = weighted_sources

	elif options.order=='distance':
		ordered_offsets = [source.offset for source in weighted_sources]
		ordered_sources = [source for offset,source in sorted(zip(ordered_offsets,weighted_sources),key=lambda pair: pair[0])]
		
		
	elif 'name' in options.order:
		name = options.order.split("=")[1]
		top_source_ind = [source.name for source in weighted_sources].index(name)
		top_source = weighted_sources[top_source_ind]
		ordered_sources = [top_source]
		for i in xrange(len(weighted_sources)):
			if i!=top_source_ind:
				ordered_sources.append(weighted_sources[i])
		print '++++++++++++++++++++++++++++++++++++++\nBase Source forced as %s with \nconvolved flux %.1fJy at a distance %.2fdeg\n---------------------------------' %(top_source.name,top_source.weighted_flux,top_source.offset)
		
	elif 'experimental' in options.order:
		#print 'here, here, here'
		
		
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
			if dist_cut>cutoff:
				print "++++++++++++++++++++++++++++++++++++++\nNo source above %.2fJy within initial cutoff distance\nNO SOURCE LIST GENERATED\n++++++++++++++++++++++++++++++++++++++" %float(flux_cut)
				sys.exit()

		##This is the brightest source within the base source cutoff distance
		brightest_close_flux = sorted(close_fluxs,reverse=True)[0]
		
		##This is where the bright close source appears in the beam weighted source list
		weighted_fluxes = [source.weighted_flux for source in weighted_sources]
		brightest_ind = weighted_fluxes.index(brightest_close_flux)
		
		##Use the positional offset as an identifier as well in case so other source has the
		##same flux
		weighted_offsets = [source.offset for source in weighted_sources]
		brightest_close_offset = weighted_offsets[brightest_ind]
		
		##Find out name of source
		weighted_names = [source.name for source in weighted_sources]
		brightest_close_name = weighted_names[brightest_ind]
		
		print "++++++++++++++++++++++++++++++++++++++\nBase source %s convolved flux is %.3fJy at a distance \nof %.3fdeg from point centre\n---------------------------------" %(brightest_close_name,brightest_close_flux,brightest_close_offset)
		
		##Put this source at the top of the ordered list, and then append all other sources after
		##NOTE - this means that apart from the top source, all other sources are flux ordered.
		ordered_sources = [weighted_sources[brightest_ind]]
		for source in weighted_sources:
			if source.weighted_flux!=brightest_close_flux and source.offset!=brightest_close_offset:
				ordered_sources.append(source)
		
	##Make a new single patch source based on the user specified number of components
	output_name = "%s_%s_patch%s.txt" %(options.srclist.split('/')[-1].split('.')[0],obsID,options.num_sources)
	out_file = open(output_name,'w+')

	##Print out the strongest source as the primary calibator
	for source in ordered_sources[:1]:
		out_file.write('SOURCE %s %.10f %.10f' %(obsID,source.ras[0],source.decs[0]))
		##If base source was a gaussian put it in:
		if len(source.gaussians) > 0 and 0 in source.gaussian_indexes:
			out_file.write('\n'+source.gaussians[0])
		##write out the fluxes
		for flux,freq in zip(source.fluxs[0],source.freqs[0]):
			out_file.write("\nFREQ %.4e %.5f 0 0 0" %(freq,flux))
		##Shapelet coeffs come after the frequencies
		##If base source was a shapelet, put in and it's coefficients in
		if len(source.shapelets) > 0 and 0 in source.shapelet_indexes:
			out_file.write('\n'+source.shapelets[0])	
			for coeff in source.shapelet_coeffs[0]:
				out_file.write('\n'+coeff)
		##Cycle through any components in that primary calibator
		for i in range(1,len(source.ras)):
			out_file.write('\nCOMPONENT %.10f %.10f' %(source.ras[i],source.decs[i]))
			##Cycle through gaussian IDs and add the gaussian line if applicable
			for gaus_ind,gaus_line in zip(source.gaussian_indexes,source.gaussians):
				if gaus_ind == i: out_file.write('\n'+gaus_line)
			##add the fluxes
			for flux,freq in zip(source.fluxs[i],source.freqs[i]):
				out_file.write("\nFREQ %.4e %.5f 0 0 0" %(freq,flux))
			##Have a look and see if there is a shapelet in this source
			##and write out if applicable
			try:
				shap_ind = source.shapelet_indexes.index(i)
				out_file.write('\n'+source.shapelets[shap_ind])
				for line in source.shapelet_coeffs[shap_ind]:
					out_file.write('\n'+line)
			except:
				pass
			out_file.write('\nENDCOMPONENT')

	##For all other sources, add all information as COMPONENTS
	for source in ordered_sources[1:int(options.num_sources)]:
		for i in xrange(len(source.ras)):
			out_file.write('\nCOMPONENT %.10f %.10f' %(source.ras[i],source.decs[i]))
			##Cycle through gaussian IDs and add the gaussian line if applicable
			for gaus_ind,gaus_line in zip(source.gaussian_indexes,source.gaussians):
				if gaus_ind == i: out_file.write('\n'+gaus_line)
			#Add the fluxes
			for flux,freq in zip(source.fluxs[i],source.freqs[i]):
				out_file.write("\nFREQ %.4e %.5f 0 0 0" %(freq,flux))
			##Have a look and see if there is a shapelet in this source
			##and write out if applicable
			try:
				shap_ind = source.shapelet_indexes.index(i)
				out_file.write('\n'+source.shapelets[shap_ind])
				for line in source.shapelet_coeffs[shap_ind]:
					out_file.write('\n'+line)
			except:
				pass
			out_file.write('\nENDCOMPONENT')
	out_file.write('\nENDSOURCE')
	out_file.close()

print "Created %s\n++++++++++++++++++++++++++++++++++++++" %output_name
##Finito!!

if options.plot:
	cmd = "./plot_srclist.py -m %s -s %s" %(options.metafits, output_name)
	subprocess.call(cmd,shell=True)
